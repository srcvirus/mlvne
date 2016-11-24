#include "cplex_solver.h"
#include <unistd.h>

MultiLayerVNESolver::MultiLayerVNESolver(
    Graph *pn_topology, Graph *vn_topology,
    std::vector<std::vector<int> > *location_constraint, int cost_new_ip_link,
    int k) {
  model_ = IloModel(env_);
  cplex_ = IloCplex(model_);
  constraints_ = IloConstraintArray(env_);
  preferences_ = IloNumArray(env_);
  objective_ = IloExpr(env_);

  pn_topology_ = pn_topology;
  vn_topology_ = vn_topology;
  location_constraint_ = location_constraint;
  k_ = k;
  cost_new_ip_link_ = cost_new_ip_link;
  x_mn_uvi_ = IloIntVar5dArray(env_, vn_topology_->node_count());
  y_m_u_ = IloIntVar2dArray(env_, vn_topology_->node_count());
  l_m_u_ = IloInt2dArray(env_, vn_topology_->node_count());

  // Decision variable initialization for virtual network.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    x_mn_uvi_[m] = IloIntVar4dArray(env_, vn_topology_->node_count());
    for (int n = 0; n < vn_topology_->node_count(); ++n) {
      x_mn_uvi_[m][n] = IloIntVar3dArray(env_, pn_topology_->node_count());
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        x_mn_uvi_[m][n][u] = IloIntVar2dArray(env_, pn_topology_->node_count());
        for (int v = 0; v < pn_topology_->node_count(); ++v) {
          // TODO: We are assuming that there can be at most 20 parallel links
          // in between a pair of nodes, in reality this can be more. Probably
          // need to fix it at some point.
          x_mn_uvi_[m][n][u][v] = IloIntVarArray(env_, 20, 0, 10);
        }
      }
    }
  }
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    y_m_u_[m] = IloIntVarArray(env_, pn_topology_->node_count(), 0, 1);
    l_m_u_[m] = IloIntArray(env_, pn_topology_->node_count(), 0, 1);
  }
  DEBUG("vnodes = %d\n", vn_topology_->node_count());
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < pn_topology_->node_count(); ++u) l_m_u_[m][u] = 0;
    const std::vector<int> &loc_constraints = location_constraint_->at(m);
    for (int i = 0; i < loc_constraints.size(); ++i) {
      DEBUG("loc = %d, pnodes = %d\n", loc_constraints[i],
            pn_topology_->node_count());
      assert(loc_constraints[i] < pn_topology_->node_count());
      l_m_u_[m][loc_constraints[i]] = 1;
    }
  }
}

void MultiLayerVNESolver::BuildModel() {
  // Constraint: Location constraint of virtual nodes.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < pn_topology_->node_count(); ++u) {
      constraints_.add(y_m_u_[m][u] <= l_m_u_[m][u]);
    }
  }

  // Constraint: Capacity constraint of physical links.
  for (int u = 0; u < pn_topology_->node_count(); ++u) {
    const std::vector<edge_endpoint> &u_neighbors =
        pn_topology_->adj_list()->at(u);
    for (int i = 0; i < u_neighbors.size(); ++i) {
      const edge_endpoint &end_point = u_neighbors[i];
      int v = end_point.node_id;
      long beta_uv = end_point.residual_bandwidth;
      assert(beta_uv > 0);
      int order = end_point.order;
      int cost_uv = end_point.cost;
      IloIntExpr sum(env_);
      for (int m = 0; m < vn_topology_->node_count(); ++m) {
        const std::vector<edge_endpoint> &m_neighbors =
            vn_topology_->adj_list()->at(m);
        for (int j = 0; j < m_neighbors.size(); ++j) {
          const edge_endpoint &vend_point = m_neighbors[j];
          int n = vend_point.node_id;
          int beta_mn = vend_point.bandwidth;
          if (!pn_topology_->IsPseudoEdge(u, v, order)) {
            beta_mn *= k_;
          }
          DEBUG("u = %d, v = %d, order = %d, m = %d, n = %d\n", u, v, order, m,
                n);
          sum += (x_mn_uvi_[m][n][u][v][order] + x_mn_uvi_[m][n][v][u][order]) *
                 beta_mn;
        }
      }
      constraints_.add(sum <= beta_uv);
    }
  }

  // Constraint: Every virtual link is mapped to one or more physical links.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const std::vector<edge_endpoint> &m_neighbors =
        vn_topology_->adj_list()->at(m);
    for (int i = 0; i < m_neighbors.size(); ++i) {
      const edge_endpoint &vend_point = m_neighbors[i];
      int n = vend_point.node_id;
      IloIntExpr sum(env_);
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        const std::vector<edge_endpoint> &u_neighbors =
            pn_topology_->adj_list()->at(u);
        for (int j = 0; j < u_neighbors.size(); ++j) {
          const edge_endpoint &end_point = u_neighbors[j];
          int v = end_point.node_id;
          int order = end_point.order;
          constraints_.add(IloIfThen(env_, x_mn_uvi_[m][n][u][v][order] == 1,
                                     x_mn_uvi_[m][n][v][u][order] == 0));
          constraints_.add(IloIfThen(env_, x_mn_uvi_[m][n][v][u][order] == 1,
                                     x_mn_uvi_[m][n][u][v][order] == 0));
          sum += x_mn_uvi_[m][n][u][v][order];
        }
      }
      constraints_.add(sum >= 1);
    }
  }

  // Constraint: Every virtual node is mapped to exactly one physical node.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    IloIntExpr sum(env_);
    for (int u = 0; u < pn_topology_->node_count(); ++u) {
      DEBUG("u = %d, m = %d,\n", u, m);
      sum += y_m_u_[m][u];
    }
    constraints_.add(sum == 1);
  }

  // Constraint: No two virtual nodes are mapped to the same physical node.
  for (int u = 0; u < pn_topology_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int m = 0; m < vn_topology_->node_count(); ++m) {
      DEBUG("u = %d, m = %d,\n", u, m);
      sum += y_m_u_[m][u];
    }
    constraints_.add(sum <= 1);
  }

  // Constraint: Flow constraint to ensure path connectivity.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const std::vector<edge_endpoint> &m_neighbors =
        vn_topology_->adj_list()->at(m);
    for (int i = 0; i < m_neighbors.size(); ++i) {
      const edge_endpoint &vend_point = m_neighbors[i];
      int n = vend_point.node_id;
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        IloIntExpr sum(env_);
        const std::vector<edge_endpoint> &u_neighbors =
            pn_topology_->adj_list()->at(u);
        for (int j = 0; j < u_neighbors.size(); ++j) {
          const edge_endpoint &end_point = u_neighbors[j];
          int v = end_point.node_id;
          int order = end_point.order;
          sum += (x_mn_uvi_[m][n][u][v][order] - x_mn_uvi_[m][n][v][u][order]);
        }
        constraints_.add(sum == (y_m_u_[m][u] - y_m_u_[n][u]));
      }
    }
  }

  // Objective function.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const std::vector<edge_endpoint> &m_neighbors =
        vn_topology_->adj_list()->at(m);
    for (int i = 0; i < m_neighbors.size(); ++i) {
      const edge_endpoint &vend_point = m_neighbors[i];
      int n = vend_point.node_id;
      if (m < n) continue;
      long beta_mn = vend_point.bandwidth;
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        const std::vector<edge_endpoint> &u_neighbors =
            pn_topology_->adj_list()->at(u);
        for (int j = 0; j < u_neighbors.size(); ++j) {
          const edge_endpoint &end_point = u_neighbors[j];
          int v = end_point.node_id;
          int order = end_point.order;
          int cost_uv = end_point.cost;
          if (!pn_topology_->IsPseudoEdge(u, v, order)) {
            objective_ += (x_mn_uvi_[m][n][u][v][order] * cost_new_ip_link_);
          }
          objective_ += (x_mn_uvi_[m][n][u][v][order] * cost_uv * beta_mn);
          DEBUG("u = %d, v = %d, order = %d, m = %d, n = %d\n", u, v, order, m,
                n);
        }
      }
    }
  }
  constraints_.add(objective_ > 0);
  model_.add(constraints_);
  model_.add(IloMinimize(env_, objective_));
}

bool MultiLayerVNESolver::Solve() {
  // TODO(shihab): Tune parameters of CPLEX solver.
  int n_threads = sysconf(_SC_NPROCESSORS_ONLN) * 2;
  if (n_threads < 64) n_threads = 64;
  cplex_.setParam(IloCplex::Threads, n_threads);
  cplex_.exportModel("drone.lp");
  bool is_success = cplex_.solve();
  // return is_success;
  if (cplex_.getStatus() == IloAlgorithm::Infeasible) {
    IloConstraintArray infeasible(env_);
    IloNumArray preferences(env_);
    infeasible.add(constraints_);
    for (int i = 0; i < infeasible.getSize(); ++i) preferences.add(1.0);
    if (cplex_.refineConflict(infeasible, preferences)) {
      IloCplex::ConflictStatusArray conflict = cplex_.getConflict(infeasible);
      env_.getImpl()->useDetailedDisplay(IloTrue);
      std::cout << "Conflict : " << std::endl;
      for (IloInt i = 0; i < infeasible.getSize(); i++) {
        // std::cout << conflict[i] << std::endl;
        if (conflict[i] == IloCplex::ConflictMember)
          std::cout << "Proved  : " << infeasible[i] << std::endl;
        if (conflict[i] == IloCplex::ConflictPossibleMember)
          std::cout << "Possible: " << infeasible[i] << std::endl;
      }
    }
  }
  return is_success;
}

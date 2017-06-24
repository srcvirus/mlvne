#include "cplex_solver.h"
#include <unistd.h>

MultiLayerVNESolver::MultiLayerVNESolver(
    Graph* ip_topology, Graph* otn_topology, Graph* vn_topology,
    OverlayMapping* ip_otn_mapping,
    std::vector<std::vector<int> >* location_constraint) {
  model_ = IloModel(env_);
  cplex_ = IloCplex(model_);
  constraints_ = IloConstraintArray(env_);
  preferences_ = IloNumArray(env_);
  objective_ = IloExpr(env_);
  ip_topology_ = ip_topology;
  otn_topology_ = otn_topology;
  vn_topology_ = vn_topology;
  location_constraint_ = location_constraint;

  // Initialize bandwidth and cost matrix.
  DEBUG("Initializing bw/cost matrix.\n");
  b_uvi_.resize(ip_topology_->node_count());
  cost_uvi_.resize(ip_topology_->node_count());
  for (int i = 0; i < ip_topology_->node_count(); ++i) {
    b_uvi_[i].resize(ip_topology_->node_count());
    cost_uvi_[i].resize(ip_topology_->node_count());
  }
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    const auto& u_neighbors = ip_topology_->adj_list()->at(u);
    for (const auto end_point : u_neighbors) {
      int v = end_point.node_id;
      int order = end_point.order;
      long bandwidth = end_point.bandwidth;
      int cost = end_point.cost;
      if (order >= b_uvi_[u][v].size()) b_uvi_[u][v].resize(order + 1);
      if (order >= cost_uvi_[u][v].size()) cost_uvi_[u][v].resize(order + 1);
      b_uvi_[u][v][order] = bandwidth;
      cost_uvi_[u][v][order] = cost;
    }
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      b_uvi_[u][v].resize(ip_topology_->GetPortCount(u),
                          std::min(ip_topology_->GetPortCapacity(u),
                                   ip_topology_->GetPortCapacity(v)));
      cost_uvi_[u][v].resize(ip_topology_->GetPortCount(u), 1);
      DEBUG("ip_topology_->GetPortCount(u) = %d\n",
            ip_topology_->GetPortCount(u));
    }
  }

  // Decision variables.
  x_mn_uvi_ = IloIntVar5dArray(env_, vn_topology_->node_count());
  y_m_u_ = IloIntVar2dArray(env_, vn_topology_->node_count());
  z_uvi_pq_ = IloIntVar5dArray(env_, ip_topology_->node_count());
  gamma_uvi_ = IloIntVar3dArray(env_, ip_topology_->node_count());

  // Input variables.
  l_m_u_ = IloInt2dArray(env_, vn_topology_->node_count());
  ip_link_uvi_ = IloInt3dArray(env_, ip_topology_->node_count());
  tau_u_p_ = IloInt2dArray(env_, ip_topology_->node_count());

  // Decision variable initialization for virtual network.
  DEBUG("Initializeing variables for VN to IP link embedding.\n");
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    x_mn_uvi_[m] = IloIntVar4dArray(env_, vn_topology_->node_count());
    for (int n = 0; n < vn_topology_->node_count(); ++n) {
      x_mn_uvi_[m][n] = IloIntVar3dArray(env_, ip_topology_->node_count());
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        x_mn_uvi_[m][n][u] = IloIntVar2dArray(env_, ip_topology_->node_count());
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          x_mn_uvi_[m][n][u][v] =
              IloIntVarArray(env_, ip_topology_->GetPortCount(u) + 1, 0, 1);
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            std::string var_name =
                "x[" + std::to_string(m) + "][" + std::to_string(n) + "][" +
                std::to_string(u) + "][" + std::to_string(v) + "][" +
                std::to_string(order) + "]";
            x_mn_uvi_[m][n][u][v][order] =
                IloIntVar(env_, 0, 1, var_name.c_str());
          }
        }
      }
    }
  }
  DEBUG("Initializing variables for VNode to IP node mapping.\n");
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    y_m_u_[m] = IloIntVarArray(env_, ip_topology_->node_count(), 0, 1);
    l_m_u_[m] = IloIntArray(env_, ip_topology_->node_count(), 0, 1);
    for (int u = 0; u < ip_topology_->node_count(); ++u) {
      std::string var_name =
          "y[" + std::to_string(m) + "][" + std::to_string(u) + "]";
      y_m_u_[m][u] = IloIntVar(env_, 0, 1, var_name.c_str());
    }
  }
  DEBUG("vnodes = %d\n", vn_topology_->node_count());
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < ip_topology_->node_count(); ++u) l_m_u_[m][u] = 0;
    const auto& loc_constraints = location_constraint_->at(m);
    for (int i = 0; i < loc_constraints.size(); ++i) {
      DEBUG("loc = %d, pnodes = %d\n", loc_constraints[i],
            ip_topology_->node_count());
      assert(loc_constraints[i] < ip_topology_->node_count());
      l_m_u_[m][loc_constraints[i]] = 1;
    }
  }

  // Initialize decision variables for creating new IP links.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    gamma_uvi_[u] = IloIntVar2dArray(env_, ip_topology_->node_count());
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      gamma_uvi_[u][v] =
          IloIntVarArray(env_, ip_topology_->GetPortCount(u) + 1, 0, 1);
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        std::string var_name = "gamma[" + std::to_string(u) + "][" +
                               std::to_string(v) + "][" +
                               std::to_string(order) + "]";
        gamma_uvi_[u][v][order] = IloIntVar(env_, 0, 1, var_name.c_str());
      }
    }
  }

  // Initialize the input variable that represents existing IP links.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    ip_link_uvi_[u] = IloInt2dArray(env_, ip_topology_->node_count());
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      DEBUG("u = %d, v = %d, port_count(u) = %d\n", u, v,
            ip_topology_->GetPortCount(u));
      ip_link_uvi_[u][v] =
          IloIntArray(env_, ip_topology_->GetPortCount(u) + 1, 0, 1);
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        ip_link_uvi_[u][v][order] = 0;
      }
    }
    const auto& u_neighbors = ip_topology_->adj_list()->at(u);
    for (const auto end_point : u_neighbors) {
      int v = end_point.node_id;
      int order = end_point.order;
      DEBUG("u = %d, v = %d, order = %d\n", u, v, order);
      ip_link_uvi_[u][v][order] = 1;
    }
  }

  // Initialize IP to OTN node attachment input variable.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    tau_u_p_[u] = IloIntArray(env_, otn_topology_->node_count(), 0, 1);
    for (int p = 0; p < otn_topology_->node_count(); ++p) tau_u_p_[u][p] = 0;
    tau_u_p_[u][ip_otn_mapping->node_map[u]] = 1;
  }

  // Initialize IP to OTN mapping decision variable.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    z_uvi_pq_[u] = IloIntVar4dArray(env_, ip_topology_->node_count());
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      z_uvi_pq_[u][v] =
          IloIntVar3dArray(env_, ip_topology_->GetPortCount(u) + 1);
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        z_uvi_pq_[u][v][order] =
            IloIntVar2dArray(env_, otn_topology_->node_count());
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          DEBUG("u = %d, v = %d, order = %d, p = %d\n", u, v, order, p);
          z_uvi_pq_[u][v][order][p] =
              IloIntVarArray(env_, otn_topology_->node_count(), 0, 1);
          for (int q = 0; q < otn_topology_->node_count(); ++q) {
            std::string var_name =
                "z[" + std::to_string(u) + "][" + std::to_string(v) + "][" +
                std::to_string(order) + "][" + std::to_string(p) + "][" +
                std::to_string(q) + "]";
            z_uvi_pq_[u][v][order][p][q] =
                IloIntVar(env_, 0, 1, var_name.c_str());
          }
        }
      }
    }
  }
}

void MultiLayerVNESolver::BuildModel() {
  DEBUG("Building model.\n");
  // Constraint (1): A virtual link can be mapped to an existing or a newly
  // created IP link.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto& m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto vend_point : m_neighbors) {
      int n = vend_point.node_id;
      if (m < n) continue;
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < std::min(ip_topology_->GetPortCount(u),
                                               ip_topology_->GetPortCount(v));
               ++order) {
            constraints_.add(x_mn_uvi_[m][n][u][v][order] <=
                             ip_link_uvi_[u][v][order] +
                                 gamma_uvi_[v][u][order] +
                                 gamma_uvi_[u][v][order]);
          }
        }
      }
    }
  }

  // Constraint (2): Every virtual link is mapped to one or more IP link(s).
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto& m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto vend_point : m_neighbors) {
      int n = vend_point.node_id;
      if (m < n) continue;
      IloIntExpr sum(env_);
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            DEBUG("m = %d, n = %d, u = %d, v = %d, order = %d\n", m, n, u, v,
                  order);
            for (int other_order = 0;
                 other_order < ip_topology_->GetPortCount(v); ++other_order) {
              constraints_.add(
                  IloIfThen(env_, x_mn_uvi_[m][n][u][v][order] == 1,
                            x_mn_uvi_[m][n][v][u][other_order] == 0));
              constraints_.add(
                  IloIfThen(env_, x_mn_uvi_[m][n][v][u][other_order] == 1,
                            x_mn_uvi_[m][n][u][v][order] == 0));
            }
            sum += x_mn_uvi_[m][n][u][v][order];
          }
        }
      }
      constraints_.add(sum >= 1);
    }
  }

  // Constraint (3): Only one instance of an IP link is selected between a pair 
  // of IP nodes.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto& m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto vend_point : m_neighbors) {
      int n = vend_point.node_id;
      if (m < n) continue;
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          IloIntExpr sum(env_);
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            sum += x_mn_uvi_[m][n][u][v][order];
            for (int order_other = 0;
                 order_other < ip_topology_->GetPortCount(v); ++order_other) {
              constraints_.add(gamma_uvi_[u][v][order] + gamma_uvi_[v][u][order_other] <= 1);
              // constraints_.add(IloIfThen(env_, gamma_uvi_[u][v][order] == 1,
              //                            gamma_uvi_[v][u][order_other] == 0));
            }
          }
          constraints_.add(sum <= 1);
        }
      }
    }
  }

  // Constraint (4): Capacity constraint for IP links.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        IloIntExpr sum(env_);
        for (int m = 0; m < vn_topology_->node_count(); ++m) {
          const auto& m_neighbors = vn_topology_->adj_list()->at(m);
          for (const auto vend_point : m_neighbors) {
            int n = vend_point.node_id;
            if (m < n) continue;
            long bandwidth = vend_point.bandwidth;
            sum += (x_mn_uvi_[m][n][u][v][order] * bandwidth);
          }
        }
        // Is this constraint correct? Verify !!
        constraints_.add(sum <= b_uvi_[u][v][order]);
      }
    }
  }

  // Constraint (5): Flow conservation constraint for the IP layer.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto& m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto vend_point : m_neighbors) {
      int n = vend_point.node_id;
      if (m < n) continue;
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        const auto& neighbors = ip_topology_->adj_list()->at(u);
        IloIntExpr sum(env_);
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < std::min(ip_topology_->GetPortCount(u),
                                               ip_topology_->GetPortCount(v));
               ++order) {
            sum +=
                (x_mn_uvi_[m][n][u][v][order] - x_mn_uvi_[m][n][v][u][order]);
          }
        }
        constraints_.add(sum == (y_m_u_[m][u] - y_m_u_[n][u]));
      }
    }
  }

  // Constraint (6): Location constraint of virtual nodes.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < ip_topology_->node_count(); ++u) {
      constraints_.add(y_m_u_[m][u] <= l_m_u_[m][u]);
    }
  }

  // Constraint (7): Every virtual node is mapped to exactly one physical node.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    IloIntExpr sum(env_);
    for (int u = 0; u < ip_topology_->node_count(); ++u) {
      sum += y_m_u_[m][u];
    }
    constraints_.add(sum == 1);
  }

  // Constraint (8): No two virtual nodes are mapped to the same physical node.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int m = 0; m < vn_topology_->node_count(); ++m) {
      sum += y_m_u_[m][u];
    }
    constraints_.add(sum <= 1);
  }

  // Constraint (9): A specific instance of IP link between a pair of IP nodes 
  // is either newly created or alredy existed, but not both.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        constraints_.add(gamma_uvi_[u][v][order] + ip_link_uvi_[u][v][order] <=
                         1);
      }
    }
  }

  // Constraint (10): Number of newly created IP links should be constrainted 
  // by the number of ports on an IP node.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < std::min(ip_topology_->GetPortCount(u),
                                           ip_topology_->GetPortCount(v));
           ++order) {
        sum += (gamma_uvi_[u][v][order] + gamma_uvi_[v][u][order] +
                ip_link_uvi_[u][v][order]);
      }
    }
    constraints_.add(sum <= ip_topology_->GetPortCount(u));
  }

  // Constraint (11): Only the newly created IP links are mapped on the OTN 
  // layer.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          const auto& p_neighbors = otn_topology_->adj_list()->at(p);
          for (const auto& end_point : p_neighbors) {
            int q = end_point.node_id;
            constraints_.add(z_uvi_pq_[u][v][order][p][q] <=
                             gamma_uvi_[u][v][order]);
          }
        }
      }
    }
  }

  // Constraint (13): Capacity constraint for the OTN layer.
  for (int p = 0; p < otn_topology_->node_count(); ++p) {
    const auto& p_neighbors = otn_topology_->adj_list()->at(p);
    for (const auto end_point : p_neighbors) {
      int q = end_point.node_id;
      long bw_pq = end_point.bandwidth;
      IloIntExpr sum(env_);
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            sum += (z_uvi_pq_[u][v][order][p][q] * b_uvi_[u][v][order]);
          }
        }
      }
      constraints_.add(sum <= bw_pq);
    }
  }

  // Constraint (14): Flow conservation constraint for the OTN layer.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          IloIntExpr sum(env_);
          const auto& p_neighbors = otn_topology_->adj_list()->at(p);
          for (int i = 0; i < p_neighbors.size(); ++i) {
            const auto& end_point = p_neighbors[i];
            int q = end_point.node_id;
            sum +=
                (z_uvi_pq_[u][v][order][p][q] - z_uvi_pq_[u][v][order][q][p]);
          }
          if (tau_u_p_[u][p] == 1) {
            constraints_.add(sum == gamma_uvi_[u][v][order]);
          } else if (tau_u_p_[v][p] == 1) {
            constraints_.add(sum == -gamma_uvi_[u][v][order]);
          } else {
            constraints_.add(sum == 0);
          }
        }
      }
    }
  }

  // Objective function.
  // First, add cost corresponding to the number of newly created IP links.
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          const auto& p_neighbors = otn_topology_->adj_list()->at(p);
          for (const auto end_point : p_neighbors) {
            int q = end_point.node_id;
            objective_ += (z_uvi_pq_[u][v][order][p][q] * b_uvi_[u][v][order] *
                           end_point.cost);
          }
        }
      }
    }
  }

  // Add cost component corresponding to virtual to IP mapping.
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto& m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto& vend_point : m_neighbors) {
      int n = vend_point.node_id;
      long bw = vend_point.bandwidth;
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            objective_ +=
                (x_mn_uvi_[m][n][u][v][order] * bw * cost_uvi_[u][v][order]);
          }
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
  int n_threads = sysconf(_SC_NPROCESSORS_ONLN);
  if (n_threads < 64) n_threads = 64;
  cplex_.setParam(IloCplex::Threads, n_threads);
  cplex_.exportModel("mlvne.lp");
  bool success = cplex_.solve();
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
  return success;
}

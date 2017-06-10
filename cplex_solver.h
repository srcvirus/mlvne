#ifndef CPLEX_SOLVER_
#define CPLEX_SOLVER_

#include "datastructure.h"
#include "util.h"

#include <ilcplex/ilocplex.h>

// Type definitions for holding upto 5-dimensional decision variables.
typedef IloArray<IloIntVarArray> IloIntVar2dArray;
typedef IloArray<IloIntVar2dArray> IloIntVar3dArray;
typedef IloArray<IloIntVar3dArray> IloIntVar4dArray;
typedef IloArray<IloIntVar4dArray> IloIntVar5dArray;

typedef IloArray<IloIntArray> IloInt2dArray;
typedef IloArray<IloInt2dArray> IloInt3dArray;

typedef IloArray<IloExprArray> IloExpr2dArray;

using std::string;

class MultiLayerVNESolver {
 public:
  MultiLayerVNESolver() {}
  MultiLayerVNESolver(Graph *ip_topology, Graph *otn_topology,
                      Graph *vn_topology, OverlayMapping *ip_otn_mapping,
                      std::vector<std::vector<int>> *location_constraint);

  IloEnv &env() { return env_; }
  IloModel &model() { return model_; }
  IloCplex &cplex() { return cplex_; }
  IloConstraintArray &constraints() { return constraints_; }
  IloNumArray &preferences() { return preferences_; }
  IloIntVar5dArray &x_mn_uvi() { return x_mn_uvi_; }
  IloIntVar2dArray &y_m_u() { return y_m_u_; }
  IloIntVar5dArray &z_uvi_pq() { return z_uvi_pq_; }
  IloIntVar3dArray &gamma_uvi() { return gamma_uvi_; }
  IloExpr &objective() { return objective_; }
  int k() { return k_; }
  void BuildModel();
  bool Solve();
  const static int kInfinity = 100000;

 private:
  // CPLEX related variables.
  IloEnv env_;
  IloModel model_;
  IloCplex cplex_;
  IloConstraintArray constraints_;
  IloNumArray preferences_;

  // Inputs.
  Graph *ip_topology_;
  Graph *otn_topology_;
  Graph *vn_topology_;
  std::vector<std::vector<int>> *location_constraint_;

  // ?
  int k_;
  int cost_new_ip_link_;

  // Decision variable for VLink to IP Link mapping.
  IloIntVar5dArray x_mn_uvi_;

  // Decision variable for creating new IP links.
  IloIntVar3dArray gamma_uvi_;

  // Decision variable for VNode to IP node mapping.
  IloIntVar2dArray y_m_u_;

  // Decision variable for IP link to OTN link mapping.
  IloIntVar5dArray z_uvi_pq_;

  // Input variable indicating location constraint.
  IloInt2dArray l_m_u_;

  // Input variable indicating existing IP links.
  IloInt3dArray ip_link_uvi_;

  // Input variable indicating IP node to OTN node attachment.
  IloInt2dArray tau_u_p_;

  // Bandwidth matrix for IP layer.
  std::vector<std::vector<std::vector<long>>> b_uvi_;

  // Cost matrix for IP layer.
  std::vector<std::vector<std::vector<int>>> cost_uvi_;

  // Objective function.
  IloExpr objective_;
};
#endif  // CPLEX_SOLVER_

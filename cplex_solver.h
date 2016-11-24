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
  MultiLayerVNESolver(Graph *pn_topology, Graph *vn_topology,
                      std::vector<std::vector<int>> *location_constraint,
                      int cost_new_ip_link, int k = 1);

  IloEnv &env() { return env_; }
  IloModel &model() { return model_; }
  IloCplex &cplex() { return cplex_; }
  IloConstraintArray &constraints() { return constraints_; }
  IloNumArray &preferences() { return preferences_; }
  IloIntVar5dArray &x_mn_uvi() { return x_mn_uvi_; }
  IloIntVar2dArray &y_m_u() { return y_m_u_; }
  IloExpr &objective() { return objective_; }
  int k() { return k_; }
  void BuildModel();
  bool Solve();
  const static int kInfinity = 100000;
 private:
  IloEnv env_;
  IloModel model_;
  IloCplex cplex_;
  IloConstraintArray constraints_;
  IloNumArray preferences_;
  Graph *pn_topology_;
  Graph *vn_topology_;
  OverlayMapping *ip_otn_mapping;
  std::vector<std::vector<int>> *location_constraint_;
  int k_;
  int cost_new_ip_link_;
  // Decision variable for edge mapping.
  IloIntVar5dArray x_mn_uvi_;
  // Decision variable for node mapping.
  IloIntVar2dArray y_m_u_;
  // Variable indicating location constraint.
  IloInt2dArray l_m_u_;
  // Objective function.
  IloExpr objective_;
};
#endif  // CPLEX_SOLVER_

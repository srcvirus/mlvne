#ifndef VNE_SOLUTION_BUILDER_H_
#define VNE_SOLUTION_BUILDER_H_

#include "cplex_solver.h"

class VNESolutionBuilder {
 public:
  VNESolutionBuilder(MultiLayerVNESolver *vne_solver_ptr, Graph *ip_topology,
                     Graph *otn_topology, Graph *vn_topology)
      : vne_solver_ptr_(vne_solver_ptr),
        ip_topology_(ip_topology),
        otn_topology_(otn_topology),
        vn_topology_(vn_topology) {}

  // Prints virtual node to IP node mapping on stdout. If filename is not NULL
  // then the same output is written to the corresponding file as well. Each
  // line in the output has the following format:
  // Virtual node <vnode_id> --> IP node <ip_node_id>
  void PrintVNodeMapping(const char *filename);

  // Prints virtual link to IP link mapping on stdout. If filename is not NULL
  // then the same output is written to the corresponding file as well. Each
  // line in the output has the following format:
  // Virtual link (<src>, <dst>) --> IP link (<src>, <dst>, <order>).
  void PrintVLinkMapping(const char *filename);

  // Prints status of running the solver, i.e., Optimal, Infeasible, etc.
  void PrintSolutionStatus(const char *filename);

  // Prints the value of objective function obtained by the solver.
  void PrintCost(const char *filename);

  // Prints new IP links and their mapping on OTN on stdout. If filename is not
  // NULL then output is written to the corresponding file as well. Each line in
  // the output has the following format:
  // New IP Link (<src>, <dst>, <order>) --> OTN link (<otn_src>, <otn_dst>).
  // Order is used to break ties for parallel links.
  void PrintNewIPLinks(const char *filename);

 private:
  MultiLayerVNESolver *vne_solver_ptr_;
  Graph *ip_topology_;
  Graph *otn_topology_;
  Graph *vn_topology_;
};

#endif  // VNE_SOLUTION_BUILDER_H_

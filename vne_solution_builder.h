#ifndef VNE_SOLUTION_BUILDER_H_
#define VNE_SOLUTION_BUILDER_H_

#include "cplex_solver.h"

class VNESolutionBuilder {
 public:
  VNESolutionBuilder(MultiLayerVNESolver *vne_solver_ptr, Graph *pn_topology,
                     Graph *vn_topology)
      : vne_solver_ptr_(vne_solver_ptr),
        pn_topology_(pn_topology),
        vn_topology_(vn_topology) {}

  void PrintMapping(const OverlayMapping *, const char *, const char *);
  void PrintEdgeMapping(const char *filename);
  void PrintNodeMapping(const char *filename);
  void PrintSolutionStatus(const char *filename);
  void PrintCost(const char *filename);
  void PrintNewIPLinks(const OverlayMapping *vne, const char *filename);
  unique_ptr<OverlayMapping> BuildVNEmbedding();
  unique_ptr<OverlayMapping> TranslateEmbeddingToIP(
      OverlayMapping *vne, Graph *ip_topology, OverlayMapping *ip_otn_mapping);

 private:
  MultiLayerVNESolver *vne_solver_ptr_;
  Graph *pn_topology_;
  Graph *vn_topology_;
};

#endif  // VNE_SOLUTION_BUILDER_H_

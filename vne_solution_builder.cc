#include "vne_solution_builder.h"
#include "util.h"

#include <algorithm>
#include <math.h>

unique_ptr<OverlayMapping> VNESolutionBuilder::BuildNewIPLinks(
    const OverlayMapping* vne, const OverlayMapping *ip_otn_mapping) {
  unique_ptr<OverlayMapping> new_ip_links(new OverlayMapping);
  edge_map_t::const_iterator vne_it;
  for (vne_it = vne->edge_map.begin(); vne_it != vne->edge_map.end(); ++vne_it) {
    const edge_t& vlink = vne_it->first;
    const path_t& mapped_path = vne_it->second;
    std::vector<path_t> paths;
    path_t current_path;
    for (int i = 0; i < mapped_path.size(); ++i) {
      if (pn_topology_->
          GetEdgeCost(
            mapped_path[i].first, 
            mapped_path[i].second, 
            mapped_path[i].order) == kInfinity) {
        current_path.push_back(mapped_path[i]);
      } else {
        if (!current_path.empty()) {
          paths.push_back(current_path);
          current_path.clear();
        }
      }
    }
    if (!current_path.empty()) 
      paths.push_back(current_path);
    for (int i = 0; i < paths.size(); ++i) {
      std::set<int> v;
      for (int j = 0; j < paths[i].size(); ++j) {
        if (pn_topology_->special_nodes()->find(paths[i][j].first) !=
            pn_topology_->special_nodes()->end()) {
          v.insert(paths[i][j].first);
        }
        if (pn_topology_->special_nodes()->find(paths[i][j].second) !=
            pn_topology_->special_nodes()->end()) {
          v.insert(paths[i][j].second);
        }
      }
      assert(v.size() == 2);
      std::set<int>::iterator it = v.begin();
      edge_t ip_link(0,0,0);
      ip_link.first = (*it++);
      ip_link.second = *it;
      ip_link.first = static_cast<int>(std::find(ip_otn_mapping->node_map.begin(),
                                ip_otn_mapping->node_map.end(), ip_link.first) - 
                      ip_otn_mapping->node_map.begin());
      ip_link.second = static_cast<int>(std::find(ip_otn_mapping->node_map.begin(),
                                ip_otn_mapping->node_map.end(), ip_link.second) -
                       ip_otn_mapping->node_map.begin());
      if (ip_link.first < ip_link.second) 
        std::swap(ip_link.first, ip_link.second);
      ip_link.order = 0;
      while (new_ip_links->edge_map.find(ip_link) != 
          new_ip_links->edge_map.end()) {
        ++ip_link.order;
      }
      new_ip_links->edge_map[ip_link] = paths[i];
    }
  }
  return boost::move(new_ip_links);
}

unique_ptr<OverlayMapping> VNESolutionBuilder::BuildVNEmbedding() {
  unique_ptr<OverlayMapping> vn_embedding(new OverlayMapping());
  vn_embedding->node_map.resize(vn_topology_->node_count());
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar5dArray &x_mn_uvi = vne_solver_ptr_->x_mn_uvi();
  const IloIntVar2dArray &y_m_u = vne_solver_ptr_->y_m_u();
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < pn_topology_->node_count(); ++u) {
      if (fabs(cplex.getValue(y_m_u[m][u]) - 1) < EPS) {
        vn_embedding->node_map[m] = u;
      }
    }
  }
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const std::vector<edge_endpoint> &m_neighbors = 
      vn_topology_->adj_list()->at(m);
    for (int i = 0; i < m_neighbors.size(); ++i) {
      const edge_endpoint &vend_point = m_neighbors[i];
      int n = vend_point.node_id;
      if (m < n) continue;
      edge_t vlink = edge_t(m, n);
      path_t mapped_path;
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        const std::vector<edge_endpoint> &u_neighbors = 
          pn_topology_->adj_list()->at(u);
        for (int j = 0; j < u_neighbors.size(); ++j) {
          const edge_endpoint &end_point = u_neighbors[j];
          int v = end_point.node_id;
          int order = end_point.order;
          if (fabs(cplex.getValue(x_mn_uvi[m][n][u][v][order]) - 1) < EPS) {
            mapped_path.push_back(edge_t(u, v, order));
          }
        }
      }
      int map_m = vn_embedding->node_map[m];
      int next_to_look = map_m;
      std::vector<edge_t> sorted_path;
      while (!mapped_path.empty()) {
        int k = -1;
        for (k = 0; k < mapped_path.size(); ++k) {
          if (mapped_path[k].first == next_to_look || 
              mapped_path[k].second == next_to_look) 
            break;
        }
        assert(k != -1 && k < mapped_path.size());
        sorted_path.push_back(mapped_path[k]);
        next_to_look = mapped_path[k].first == next_to_look? mapped_path[k].second : mapped_path[k].first;
        mapped_path.erase(mapped_path.begin() + k);
      }
      vn_embedding->edge_map[vlink] = sorted_path;
    }
  }
  return boost::move(vn_embedding);
}

void VNESolutionBuilder::PrintEdgeMapping(const char *filename) {
  FILE *outfile = NULL;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar5dArray &x_mn_uvi = vne_solver_ptr_->x_mn_uvi();
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const std::vector<edge_endpoint> &m_neighbors = 
      vn_topology_->adj_list()->at(m);
    for (int i = 0; i < m_neighbors.size(); ++i) {
      const edge_endpoint &vend_point = m_neighbors[i];
      int n = vend_point.node_id;
      if (m < n) continue;
      for (int u = 0; u < pn_topology_->node_count(); ++u) {
        const std::vector<edge_endpoint> &u_neighbors = 
          pn_topology_->adj_list()->at(u);
        for (int j = 0; j < u_neighbors.size(); ++j) {
          const edge_endpoint &end_point = u_neighbors[j];
          int v = end_point.node_id;
          int order = end_point.order;
          if (fabs(cplex.getValue(x_mn_uvi[m][n][u][v][order]) - 1) < EPS) {
            printf("Virtual edge (%d, %d) --> pn edge (%d, %d, %d)\n", m, n,
                   u, v, order);
            if (outfile) {
              fprintf(outfile,
                      "Virtual edge (%d, %d) --> pn edge (%d, %d, %d)\n", m,
                      n, u, v, order);
            }
          }
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintNodeMapping(const char *filename) {
  FILE *outfile = NULL;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar2dArray &y_m_u = vne_solver_ptr_->y_m_u();
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < pn_topology_->node_count(); ++u) {
      if (fabs(cplex.getValue(y_m_u[m][u]) - 1) < EPS) {
        printf("Virtual node %d --> pn node %d\n", m, u);
        if (outfile) {
          fprintf(outfile, "Virtual node %d --> pn node %d\n", m, u);
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintNewIPLinks(
    const OverlayMapping *vne, const char *filename) {
  FILE* outfile;
  if (filename) outfile = fopen(filename, "w");
  if (vne == NULL) return;
  std::map<edge_t, path_t>::const_iterator it;
  for (it = vne->edge_map.begin(); it != vne->edge_map.end(); ++it) {
    const edge_t& ip_link = it->first;
    const path_t& otn_path = it->second;
    for (int i = 0; i < otn_path.size(); ++i) {
      printf("New IP Link (%d, %d, %d) --> OTN link (%d, %d, %d)\n",
          ip_link.first, ip_link.second, ip_link.order, otn_path[i].first,
          otn_path[i].second, otn_path[i].order);
      if (outfile) {
        fprintf(outfile, "New IP Link (%d, %d, %d) --> OTN link (%d, %d, %d)\n",
            ip_link.first, ip_link.second, ip_link.order, otn_path[i].first,
            otn_path[i].second, otn_path[i].order);
      }
    }
  }
}

void VNESolutionBuilder::PrintSolutionStatus(const char *filename) {
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  std::cout << "Solution status = " << cplex.getStatus() << std::endl;
  if (filename) {
    std::ofstream ofs(filename);
    ofs << cplex.getStatus();
    ofs.close();
  }
}

void VNESolutionBuilder::PrintCost(const char *filename) {
  FILE *outfile = NULL;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  printf("Cost = %lf\n", cplex.getObjValue());
  if (outfile) {
    fprintf(outfile, "%lf\n", cplex.getObjValue());
  }
}

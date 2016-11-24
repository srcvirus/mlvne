#include "cplex_solver.h"
#include "datastructure.h"
#include "io.h"
#include "util.h"
#include "vne_solution_builder.h"


const std::string kUsage =
    "./mlvn "
    "--otn_topology_file=<otn_topology_file>\n"
    "\t--ip_topology_file=<ip_topology_file>\n"
    "\t--ip_node_mapping_file=<ip_mapping_file>\n"
    "\t--ip_link_mapping_file=<ip_link_mapping_file>\n"
    "\t--ip_port_info_file=<ip_port_info_file>\n"
    "\t--vn_topology_file=<vn_topology_file>\n"
    "\t--vn_location_file=<vn_location_file>\n";

unique_ptr<Graph> TransformSubstrateTopology(
    const Graph* otn_topology, const Graph* ip_topology,
    const OverlayMapping* ip_otn_mapping) {
  unique_ptr<Graph> pn_topology(new Graph());
  for (int i = 0; i < otn_topology->node_count(); ++i) {
    const std::vector<edge_endpoint>& adj_list =
        otn_topology->adj_list()->at(i);
    for (int j = 0; j < adj_list.size(); ++j) {
      const edge_endpoint& endpoint = adj_list[j];
      if (i > endpoint.node_id)
        continue;
      pn_topology->AddEdge(i, endpoint.node_id, endpoint.bandwidth,
                           endpoint.delay, endpoint.cost, false);
    }
  }
  for (int i = 0; i < ip_topology->node_count(); ++i) {
    const std::vector<edge_endpoint>& adj_list = ip_topology->adj_list()->at(i);
    for (int j = 0; j < adj_list.size(); ++j) {
      const edge_endpoint& endpoint = adj_list[j];
      if (i > endpoint.node_id) 
        continue;
      int otn_u = ip_otn_mapping->node_map[i];
      int otn_v = ip_otn_mapping->node_map[endpoint.node_id];
      pn_topology->AddSpecialNode(otn_u);
      pn_topology->AddSpecialNode(otn_v);

      // Update residual bandwidth on the collapsed graph. Substract the IP
      // link's bandwidth from the residual bandwidth of the links on the mapped
      // OTN path.
      const path_t& mapped_otn_path =
          ip_otn_mapping->edge_map.find(edge_t(i, endpoint.node_id))->second;
      DEBUG("mapped path length for %d, %d = %d\n", i, endpoint.node_id, mapped_otn_path.size());
      int ip_link_cost = 0;
      for (int k = 0; k < mapped_otn_path.size(); ++k) {
        int u = mapped_otn_path[k].first, v = mapped_otn_path[k].second;
        if ((u == 0 || u == 49) && (v == 0 || v == 49)) {
            printf("BINGO\n");
        }
        ip_link_cost += otn_topology->GetEdgeCost(u, v);
        long bw = pn_topology->GetEdgeResidualBandwidth(u, v);
        pn_topology->SetEdgeResidualBandwidth(u, v, bw - endpoint.bandwidth);
        pn_topology->SetEdgeResidualBandwidth(v, u, bw - endpoint.bandwidth);
      }
      // TODO: This will add both otn_u, otn_v and otn_v, otn_u to the graph.
      // Need to handle this.
      pn_topology->AddEdge(otn_u, otn_v, endpoint.bandwidth, endpoint.delay,
                           ip_link_cost, true);
    }
  }
  return boost::move(pn_topology);
}

int main(int argc, char* argv[]) {
  if (argc < 8) {
    printf("Not sufficient arguments. Expected 8, provided %d\n", argc);
    printf("%s\n", kUsage.c_str());
    return 1;
  }
  // Parse the command line arguments.
  using std::string;
  unique_ptr<std::map<string, string>> arg_map(ParseArgs(argc, argv).release());
  string otn_topology_file = "", ip_topology_file = "", ip_port_info_file = "",
         ip_node_mapping_file = "", ip_link_mapping_file = "",
         vn_topology_file = "", vn_location_file = "";
  std::map<string, string>::iterator arg_map_it;
  for (arg_map_it = arg_map->begin(); arg_map_it != arg_map->end();
       ++arg_map_it) {
    if (arg_map_it->first == "--otn_topology_file") {
      otn_topology_file = arg_map_it->second;
    } else if (arg_map_it->first == "--ip_topology_file") {
      ip_topology_file = arg_map_it->second;
    } else if (arg_map_it->first == "--ip_node_mapping_file") {
      ip_node_mapping_file = arg_map_it->second;
    } else if (arg_map_it->first == "--ip_link_mapping_file") {
      ip_link_mapping_file = arg_map_it->second;
    } else if (arg_map_it->first == "--ip_port_info_file") {
      ip_port_info_file = arg_map_it->second;
    } else if (arg_map_it->first == "--vn_topology_file") {
      vn_topology_file = arg_map_it->second;
    } else if (arg_map_it->first == "--vn_location_file") {
      vn_location_file = arg_map_it->second;
    } else {
      printf("Unrecognized command line option: %s\n", arg_map_it->first.c_str());
      printf("%s\n", kUsage.c_str());
    }
  }

  // Initialize the inputs.
  unique_ptr<Graph> otn_topology(
      InitializeTopologyFromFile(otn_topology_file.c_str()).release());
  unique_ptr<Graph> ip_topology(
      InitializeTopologyFromFile(ip_topology_file.c_str()).release());
  unique_ptr<std::vector<std::vector<int>>> ip_port_info(
      InitializePortInfoFromFile(ip_port_info_file.c_str()).release());
  unique_ptr<OverlayMapping> ip_otn_mapping(
      InitializeOverlayMappingFromFile(ip_node_mapping_file.c_str(),
                                       ip_link_mapping_file.c_str())
          .release());
  unique_ptr<Graph> vn_topology(
      InitializeTopologyFromFile(vn_topology_file.c_str()).release());
  unique_ptr<std::vector<std::vector<int>>> location_constraint(
      InitializeVNLocationsFromFile(vn_location_file.c_str(),
                                    vn_topology->node_count())
          .release());

  for (int i = 0; i < ip_port_info->size(); ++i) {
    int u = ip_port_info->at(i)[0];
    int port_count = ip_port_info->at(i)[1];
    int port_capacity = ip_port_info->at(i)[2];
    ip_topology->SetPortCount(u, port_count);
    ip_topology->SetPortCapacity(u, port_capacity);
  }

  // Collaps IP and OTN into one layer.
  unique_ptr<Graph> phys_topology(
      TransformSubstrateTopology(otn_topology.get(), ip_topology.get(),
                                 ip_otn_mapping.get())
          .release());

  // Create a new location constraint set w.r.t the collapsed SN.
  unique_ptr<std::vector<std::vector<int> > > xlocation_constraint(
      new std::vector<std::vector<int> >(location_constraint->size()));
  for (int i = 0; i < location_constraint->size(); ++i) {
    for (int j = 0; j < location_constraint->at(i).size(); ++j) {
      int loc = location_constraint->at(i)[j];
      DEBUG("prev_loc = %d, new_loc = %d\n", loc, ip_otn_mapping->node_map[loc]);
      xlocation_constraint->at(i).push_back(ip_otn_mapping->node_map[loc]);
    }
  }
  // printf("%s\n", phys_topology->GetDebugString().c_str());
  unique_ptr<MultiLayerVNESolver> mlvne_solver(
      new MultiLayerVNESolver(phys_topology.get(), vn_topology.get(), 
        xlocation_constraint.get(), kInfinity, 5));
  VNESolutionBuilder vne_sbuilder(mlvne_solver.get(), phys_topology.get(), vn_topology.get());
  mlvne_solver->BuildModel();
  if (mlvne_solver->Solve()) {
    unique_ptr<OverlayMapping> vne(vne_sbuilder.BuildVNEmbedding().release());
    unique_ptr<OverlayMapping> new_ip_links(
        vne_sbuilder.BuildNewIPLinks(vne.get(), ip_otn_mapping.get()));
    vne_sbuilder.PrintNodeMapping((vn_topology_file + ".nmap").c_str());
    vne_sbuilder.PrintEdgeMapping((vn_topology_file + ".emap").c_str());
    vne_sbuilder.PrintCost((vn_topology_file + ".cost").c_str());
    vne_sbuilder.PrintNewIPLinks(new_ip_links.get(), (vn_topology_file + ".new_ip").c_str());
  }
  vne_sbuilder.PrintSolutionStatus((vn_topology_file + ".status").c_str());
  return 0;
}

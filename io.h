#ifndef IO_H_
#define IO_H_

#include "datastructure.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <memory>
#include <string>

using std::unique_ptr;

typedef std::vector<std::vector<std::string> > csv_vector_t;
typedef unique_ptr<csv_vector_t> csv_vector_ptr_t;

unique_ptr<std::map<std::string, std::string> > ParseArgs(int argc,
                                                         char* argv[]) {
  unique_ptr<std::map<std::string, std::string> > arg_map(
      new std::map<std::string, std::string>());
  for (int i = 1; i < argc; ++i) {
    char* key = strtok(argv[i], "=");
    char* value = strtok(NULL, "=");
    DEBUG(" [%s] => [%s]\n", key, value);
    arg_map->insert(std::make_pair(key, value));
  }
  return std::move(arg_map);
}

csv_vector_ptr_t ReadCSVFile(const char* filename) {
  DEBUG("[Parsing %s]\n", filename);
  FILE* file_ptr = fopen(filename, "r");
  if (!file_ptr) {
    DEBUG("Invalid file %s\n", filename);
    return nullptr;
  }
  const static int kBufferSize = 1024;
  char line_buffer[kBufferSize];
  csv_vector_ptr_t ret_vector(new csv_vector_t());
  std::vector<std::string> current_line;
  int row_number = 0;
  while (fgets(line_buffer, kBufferSize, file_ptr)) {
    DEBUG("Read %d characters\n", strlen(line_buffer));
    if (strlen(line_buffer) <= 0) continue;
    if (line_buffer[0] == '\n' || line_buffer[0] == '\r') continue;
    current_line.clear();
    char* token = strtok(line_buffer, ",\n\r");
    current_line.push_back(token);
    while ((token = strtok(NULL, ",\n"))) {
      current_line.push_back(token);
    }
    ret_vector->push_back(current_line);
  }
  fclose(file_ptr);
  DEBUG("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
  return std::move(ret_vector);
}

unique_ptr<Graph> InitializeTopologyFromFile(const char* filename) {
  int node_count = 0, edge_count = 0;
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return nullptr; 
  }
  unique_ptr<Graph> graph(new Graph());
  node_count = std::stoi(csv_vector->at(0)[0]);
  for (int i = 1; i < csv_vector->size(); ++i) {
    const auto& row = csv_vector->at(i);

    // Each line has the following format:
    // SourceID, DestinationID, Cost, Bandwidth.
    int u = std::stoi(row[0]);
    int v = std::stoi(row[1]);
    int cost = std::stoi(row[2]);
    long bw = std::stol(row[3]);

    DEBUG("Line[%d]: u = %d, v = %d, cost = %d, bw = %ld\n", i, u, v, cost,
          bw);
    graph->AddEdge(u, v, bw, cost, false);
  }
  return std::move(graph);
}

unique_ptr<std::vector<std::vector<int> > > InitializeVNLocationsFromFile(
    const char* filename, int num_virtual_nodes) {
  DEBUG("Parsing %s\n", filename);
  unique_ptr<std::vector<std::vector<int> > > ret_vector(
      new std::vector<std::vector<int> >(num_virtual_nodes));
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return nullptr;
  }
  DEBUG("Parsing %s successful\n", filename);
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);
    int vnode_id = atoi(row[0].c_str());
    for (int j = 1; j < row.size(); ++j) {
      ret_vector->at(vnode_id).push_back(atoi(row[j].c_str()));
    }
  }
  return std::move(ret_vector);
}

unique_ptr<OverlayMapping> InitializeOverlayMappingFromFile(
    const char* nmap_file, const char* emap_file) {
  DEBUG("Reading node embedding from: %s\n", nmap_file);
  DEBUG("Reading edge embedding from: %s\n", emap_file);
  unique_ptr<OverlayMapping> mapping(new OverlayMapping());

  // Initialize node mapping.
  csv_vector_ptr_t nmap_csv_vector = ReadCSVFile(nmap_file);
  if (nmap_csv_vector.get() == NULL) {
    return nullptr;
  }
  for (int i = 0; i < nmap_csv_vector->size(); ++i) {
    const std::vector<std::string>& row = nmap_csv_vector->at(i);
    int vnode = std::stoi(row[0].c_str());
    int vnode_map = std::stoi(row[1].c_str());
    if (vnode > static_cast<int>(mapping->node_map.size()) - 1)
      mapping->node_map.resize(vnode + 1);
    mapping->node_map[vnode] = vnode_map;
  }

  // Initialize edge mapping.
  csv_vector_ptr_t emap_csv_vector = ReadCSVFile(emap_file);
  if (emap_csv_vector.get() == NULL) {
    return nullptr;
  }
  for (int i = 0; i < emap_csv_vector->size(); ++i) {
    const std::vector<std::string>& row = emap_csv_vector->at(i);
    int m = std::stoi(row[0].c_str());
    int n = std::stoi(row[1].c_str());
    int u = std::stoi(row[2].c_str());
    int v = std::stoi(row[3].c_str());
    edge_t overlay_link(m, n), underlay_link(u, v);
    if (mapping->edge_map.find(overlay_link) == mapping->edge_map.end()) {
      mapping->edge_map[overlay_link] = path_t();
    }
    mapping->edge_map[overlay_link].push_back(underlay_link);
    DEBUG("Current embedding path length of (%d, %d) is %u\n", m, n,
          mapping->edge_map[overlay_link].size());
  }
  DEBUG("Embedding of %d links read successfully\n", mapping->edge_map.size());
  return std::move(mapping);
}

unique_ptr<std::vector<std::vector<int> > > InitializePortInfoFromFile(
    const char* port_info_file) {
  unique_ptr<std::vector<std::vector<int> > > port_info(
      new std::vector<std::vector<int> >());
  csv_vector_ptr_t csv_vector = ReadCSVFile(port_info_file);
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);
    int u = std::stoi(row[0].c_str());
    int num_ports = std::stoi(row[1].c_str());
    int port_capacity = std::stoi(row[2].c_str());
    std::vector<int> v;
    v.push_back(u);
    v.push_back(num_ports);
    v.push_back(port_capacity);
    port_info->push_back(v);
  }
  return std::move(port_info);
}

#endif  // IO_H_

#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <math.h>
#include <stdlib.h>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define INF 99999999
#define MAXN 1000
#define MAX_PARALLEL_LINKS 10
#define NIL -1

using std::unique_ptr;

// Type definition for convenience.
typedef struct edge_t {
  int first, second, order;
  edge_t(int f, int s, int o = 0) : first(f), second(s), order(o) {}
  bool operator<(const edge_t& e) const {
    if (first != e.first) return first < e.first;
    if (second != e.second) return second < e.second;
    return order < e.order;
  }
  bool operator==(const edge_t& e) const {
    return first == e.first && second == e.second && order == e.order;
  }
} edge_t;

typedef std::vector<edge_t> path_t;
typedef std::set<edge_t> edge_set_t;
typedef std::map<edge_t, path_t> edge_map_t;

inline edge_t ConstructEdge(int u, int v) {
  if (u > v) std::swap(u, v);
  return edge_t(u, v);
}

template <typename T>
struct matrix_t {
  std::vector<std::vector<T> > matrix;
  matrix_t() {}
  matrix_t(int rows, int columns, T fill_value = T())
      : matrix(rows, std::vector<T>(columns, fill_value)) {}
};

// An entry in an adjacent list. An entry contains the node_id of the endpoint.
// The entry contains bandwidth, residual bandwidth and cost of the
// corresponding edge.
struct edge_endpoint {
  int node_id;
  int order;  // In case of parallel links, order is a tiebreaker.
  long bandwidth;
  long residual_bandwidth;
  int cost;
  bool is_pseudo_endpoint;
  edge_endpoint()
      : node_id(NIL),
        bandwidth(0),
        residual_bandwidth(0),
        cost(INF),
        is_pseudo_endpoint(false) {}
  edge_endpoint(int node_id, int order, long bw, int cost,
                bool is_pseudo_endpoint)
      : node_id(node_id),
        order(order),
        bandwidth(bw),
        residual_bandwidth(bw),
        cost(cost),
        is_pseudo_endpoint(is_pseudo_endpoint) {}
  std::string GetDebugString() const {
    return "node_id = " + std::to_string(node_id) +
           ", order = " + std::to_string(order) +
           ", bandwidth = " + std::to_string(bandwidth) +
           ", residual_bandwidth = " +
           std::to_string(residual_bandwidth) + ", cost = " +
           std::to_string(cost) + ", is_pseudo_endpoint = " +
           std::to_string(is_pseudo_endpoint);
  }
};

// Graph data structure. Represents an undirected graph. Internally, the graph
// is represented as an adjacency list.
class Graph {
 public:
  Graph() {
    adj_list_ = unique_ptr<std::vector<std::vector<edge_endpoint> > >(
        new std::vector<std::vector<edge_endpoint> >);
    port_counts_ = unique_ptr<std::vector<int> >(new std::vector<int>());
    port_capacities_ = unique_ptr<std::vector<int> >(new std::vector<int>());
    special_nodes_ = unique_ptr<std::set<int> >(new std::set<int>());
    node_count_ = edge_count_ = 0;
  }

  Graph(const Graph& g) {
    this->node_count_ = g.node_count_;
    this->edge_count_ = g.edge_count_;
    this->adj_list_ = unique_ptr<std::vector<std::vector<edge_endpoint> > >(
        new std::vector<std::vector<edge_endpoint> >(*g.adj_list_.get()));
    this->port_counts_ = unique_ptr<std::vector<int> >(
        new std::vector<int>(*g.port_counts_.get()));
    this->port_capacities_ = unique_ptr<std::vector<int> >(
        new std::vector<int>(*g.port_capacities_.get()));
    this->special_nodes_ =
        unique_ptr<std::set<int> >(new std::set<int>(*g.special_nodes_.get()));
  }

  // Accessor methods.
  int node_count() const { return node_count_; }
  int edge_count() const { return edge_count_; }
  const std::vector<std::vector<edge_endpoint> >* adj_list() const {
    return static_cast<const std::vector<std::vector<edge_endpoint> >* >(
        adj_list_.get());
  }
  const std::set<int>* special_nodes() const {
    return static_cast<const std::set<int>*>(special_nodes_.get());
  }

  void AddSpecialNode(int u) { special_nodes_->insert(u); }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::AddEdge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  void AddEdge(int u, int v, long bw, int cost, bool is_pseudo_endpoint) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    int order = 0;
    const auto& neighbors = adj_list_->at(u);
    for (const auto end_point : neighbors) {
      if (end_point.node_id == v) ++order;
    }
    adj_list_->at(u).push_back(
        edge_endpoint(v, order, bw, cost, is_pseudo_endpoint));
    adj_list_->at(v).push_back(
        edge_endpoint(u, order, bw, cost, is_pseudo_endpoint));
    ++edge_count_;
    node_count_ = adj_list_->size();
  }

  bool IsPseudoEdge(int u, int v, int order = 0) const {
    const auto& neighbors = adj_list_->at(u);
    for (const auto end_point : neighbors) {
      if (end_point.node_id == v && end_point.order == order)
        return end_point.is_pseudo_endpoint;
    }
  }

  int GetEdgeCost(int u, int v, int order = 0) const {
    const auto& neighbors = adj_list_->at(u);
    for (const auto end_point : neighbors) {
      if (end_point.node_id == v && end_point.order == order)
        return end_point.cost;
    }
    return -1;
  }

  long GetEdgeBandwidth(int u, int v, int order = 0) const {
    const auto& neighbors = adj_list_->at(u);
    for (const auto end_point : neighbors) {
      if (end_point.node_id == v && end_point.order == order)
        return end_point.bandwidth;
    }
    return -1;
  }

  void SetEdgeBandwidth(int u, int v, long bw, int order = 0) {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id == v && end_point_it->order == order) {
        end_point_it->bandwidth = bw;
        break;
      }
    }
  }

  long GetEdgeResidualBandwidth(int u, int v, int order = 0) const {
    const auto& neighbors = adj_list_->at(u);
    for (const auto end_point : neighbors) {
      if (end_point.node_id == v && end_point.order == order)
        return end_point.residual_bandwidth;
    }
    return -1;
  }

  void SetEdgeResidualBandwidth(int u, int v, long rbw, int order = 0) {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point_it = neighbors.begin(); end_point_it != neighbors.end();
        ++end_point_it) {
      if (end_point_it->node_id == v && end_point_it->order == order) {
        end_point_it->residual_bandwidth = rbw;
        break;
      }
    }
  }

  // Returns the cumulative bandwidth of all links originating at node u.
  long GetTotalNodeBandwidth(int u) const {
    const auto& neighbors = adj_list_->at(u);
    long total_bw = 0;
    for (const auto end_point : neighbors) {
      total_bw += end_point.bandwidth;
    }
    return total_bw;
  }

  inline int GetNodeDegree(int u) const { return adj_list_->at(u).size(); }
  inline int GetPortCapacity(int u) const { return port_capacities_->at(u); }
  inline int GetPortCount(int u) const { port_counts_->at(u); }
  inline int GetResidualPortCount(int u) const {
    return port_counts_->at(u) - GetNodeDegree(u);
  }

  void SetPortCount(int u, int port_count) {
    if (port_counts_->size() <= u) {
      port_counts_->resize(u + 1);
    }
    port_counts_->at(u) = port_count;
  }

  void SetPortCapacity(int u, int port_capacity) {
    if (port_capacities_->size() <= u) {
      port_capacities_->resize(u + 1);
    }
    port_capacities_->at(u) = port_capacity;
  }

  std::string GetDebugString() const {
    std::string ret_string =
        "node_count = " + std::to_string(node_count_);
    ret_string += ", edge_count = " +
                  std::to_string(edge_count_) + "\n";
    for (int i = 0; i < node_count_; ++i) {
      const auto& neighbors = adj_list_->at(i);
      ret_string += std::to_string(i) + " --> ";
      for (const auto neighbor : neighbors) {
        ret_string += "\t(" + neighbor.GetDebugString() + ")\n";
      }
      ret_string += "\n";
    }
    return ret_string;
  }

  virtual ~Graph() { adj_list_.reset(); }

 private:
  unique_ptr<std::vector<std::vector<edge_endpoint> > > adj_list_;
  int node_count_, edge_count_, total_port_count_;
  unique_ptr<std::vector<int> > port_capacities_;
  unique_ptr<std::vector<int> > port_counts_;
  unique_ptr<std::set<int> > special_nodes_;
};

struct OverlayMapping {
  std::vector<int> node_map;
  edge_map_t edge_map;
  OverlayMapping() {}
  OverlayMapping(const std::vector<int>& nmap, const edge_map_t& emap)
      : node_map(nmap), edge_map(emap) {}
  virtual ~OverlayMapping() {}
};

#endif  // DATASTRUCTURE_H_

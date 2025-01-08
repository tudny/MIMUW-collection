#include "tbb/tbb.h"
#include <iostream>
#include <math.h>
#include <map>
#include <set>
#include <random>


typedef int node_t;
typedef long weight_t;
typedef std::pair<node_t, weight_t> edge_t;


// Our graph just maps a node onto its neighbors.
typedef std::map<node_t, std::set<edge_t>> graph_t;


// Initialize a DAG with random adjacencies.
void rand_init_DAG_graph(graph_t& graph, int node_count,
                         double edge_probability) {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> dis{0, 1};
    for (int i = 0; i < node_count; ++i) {
        auto neighbors = std::set<edge_t>();
        for (int j = i+1; j < node_count; ++j) {
            if (dis(gen) < edge_probability) {
                auto weight = static_cast<weight_t>(dis(gen) * 100);
                neighbors.insert({j, weight});
            }
        }
        graph[i] = neighbors;
    }
}

void dijkstra(
    graph_t &graph,
    node_t source,
    std::map<node_t, weight_t> &min_distance_result
) {
    tbb::concurrent_hash_map<node_t, weight_t> min_distance;
    char c = 0;

    tbb::concurrent_priority_queue<std::pair<weight_t, node_t>> q;
    q.push({-0, source});

    tbb::parallel_for_each(&c, &c + 1, [&](char, tbb::feeder<char>& feeder) {
        std::pair<weight_t, node_t> current_element;
        q.try_pop(current_element);

        auto current_node = current_element.second;
        auto current_dis = -current_element.first;

        decltype(min_distance)::const_accessor real_current_dis;
        if (min_distance.find(real_current_dis, current_node) && real_current_dis->second < current_dis) {
            return;
        }

        for (const auto &kid : graph[current_node]) {
            const auto &[neighbor, weight] = kid;
            auto new_dis = current_dis + weight;
            decltype(min_distance)::accessor real_new_dis;
            if (!min_distance.insert(real_new_dis, neighbor)) {
                if (real_new_dis->second > new_dis) {
                    real_new_dis->second = new_dis;
                    q.push({-new_dis, neighbor});
                    feeder.add(0);
                }
            } else {
                real_new_dis->second = new_dis;
                q.push({-new_dis, neighbor});
                feeder.add(0);
            }
        }
    });

    for (const auto &pair : min_distance) {
        min_distance_result[pair.first] = pair.second;
    }
}



int main(int argc, char* argv[]) {
    const int node_count = 10;
    graph_t graph;
    rand_init_DAG_graph(graph, node_count, 0.5);
    node_t node = 0;

    for (const auto &node : graph) {
        std::cout << "Node " << node.first << " has neighbors: ";
        for (const auto &neighbor : node.second) {
            std::cout << "(" << neighbor.first << ", " << neighbor.second << ") ";
        }
        std::cout << std::endl;
    }

    std::map<node_t, weight_t> min_distance;
    dijkstra(graph, node, min_distance);

    for (const auto &node : graph) {
        std::cout << "Node " << node.first << " has distance " << min_distance[node.first] << std::endl;
    }

    return 0;
}

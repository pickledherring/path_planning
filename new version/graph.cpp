#include <vector>
#include <random>
#include "graph.h"
using namespace std;

struct Edge;
struct Node {
    float coords[2];
    std::vector<Edge> conns;
};

struct Edge {
    // 0 - 1 scalar
    float weight;
    Node* start, * end;
};

Graph::Graph(float corner_coords[4])
    :x1{corner_coords[0]},
    y1{corner_coords[1]},
    x2{corner_coords[2]},
    y2{corner_coords[3]},
    nodes(0), size{0}
    {};

void Graph::add_nodes(int n, float min_dist, float min_angle) {
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<mt19937::result_type> lat(y1, y2);
    uniform_real_distribution<mt19937::result_type> lon(x1, x2);
    
    for (int i = 0; i < n; i++) {
        Node point;
        point.coords[0] = lat(rng);
        point.coords[1] = lon(rng);
        // check distances


        nodes.push_back(point);
    }

    for (int i = 0; i < n; i++) {
        // add edges

    }

}
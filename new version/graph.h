#include <vector>

struct Edge;
struct Node;

class Graph {
    public:
        // square board determined by corner_coords (x1, y1, x2, y2)
        Graph(float corner_coords[4]);
        // n: number of nodes
        // min_dist: min euclidean distance between nodes
        // min_angle: min angular distance to add an edge
        void add_nodes(int n, float min_dist, float min_angle);
    private:
        float x1, y1, x2, y2;
        std::vector<Node> nodes;
        int size;
};
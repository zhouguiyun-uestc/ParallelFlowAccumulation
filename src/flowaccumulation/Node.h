#ifndef NODE_HEAD_H
#define NODE_HEAD_H

#include <functional>

class Node {
public:
    int row;
    int col;
    double elevation;

    Node() {
        row = 0;
        col = 0;
        elevation = -9999.0;
    }
    Node(int row, int col, float elevation) {
        this->row = row;
        this->col = col;
        this->elevation = elevation;
    }

    struct Greater : public std::binary_function<Node, Node, bool> {
        bool operator()(const Node n1, const Node n2) const {
            return n1.elevation > n2.elevation;
        }
    };

    bool operator==(const Node& a) {
        return (this->col == a.col) && (this->row == a.row);
    }
    bool operator!=(const Node& a) {
        return (this->col != a.col) || (this->row != a.row);
    }
    bool operator<(const Node& a) {
        return this->elevation < a.elevation;
    }
    bool operator>(const Node& a) {
        return this->elevation > a.elevation;
    }
    bool operator>=(const Node& a) {
        return this->elevation >= a.elevation;
    }
    bool operator<=(const Node& a) {
        return this->elevation <= a.elevation;
    }
};

#endif
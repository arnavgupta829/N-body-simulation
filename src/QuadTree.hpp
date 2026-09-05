#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include "StructDefs.hpp"
#include <vector>
#include <atomic>

class QuadTree {
public:
    double theta;
    void build(std::vector<Particle>& particles, double center[2], double halfWidth);
    void buildParallel(std::vector<Particle>& particles, double center[2], double halfWidth);
    void computeForce(Particle& particle) const;
    static void computeForceBrute(Particle& target, const std::vector<Particle>& particles);
    void clear();
    const QuadTreeNode* getRoot() const { return root; }

private:
    QuadTreeNode* root = nullptr;
    std::vector<QuadTreeNode> node_pool;
    // Thread-safe node pool: atomic counter allows lock-free allocation
    std::atomic<std::size_t> pool_index{0};

    QuadTreeNode* allocateNode(double center[2], double halfWidth);
    void insert(QuadTreeNode* node, Particle* particle, int depth);
    void computeMassDistribution(QuadTreeNode* node);
    void computeForceRecursive(const QuadTreeNode* node, Particle& target) const;
    static int getQuadrant(const double pos[2], const double center[2]);
    void buildSubtreeParallel(QuadTreeNode* node, Particle* particles, int count, int depth, Particle* tempBuffer);
    void buildSubtreeSequential(QuadTreeNode* node, Particle* particles, int count, int depth);
    static void partitionByQuadrant(Particle* particles, int count, const double center[2], int quadrantCounts[4], Particle* tempBuffer);
    std::vector<Particle> partitionBuffer;
};

#endif

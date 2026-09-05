#include "QuadTree.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>

/**
 * In case we build the tree parallelly, we want a thread-safe way of adding 
 * 
 */
QuadTreeNode* QuadTree::allocateNode(double center[2], double halfWidth) {
    size_t idx = pool_index.fetch_add(1, std::memory_order_relaxed);

    // Ideally, code should never reach line 12 since we have pre-allocated sufficient
    // memory. This is just added for graceful handling
    if (idx >= node_pool.size()) {
        return nullptr;
    }

    /**
     * Placement new here should save use time
     */
    QuadTreeNode* node = &node_pool[idx];
    node->center[0] = center[0];
    node->center[1] = center[1];
    node->halfWidth = halfWidth;
    node->totalMass = 0.0;
    node->com[0] = 0.0;
    node->com[1] = 0.0;
    node->children[0] = nullptr;
    node->children[1] = nullptr;
    node->children[2] = nullptr;
    node->children[3] = nullptr;
    node->particle = nullptr;
    node->isLeaf = true;

    return node;
}

/**
 * We start from quadrant neg-x, pos-y as 0 (left-most and top-most) and
 * move along the space based on the x and y positions
 */
int QuadTree::getQuadrant(const double pos[2], const double center[2]) {
    int quadrant = 0;
    if (pos[0] >= center[0]) quadrant += 1;
    if (pos[1] < center[1])  quadrant += 2;
    return quadrant;
}

void QuadTree::build(std::vector<Particle>& particles, double center[2], double half_width) {
    clear();
    node_pool.resize(std::max((size_t)(4 * particles.size()), (size_t)100));

    root = allocateNode(center, half_width);

    for (auto& p : particles) {
        insert(root, &p, 0);
    }

    computeMassDistribution(root);
}

void QuadTree::clear() {
    root = nullptr;

    // point first index to fetch back to 0
    pool_index.store(0, std::memory_order_relaxed);
}

void QuadTree::insert(QuadTreeNode* node, Particle* particle, int depth) {
    if (depth > Constants::MAX_TREE_DEPTH) {
        node->particle = particle;
        return;
    }

    if (node->isLeaf && node->particle == nullptr) {
        node->particle = particle;
        return;
    }

    /**
     * What we previously thought to be a leaf, now turns out to have two particles
     * in that same quadrant.
     * 
     * In this case, convert the leaf into an internal node, split it further, send
     * the initial particle deeper into the tree
     */
    if (node->isLeaf && node->particle != nullptr) {
        Particle* existing = node->particle;
        node->particle = nullptr;
        node->isLeaf = false;

        int eq = getQuadrant(existing->pos, node->center);
        if (node->children[eq] == nullptr) {
            double childHalfWidth = node->halfWidth / 2.0;
            double childCenter[2];
            double xOff = (eq & 1) ? childHalfWidth : -childHalfWidth;
            double yOff = (eq & 2) ? -childHalfWidth : childHalfWidth;
            childCenter[0] = node->center[0] + xOff;
            childCenter[1] = node->center[1] + yOff;
            node->children[eq] = allocateNode(childCenter, childHalfWidth);
        }
        insert(node->children[eq], existing, depth + 1);
    }

    /**
     * Now insert the new particle
     */
    int q = getQuadrant(particle->pos, node->center);
    if (node->children[q] == nullptr) {
        double childHalfWidth = node->halfWidth / 2.0;
        double childCenter[2];
        double xOff = (q & 1) ? childHalfWidth : -childHalfWidth;
        double yOff = (q & 2) ? -childHalfWidth : childHalfWidth;
        childCenter[0] = node->center[0] + xOff;
        childCenter[1] = node->center[1] + yOff;
        node->children[q] = allocateNode(childCenter, childHalfWidth);
    }
    insert(node->children[q], particle, depth + 1);
}

void QuadTree::buildParallel(std::vector<Particle>& particles, double center[2], double halfWidth) {
    clear();
    node_pool.resize(std::max((size_t)(6 * particles.size()), (size_t)100));
    root = allocateNode(center, halfWidth);

    int n = particles.size();
    if (n == 0) return;

    // Pre-allocate partition buffer ONCE — reused at every level.
    // This eliminates the malloc/free that was happening at every
    // partition_by_quadrant call. The buffer is N particles large,
    // which is always enough since the root-level partition is the
    // largest and handles all N particles.
    partitionBuffer.resize(n);

    buildSubtreeParallel(root, particles.data(), n, 0, partitionBuffer.data());

    computeMassDistribution(root);
}

/**
 * 
 */
void QuadTree::partitionByQuadrant(
    Particle* particles,
    int count,
    const double center[2],
    int quadrantCounts[4],
    Particle* tempBuffer) 
    {
        for (int i{0}; i < 4; i++) {
            quadrantCounts[i] = 0;
        }

        for (int i = 0; i < count; i++) {
            int q = getQuadrant(particles[i].pos, center);
            quadrantCounts[q]++;
        }

        // Start offsets for each quadrant
        int offsets[4];
        offsets[0] = 0;
        offsets[1] = quadrantCounts[0];
        offsets[2] = offsets[1] + quadrantCounts[1];
        offsets[3] = offsets[2] + quadrantCounts[2];

        int pos[4] = {offsets[0], offsets[1], offsets[2], offsets[3]};

        // Insert into temp buffer continuously starting from 0
        for (int i = 0; i < count; i++) {
            int q = getQuadrant(particles[i].pos, center);
            tempBuffer[pos[q]++] = particles[i];
        }

        // Swap back
        std::copy(tempBuffer, tempBuffer + count, particles);
    }

void QuadTree::buildSubtreeParallel(QuadTreeNode* node, Particle* particles, int count, int depth, Particle* tempBuffer) {
    if (count == 0) return;

    if (count == 1) {
        node->particle = &particles[0];
        node->isLeaf = true;
        return;
    }

    if (depth >= Constants::PARALLEL_TREE_DEPTH) {
        buildSubtreeSequential(node, particles, count, depth);
        return;
    }

    // Partition particles into four quadrants using pre-allocated buffer
    node->isLeaf = false;
    int quadrant_counts[4];
    partitionByQuadrant(particles, count, node->center, quadrant_counts, tempBuffer);

    double childHalfWidth = node->halfWidth / 2.0;
    int offset = 0;

    for (int q = 0; q < 4; q++) {
        if (quadrant_counts[q] > 0) {
            double childCenter[2];
            double xOff = (q & 1) ? childHalfWidth : -childHalfWidth;
            double yOff = (q & 2) ? -childHalfWidth : childHalfWidth;
            childCenter[0] = node->center[0] + xOff;
            childCenter[1] = node->center[1] + yOff;
            node->children[q] = allocateNode(childCenter, childHalfWidth);

            QuadTreeNode* child = node->children[q];
            Particle* child_particles = particles + offset;
            int child_count = quadrant_counts[q];
            Particle* child_temp = tempBuffer + offset;

            #pragma omp task shared(node_pool) firstprivate(child, child_particles, child_count, depth, child_temp)
            {
                buildSubtreeParallel(child, child_particles, child_count, depth + 1, child_temp);
            }
        }
        offset += quadrant_counts[q];
    }

    #pragma omp taskwait
}

/**
 * If parallelism no longer benefits while building the tree, just use sequential
 */
void QuadTree::buildSubtreeSequential(QuadTreeNode* node, Particle* particles,
                                         int count, int depth) {
    for (int i = 0; i < count; i++) {
        insert(node, &particles[i], depth);
    }
}

/**
 * Once we finish building the tree based on position on particles,
 * calculate all mass distributions bottom up, i.e., traverse all 
 * the way to the leaf, use it as a single particle cluster, then
 * move back up, at each step, we take the sum of masses of all 
 * four subtrees and average of COMs
 */
void QuadTree::computeMassDistribution(QuadTreeNode* node) {
    if (node == nullptr) return;

    if (node->isLeaf) {
        if (node->particle != nullptr) {
            node->totalMass = node->particle->mass;
            node->com[0] = node->particle->pos[0];
            node->com[1] = node->particle->pos[1];
        }
        return;
    }

    node->totalMass = 0.0;
    node->com[0] = 0.0;
    node->com[1] = 0.0;

    for (int i = 0; i < 4; i++) {
        if (node->children[i] != nullptr) {
            computeMassDistribution(node->children[i]);
            double child_mass = node->children[i]->totalMass;
            node->totalMass += child_mass;
            node->com[0] += child_mass * node->children[i]->com[0];
            node->com[1] += child_mass * node->children[i]->com[1];
        }
    }

    if (node->totalMass > 0.0) {
        node->com[0] /= node->totalMass;
        node->com[1] /= node->totalMass;
    }
}


void QuadTree::computeForce(Particle& target) const {
    target.acc[0] = 0.0;
    target.acc[1] = 0.0;

    if (root != nullptr) {
        computeForceRecursive(root, target);
    }
}

/**
 * Compute the force on a certain particle based on the Barnes-Hut approach.
 * If a particle is 'far enough' from a cluster of particles represented by a node
 * by checking if s/d < theta, then we can compute the force due to the cluster
 * by assuming a single particle of total mass equal to all particles in the cluster,
 * centred at the cluster COM
 */
void QuadTree::computeForceRecursive(const QuadTreeNode* node, Particle& target) const {
    if (node == nullptr) {
        return;
    }

    if (node->isLeaf && node->particle == nullptr) {
        return;
    }

    // If we have reached the leaf, we have maximum granuality, calculate forces normally
    if (node->isLeaf && node->particle != nullptr) {
        if (node->particle == &target) {
            return;
        }

        double dx = node->particle->pos[0] - target.pos[0];
        double dy = node->particle->pos[1] - target.pos[1];
        double distSq = dx * dx + dy * dy + Constants::SOFTENING * Constants::SOFTENING;
        double invDist = 1.0 / std::sqrt(distSq);
        double invDistCubed = invDist * invDist * invDist;

        target.acc[0] += Constants::G * node->particle->mass * dx * invDistCubed;
        target.acc[1] += Constants::G * node->particle->mass * dy * invDistCubed;
        return;
    }

    double dx = node->com[0] - target.pos[0];
    double dy = node->com[1] - target.pos[1];
    double distSq = dx * dx + dy * dy;
    double cellSize = 2.0 * node->halfWidth;
    
    // If (s / d) < theta, do not recurse further
    if (cellSize * cellSize < this->theta * this->theta * distSq) {
        double softDistSq = distSq + Constants::SOFTENING * Constants::SOFTENING;
        double invDist = 1.0 / std::sqrt(softDistSq);
        double invDistCubed = invDist * invDist * invDist;

        target.acc[0] += Constants::G * node->totalMass * dx * invDistCubed;
        target.acc[1] += Constants::G * node->totalMass * dy * invDistCubed;
    } else {
        for (int i = 0; i < 4; i++) {
            if (node->children[i] != nullptr) {
                computeForceRecursive(node->children[i], target);
            }
        }
    }
}

/**
 * Compute the total force on each particle brute force, i.e., FROM each particle TO each particle
 */
void QuadTree::computeForceBrute(Particle& target, const std::vector<Particle>& particles) {
    target.acc[0] = 0.0;
    target.acc[1] = 0.0;

    for (const auto& p : particles) {
        if (&p == &target) continue;

        double dx = p.pos[0] - target.pos[0];
        double dy = p.pos[1] - target.pos[1];
        double distSq = dx * dx + dy * dy + Constants::SOFTENING * Constants::SOFTENING;
        double invDist = 1.0 / std::sqrt(distSq);
        double invDistCubed = invDist * invDist * invDist;

        target.acc[0] += Constants::G * p.mass * dx * invDistCubed;
        target.acc[1] += Constants::G * p.mass * dy * invDistCubed;
    }
}

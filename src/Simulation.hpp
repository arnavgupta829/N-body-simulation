#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "StructDefs.hpp"
#include "QuadTree.hpp"
#include <vector>

class Simulation {
  public:
    struct TimingData {
      double total = 0.0;
      double tree_build = 0.0;       // includes Morton sort
      double morton_sort = 0.0;      // Morton sort only (subset of tree_build)
      double force_compute = 0.0;    // travel through tree, compute forces
      double integration = 0.0;      // vel + pos updates
      double energy = 0.0;           // diagnostic energy computation (to verify correctness)
      double domain_bounds = 0.0;    // bounding box computation
    };

    Config config;
    
    // Simulation stats
    TimingData timing{};

    // B-H quadtree
    QuadTree tree;

    std::vector<Particle> particles;

    /**
     * For morton key sorting, we want to avoid repeatedly allocating memory
     * for each tree recalculation. We define all temp vectors we need once
     * during initialization and pass the same objects along each time
     */
    std::vector<Particle> sorted;
    std::vector<uint64_t> keys;
    std::vector<int> indices;

    // For convenience, storing whether our execution is sequential or parallel separately
    bool isParallel;

    virtual ~Simulation() = default;

    // Common methods
    void initialize(const Config& config);
    void initPlummer(int n);
    void initUniform(int n);
    void computeDomainBounds(double center[2], double& halfWidth);

    // Methods to be implemented by child classes
    virtual void buildTree(double center[2], double halfWidth) = 0;
    virtual void computeForce() = 0;
    virtual void stepVelocity(double dt) = 0;
    virtual void stepPosition(double dt) = 0;
    virtual void step(double dt) = 0;
    virtual double computeKineticEnergy() = 0;
    virtual double computerPotentialEnergy() = 0;
    virtual void run() = 0;
};

#endif
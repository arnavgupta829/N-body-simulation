#ifndef PARALLEL_SIMULATION_HPP
#define PARALLEL_SIMULATION_HPP

#include "Simulation.hpp"

class ParallelSimulation : public Simulation {
  public:
    void buildTree(double center[2], double halfWidth) override;
    void computeForce() override;
    void stepVelocity(double dt) override;
    void stepPosition(double dt) override;
    void step(double dt) override;
    double computeKineticEnergy() override;
    double computerPotentialEnergy() override;
    void run() override;
};

#endif
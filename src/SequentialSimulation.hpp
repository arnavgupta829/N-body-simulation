#ifndef SEQUENTIAL_SIMULATION_HPP
#define SEQUENTIAL_SIMULATION_HPP

#include "Simulation.hpp"

class SequentialSimulation : public Simulation {
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
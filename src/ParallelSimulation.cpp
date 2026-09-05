#include "Constants.hpp"
#include "StructDefs.hpp"
#include "Morton.hpp"
#include "QuadTree.hpp"
#include "ParallelSimulation.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <omp.h>

void ParallelSimulation::buildTree(double center[2], double halfWidth) {
  if (config.treeBuilder == PARALLEL) {
    double sort_t0 = omp_get_wtime();
    Morton::sortParticles(particles, center, halfWidth, keys, sorted, indices);
    double sort_t1 = omp_get_wtime();
    timing.morton_sort += (sort_t1 - sort_t0);
    tree.buildParallel(particles, center, halfWidth);
  } else {
    tree.build(particles, center, halfWidth);
  }
}

void ParallelSimulation::computeForce() {
  #pragma omp for schedule(dynamic, 64) 
  for (int i = 0; i < this->config.numParticles; i++) {
    tree.computeForce(particles[i]);
  }
}

void ParallelSimulation::stepVelocity(double dt) {
  #pragma omp for schedule(static) 
  for (int i = 0; i < this->config.numParticles; i++) {
      particles[i].vel[0] += 0.5 * dt * particles[i].acc[0];
      particles[i].vel[1] += 0.5 * dt * particles[i].acc[1];
  }
}

void ParallelSimulation::stepPosition(double dt) {
  #pragma omp for schedule(static) 
  for (int i = 0; i < this->config.numParticles; i++) {
      particles[i].pos[0] += dt * particles[i].vel[0];
      particles[i].pos[1] += dt * particles[i].vel[1];
  }
}

void ParallelSimulation::step(double dt) {
  double t0, t1;
  double center[2];
  double halfWidth;

  #pragma omp parallel 
  {

    #pragma omp single
    {
      t0 = omp_get_wtime();
    }
    /**
     * Step 1 - Calculate new velocity at half timestep based on initial 
     *          position of particles
     */
    stepVelocity(dt);

    /**
     * Step 2 - Calculate new position at full timestep based on velocity
     *          at half timestep
     */
    stepPosition(dt);

    #pragma omp single 
    {
      t1 = omp_get_wtime();
      timing.integration += (t1 - t0);
    }

    /**
     * Step 3 - Compute new root level bounding box. Recalculate tree and
     *          Calculate new forces based on position on particles at 
     *          full timestep. From our experiments we found that
     *          computing domain bounds does not need to be parallelized
     */
    #pragma omp single 
    {
      t0 = omp_get_wtime();
      computeDomainBounds(center, halfWidth);
      t1 = omp_get_wtime();
      timing.domain_bounds += (t1 - t0);
    }
    
    #pragma omp single 
    {
      t0 = omp_get_wtime();
      buildTree(center, halfWidth);
      t1 = omp_get_wtime();
      timing.tree_build += (t1 - t0);
    }

    /**
     * Step 4 - Recomputer new forces based on the new tree
     */
    #pragma omp single
    {
      t0 = omp_get_wtime();
    }
    computeForce();
    #pragma omp single
    {
      t1 = omp_get_wtime();
      timing.force_compute += (t1 - t0);
    }

    /**
     * Step 5 - Calculate velocity at full timestep using new forces calculated
     *          for full timestep
     */
    #pragma omp single
    {
      t0 = omp_get_wtime();
    }
    stepVelocity(dt);
    #pragma omp single
    {
      t1 = omp_get_wtime();
      timing.integration += (t1 - t0);
    }

  }
}

double ParallelSimulation::computeKineticEnergy() {
  int n{config.numParticles};
  double ke{0.0};

  #pragma omp parallel for reduction(+:ke) schedule(static)
  for (int i = 0; i < n; i++) {
      ke += 0.5 * particles[i].mass 
        * (particles[i].vel[0] * particles[i].vel[0]
          + particles[i].vel[1] * particles[i].vel[1]);
  }
  
  return ke;
}

double ParallelSimulation::computerPotentialEnergy() {
  int n{config.numParticles};
  double pe{0.0};

  #pragma omp parallel for reduction(+:pe) schedule(dynamic, 16) 
  for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
          double dx = particles[j].pos[0] - particles[i].pos[0];
          double dy = particles[j].pos[1] - particles[i].pos[1];
          double dist = std::sqrt(dx * dx + dy * dy + Constants::SOFTENING * Constants::SOFTENING);
          pe -= Constants::G * particles[i].mass * particles[j].mass / dist;
      }
  }

  return pe;
}

void ParallelSimulation::run() {
  double run_start = omp_get_wtime();

  double center[2];
  double half_width;
  computeDomainBounds(center, half_width);

  #pragma omp parallel
  {
    #pragma omp single 
    {
      buildTree(center, half_width);
    }
    computeForce();
  }

  /**
   * Verbose tools
   */
  std::cout << std::setw(8) << "Step"
            << std::setw(16) << "KE"
            << std::setw(16) << "PE"
            << std::setw(16) << "Total E"
            << std::setw(16) << "dE/E0 (%)"
            << "\n";
  std::cout << std::string(72, '-') << "\n";

  double t_energy_start = omp_get_wtime();
  // double KE = computeKineticEnergy();
  // double PE = computerPotentialEnergy();
  // double E0 = KE + PE;
  double E0{0.0};
  timing.energy += omp_get_wtime() - t_energy_start;

  // std::cout << std::setw(8) << 0
  //           << std::setw(16) << std::scientific << std::setprecision(6) << KE
  //           << std::setw(16) << PE
  //           << std::setw(16) << E0
  //           << std::setw(16) << 0.0
  //           << "\n";
  std::cout << 0 << "\n";

  for (int s = 1; s <= config.numSteps; s++) {
  // Take step
  step(Constants::DT);

  // Print debug information at every 10th step
  if (s % 10 == 0 || s == config.numSteps) {
    double t_e0 = omp_get_wtime();
    // KE = computeKineticEnergy();
    // PE = computerPotentialEnergy();
    // double E = KE + PE;
    double E{0.0};
    timing.energy += omp_get_wtime() - t_e0;

    // double dE_rel = (E - E0) / std::abs(E0) * 100.0;
    // std::cout << std::setw(8) << s
    //           << std::setw(16) << std::scientific << std::setprecision(6) << KE
    //           << std::setw(16) << PE
    //           << std::setw(16) << E
    //           << std::setw(16) << std::fixed << std::setprecision(6) << dE_rel
    //           << "\n";
    std::cout << s << "\n";
    }
  }

  timing.total = omp_get_wtime() - run_start;
}
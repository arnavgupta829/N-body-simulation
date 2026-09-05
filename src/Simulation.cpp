#define _USE_MATH_DEFINES

#include "Constants.hpp"
#include "StructDefs.hpp"
#include "Morton.hpp"
#include "QuadTree.hpp"
#include "Simulation.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <omp.h>

void Simulation::initialize(const Config& config) {
  this->config = config;
  isParallel = config.isParallel;

  int numParticles{config.numParticles};
  particles.resize(numParticles);
  sorted.resize(numParticles);
  keys.resize(numParticles);
  indices.resize(numParticles);

  tree.theta = config.theta;

  switch (this->config.loaderStyle) {
      case PLUMMER:
          initPlummer(this->config.numParticles);
          break;

      case UNIFORM:
      default:
          initUniform(this->config.numParticles);
          break;
  }
}

void Simulation::initPlummer(int n) {
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  double plummerRadius{1.0};
  double totalMass{1.0};
  double particleMass{totalMass / n};

  for (int i{0}; i < n; i++) {
      particles[i].mass = particleMass;

      double u{uniform(rng)};
      u = std::max(1e-10, std::min(u, 1.0 - 1e-10));
      double r{plummerRadius / std::sqrt(std::pow(u, -2.0 / 3.0) - 1.0)};

      double angle{2.0 * M_PI * uniform(rng)};
      particles[i].pos[0] = r * std::cos(angle);
      particles[i].pos[1] = r * std::sin(angle);

      double vEsc{std::sqrt(2.0 * Constants::G * totalMass / std::sqrt(r * r + plummerRadius * plummerRadius))};

      double xSample, ySample;
      do {
          xSample = uniform(rng);
          ySample = uniform(rng) * 0.1;
      } while (ySample > xSample * xSample * std::pow(1.0 - xSample * xSample, 3.5));

      double v{xSample * vEsc};
      double vAngle = 2.0 * M_PI * uniform(rng);
      particles[i].vel[0] = v * std::cos(vAngle);
      particles[i].vel[1] = v * std::sin(vAngle);

      particles[i].acc[0] = 0.0;
      particles[i].acc[1] = 0.0;
  }

  double vxAvg{0.0}, vyAvg{0.0};
  for (const auto& p : particles) {
      vxAvg += p.mass * p.vel[0];
      vyAvg += p.mass * p.vel[1];
  }
  vxAvg /= totalMass;
  vyAvg /= totalMass;
  for (auto& p : particles) {
      p.vel[0] -= vxAvg;
      p.vel[1] -= vyAvg;
  }
}

void Simulation::initUniform(int n) {
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);

  double totalMass = 1.0;
  double particleMass = totalMass / n;
  double radius = 10.0;

  for (int i = 0; i < n; i++) {
      particles[i].mass = particleMass;
      double x, y;
      do {
          x = uniform(rng) * radius;
          y = uniform(rng) * radius;
      } while (x * x + y * y > radius * radius);
      particles[i].pos[0] = x;
      particles[i].pos[1] = y;
      particles[i].vel[0] = uniform(rng) * 0.01;
      particles[i].vel[1] = uniform(rng) * 0.01;
      particles[i].acc[0] = 0.0;
      particles[i].acc[1] = 0.0;
  }
}

void Simulation::computeDomainBounds(double center[2], double& halfWidth) {
  int n = particles.size();

  double xMin = particles[0].pos[0];
  double xMax = particles[0].pos[0];
  double yMin = particles[0].pos[1];
  double yMax = particles[0].pos[1];

  for (int i = 0; i < n; i++) {
      xMin = std::min(xMin, particles[i].pos[0]);
      xMax = std::max(xMax, particles[i].pos[0]);
      yMin = std::min(yMin, particles[i].pos[1]);
      yMax = std::max(yMax, particles[i].pos[1]);
  }

  center[0] = (xMin + xMax) / 2.0;
  center[1] = (yMin + yMax) / 2.0;

  double dx = (xMax - xMin) / 2.0;
  double dy = (yMax - yMin) / 2.0;
  // 1% margin to avoid boundary issues
  halfWidth = std::max(dx, dy) * 1.01;
}
#include "Simulation.hpp"
#include "ParallelSimulation.hpp"
#include "SequentialSimulation.hpp"
#include "Constants.hpp"
#include "StructDefs.hpp"
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <omp.h>

int main(int argc, char* argv[]) {
    Config config{};

    std::string tempArg;

    int opt;
    while ((opt = getopt(argc, argv, "hn:s:p:t:m:T:")) != -1) {
        switch(opt) {
            case 'h':
                printf("Usage: %s [-n<numparticles>] [-s<numsteps>] [-p<particleloaderstyle>] [-t<numthreads>] [-m<treebuilderstyle>] [-T<BHTheta>]\n", argv[0]);
                printf("\t-n<numparticles> number of particles for the simulation\n");
                printf("\t-s<numsteps> total number of time-steps for the simulation\n");
                printf("\t-p<particleloaderstyle> particles loaded with uniform or plummer distribution\n");
                printf("\t-t<numthreads> number of threads (defaults to sequential)\n");
                printf("\t-m<treebuilderstyle> Barnes-Hut tree built parallelly or sequentially\n");
                printf("\t-T<BHTheta> parameter theta for the simulation\n");
                std::exit(0);
                break;
            
            case 'n':
                config.numParticles = std::atoi(optarg);
                break;

            case 's':
                config.numSteps = std::atoi(optarg);
                break;

            case 'p':
                tempArg = optarg;
                config.loaderStyle = (tempArg == "p" || tempArg == "P") ? PLUMMER : UNIFORM;
                break;

            case 't':
                config.numThreads = std::atoi(optarg);
                if (config.numThreads > 1) {
                    config.isParallel = true;
                }
                break;

            case 'm':
                tempArg = optarg;
                config.treeBuilder = (tempArg == "p" || tempArg == "P") ? PARALLEL : SEQUENTIAL;
                break;

            case 'T':
                config.theta = std::atof(optarg);
                break;
        }
    }

    // Set thread count
    if (config.isParallel) {
        omp_set_num_threads(config.numThreads);
    } else {
        omp_set_num_threads(1);
    }

    std::cout << "### Barnes-Hut N-Body Simulation ###\n";
    std::cout << "Particles:    " << config.numParticles << "\n";
    std::cout << "Steps:        " << config.numSteps << "\n";
    std::cout << "Mode:         " << (config.isParallel ? "PARALLEL" : "SEQUENTIAL") << "\n";
    std::cout << "Threads:      " << config.numThreads << "\n";
    std::cout << "Particle Load " << ((config.loaderStyle == PLUMMER) ? "PLUMMER" : "UNIFORM") << "\n";
    std::cout << "Tree Builder  " << ((config.treeBuilder == PARALLEL) ? "PARALLEL" : "SEQUENATIAL") << "\n";
    std::cout << "Theta:        " << config.theta << "\n";
    std::cout << "dt:           " << Constants::DT << "\n";
    std::cout << "Softening:    " << Constants::SOFTENING << "\n";
    std::cout << "\n";

    // Run simulation
    Simulation* sim = nullptr;
    if (config.isParallel) {
        sim = new ParallelSimulation();
    } else {
        sim = new SequentialSimulation();
    }
    sim->initialize(config);
    sim->run();

    // Overall stats
    const auto& t = sim->timing;
    std::cout << "\n";
    std::cout << "### Timing Breakdown ###" << "\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total time:       " << t.total << " s\n";
    std::cout << "  Tree build:     " << t.tree_build << " s (" << (t.tree_build / t.total * 100.0) << "%)\n";
    if (config.treeBuilder == PARALLEL) {
        std::cout << "    Morton sort:  " << t.morton_sort << " s (" << (t.morton_sort / t.total * 100.0) << "%)\n";
    }
    std::cout << "  Force compute:  " << t.force_compute << " s (" << (t.force_compute / t.total * 100.0) << "%)\n";
    std::cout << "  Integration:    " << t.integration << " s (" << (t.integration / t.total * 100.0) << "%)\n";
    std::cout << "  Domain bounds:  " << t.domain_bounds << " s (" << (t.domain_bounds / t.total * 100.0) << "%)\n";
    std::cout << "  Energy diag:    " << t.energy << " s (" << (t.energy / t.total * 100.0) << "%)\n";
    std::cout << "\n";
    std::cout << "### Performance ###\n";
    std::cout << "Time per step:     " << t.total / config.numSteps << " s\n";

    return 0;
}

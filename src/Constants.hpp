#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace Constants {

    constexpr double G = 1.0;
    constexpr double DEFAULT_THETA = 0.5;
    constexpr double SOFTENING = 0.025;
    constexpr double DT = 0.005;
    constexpr int NUM_STEPS = 100;
    constexpr int DEFAULT_NUM_PARTICLES = 1000;
    constexpr double DOMAIN_SIZE = 100.0;
    constexpr int MAX_TREE_DEPTH = 64;
    constexpr int PARALLEL_TREE_DEPTH = 4;
    
}

#endif

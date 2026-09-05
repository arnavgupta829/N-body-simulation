#ifndef MORTON_HPP
#define MORTON_HPP

#include "StructDefs.hpp"
#include <cstdint>
#include <vector>

namespace Morton {
    uint64_t computeKey(
        double px, 
        double py,           
        double minXDomain, 
        double minYDomain,
        double domainSize, 
        int bits = 21);

    uint64_t interleaveBits(uint32_t x, uint32_t y);

    void sortParticles(
        std::vector<Particle>& particles,
        const double domainCenter[2],
        double domainHalfWidth,
        std::vector<uint64_t>& keys,
        std::vector<Particle>& sorted,
        std::vector<int>& indices);

}

#endif

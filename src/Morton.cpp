#include "Morton.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>

namespace Morton {

/**
 * Spread 32 bits of x so that each bit has a gap of 1 bit in between
 */
static uint64_t spreadBits(uint32_t x) {
    uint64_t y = x;
    y = (y | (y << 16)) & 0x0000FFFF0000FFFF;
    y = (y | (y <<  8)) & 0x00FF00FF00FF00FF;
    y = (y | (y <<  4)) & 0x0F0F0F0F0F0F0F0F;
    y = (y | (y <<  2)) & 0x3333333333333333;
    y = (y | (y <<  1)) & 0x5555555555555555;
    return y;
}

uint64_t interleaveBits(uint32_t x, uint32_t y) {
    // x goes into even bits (0, 2, 4, ...), y goes into odd bits (1, 3, 5, ...)
    return spreadBits(x) | (spreadBits(y) << 1);
}

uint64_t computeKey(
    double px, 
    double py,                 
    double domainMinX, 
    double domainMinY,
    double domainSize, 
    int bits) {
       
        // Normalize position
        double nx = (px - domainMinX) / domainSize;
        double ny = (py - domainMinY) / domainSize;

        // Handle edge case, subtract delta
        nx = std::max(0.0, std::min(nx, 1.0 - 1e-10));
        ny = std::max(0.0, std::min(ny, 1.0 - 1e-10));

        uint32_t max_val = (1u << bits) - 1;
        uint32_t ix = static_cast<uint32_t>(nx * (1u << bits));
        uint32_t iy = static_cast<uint32_t>(ny * (1u << bits));

        // Clamp to valid range
        ix = std::min(ix, max_val);
        iy = std::min(iy, max_val);

        return interleaveBits(ix, iy);
    }

/**
 * Design choice: during the simulations, we realized that several temp vectors
 * were being created inside this function. TO try and save some time we pre-
 * allocated the temp vectors needed for sorting, but turns out sorting particles
 * for N < 1,000,000 parallely does not seem to be very useful anyways (in fact
 * it seems to worsen performance as shown in the report). In that case
 * we do not think it REALLY makes a difference since the malloc would be outside
 * of any parallel region.
 */
void sortParticles(
    std::vector<Particle>& particles, 
    const double domainCenter[2],  
    double domainHalfWidth, 
    std::vector<uint64_t>& keys, 
    std::vector<Particle>& sorted, 
    std::vector<int>& indices) {
        int n = particles.size();
        if (n <= 1) return;

        // Domain bounds
        double domain_min_x = domainCenter[0] - domainHalfWidth;
        double domain_min_y = domainCenter[1] - domainHalfWidth;
        double domain_size = 2.0 * domainHalfWidth;

        // Compute morton keys
        for (int i = 0; i < n; i++) {
            keys[i] = computeKey(particles[i].pos[0], particles[i].pos[1], domain_min_x, domain_min_y, domain_size);
            indices[i] = i;
        }

        // Sort indices by their Morton key
        std::sort(indices.begin(), indices.end(), [&keys](int a, int b) { return keys[a] < keys[b]; });

        // Reorder particles according to sorted indices
        for (int i = 0; i < n; i++) {
            sorted[i] = particles[indices[i]];
        }

        // Swap vectors
        particles.swap(sorted);
    }
}

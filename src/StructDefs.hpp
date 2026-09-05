#ifndef STRUCT_DEFS_HPP
#define STRUCT_DEFS_HPP

/**
 * Enums/Structs defined for the configuration object 
 */
enum LoaderStyle {
  UNIFORM,
  PLUMMER
};

enum TreeBuilder {
  SEQUENTIAL,
  PARALLEL
};

struct Config {
  int numParticles{10000};
  int numSteps{100};
  
  LoaderStyle loaderStyle{UNIFORM};

  bool isParallel{false};
  int numThreads{1};
  
  TreeBuilder treeBuilder{SEQUENTIAL};

  double theta{0.5};
};

/**
 * Particle Struct
 */
struct Particle {
    double mass;
    double pos[2];
    double vel[2];
    double acc[2];
};

/**
 * QuadTree Struct - One cell of the BH tree. Can be subdivided further
 */
struct QuadTreeNode {
    double totalMass;
    double com[2];

    double center[2];
    double halfWidth;

    QuadTreeNode* children[4];
    Particle* particle;
    bool isLeaf;
};

#endif
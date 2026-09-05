Running the code - 

1. Copy the code into any of the crunchy machines (I have tested on crunchy3 and 5)
  i.  The makefile is in the project root along with this readme
  ii. The code is inside the `src' dir

2. use the command make to create the executable `nbody'

```
[ag11163@crunchy3 ag11163_project]$ make
```

3. The following args are used to set up different parameters for the simulation
  i.    -n<number of particles> 
  ii.   -s<number of steps>
  iii.  -p<particle initialization - either `U' or `P' for uniform or plummer respectively>
  iv.   -t<number of threads. passing in `1' here makes the simulation sequential automatically>
  v.    -m<build tree parallelly using `P' or sequentially using `S'>
  vi.   -T<value > 0.1 which is used as the barnes hut relaxation parameter theta>

4. Example of a valid command being run

```
[ag11163@crunchy3 ag11163]$ ./nbody -n100000 -s100 -pP -t64 -mP -T0.5
```

Here we have
  i.    initialized with 100,000 particles
  ii.   number of iterations is 100
  iii.  used the plummer distribution to initalize the particles
  iv.   number of threads is 64
  v.    build tree parallelly, i.e., work done to create tree is spread across threads
  vi.   barnes hut parameter theta is 0.5

5. The simulation shows you the following information
  i.    initial parameters we have launched with
  ii.   debugging stats at the end of every tenth iteration. We verify the simulation
        is correct by noting that energy is conserved to within 0.01% throughout
  iii.  final stats which include time taken for each part of the code (sequential AND parallel)
    a.  building the tree
    b.  computing forces
    c.  half and full step integrations to solve for positions and velocities
    d.  energy calculations
    e.  recomputation of the bounding box of the simulation whenever the particles move

Here is a sample output

```
### Barnes-Hut N-Body Simulation ###
Particles:    100000
Steps:        100
Mode:         PARALLEL
Threads:      64
Particle Load PLUMMER
Tree Builder  SEQUENATIAL
Theta:        0.5
dt:           0.005
Softening:    0.025

    Step              KE              PE         Total E       dE/E0 (%)
------------------------------------------------------------------------
       0    1.470494e-01   -3.477144e-01   -2.006649e-01    0.000000e+00
      10    1.470602e-01   -3.477251e-01   -2.006649e-01        0.000002
      20    1.470903e-01   -3.477556e-01   -2.006653e-01       -0.000199
      30    1.471428e-01   -3.478089e-01   -2.006662e-01       -0.000628
      40    1.472244e-01   -3.478920e-01   -2.006676e-01       -0.001324
      50    1.473344e-01   -3.480038e-01   -2.006694e-01       -0.002252
      60    1.474804e-01   -3.481522e-01   -2.006719e-01       -0.003473
      70    1.476584e-01   -3.483332e-01   -2.006747e-01       -0.004902
      80    1.478716e-01   -3.485499e-01   -2.006783e-01       -0.006657
      90    1.481226e-01   -3.488049e-01   -2.006823e-01       -0.008687
     100    1.484064e-01   -3.490933e-01   -2.006870e-01       -0.010991

### Timing Breakdown ###
Total time:       83.9061 s
  Tree build:     34.4988 s (41.1160%)
  Force compute:  18.4997 s (22.0481%)
  Integration:    2.2943 s (2.7343%)
  Domain bounds:  0.7308 s (0.8710%)
  Energy diag:    25.4857 s (30.3741%)

### Performance ###
Time per step:     0.8391 s
```

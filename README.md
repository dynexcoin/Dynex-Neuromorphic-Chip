# Dynex-Neuromorphic-Chip

Reference implementations of the Dynex Neuromorphic Chip as described [on our website](https://dynexcoin.org/dynex-neuromorphic-chip/).

## C++ Reference Implementation (CPU)

File "dynex.cc" is a reference implementation of the neuromorphic circuit. It integrates the system of equations (ODE) derived from the circuit model on a CPU. The reference design solves the Boolean satisfiability problem (sometimes called propositional satisfiability problem and abbreviated SATISFIABILITY, SAT or B-SAT). A propositional logic formula, also called Boolean expression, is built from variables, operators AND (conjunction, also denoted by ∧), OR (disjunction, ∨), NOT (negation, ¬), and parentheses. A formula is said to be satisfiable if it can be made TRUE by assigning appropriate logical values (i.e. TRUE, FALSE) to its variables. The Boolean satisfiability problem (SAT) is, given a formula, to check whether it is satisfiable. This decision problem is of central importance in many areas of computer science, including theoretical computer science, complexity theory, algorithmics, cryptography and artificial intelligence. Please note that this reference implementation proves the mathematical model and is not as efficient as the implementation by commerical miners. 

To illustrate the performance of neuromorphic computing, the following example showcases an implementation of a constraint satisfaction problem, where a problem formulation with complexity O(n^100,000) is being solved using the Dynex Neuromorphic Chip. The problem consists of 100,000 unique variables. Existing methodologies based on current and Quantum technology (reducing the complexity with Shor’s algorithm to O(n^50,000) cannot solve this problem class efficiently today. The Dynex Neuromorphic Chip solves the problem in a few seconds because of its inherent parallelization, it’s long-range order and its capability to utilize instantons.

This reference implementation proves the mathematical model and has not been optimised for speed and performance.

### Requirements:
Please note that it is required to have the [Boost library](https://www.boost.org) (Version 1.74.0 or better) installed: 

```
sudo apt-get install libboost-all-dev (Ubuntu Linux)
brew install boost (MacOS)
```

### Build from source:

```
g++ dynex.cc -o dynex -std=c++17 -Ofast -lpthread -fpermissive (Linux)
g++ dynex.cc -o dynex -std=c++17 -Ofast -I /opt/homebrew/cellar/boost/1.78.0/include -L /opt/homebrew/cellar/boost/1.78.0/lib (MacOS)
```

### Run the neuromorphic sat solver:

```
./dynex -i cnf/transformed_barthel_n_100_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_200_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_500_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_1000_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_10000_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_002.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_004.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_005.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_007.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_008.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_014.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_016.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_018.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_020.cnf
./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_024.cnf
```

### Benchmark Comparison 

Comparison with Kissat_MAB-HyWalk, winner of the SAT Competition 2023 [http://www.satcompetition.org/], a single threaded solver.

```
PROBLEM INSTANCE                                                  Max.Walltime     Kissat_MAB-HyWalk         Dynex 1 Core
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf    15 minutes       no solution          10.67s (98 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_002.cnf    15 minutes       no solution          13.90s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_004.cnf    15 minutes       no solution          17.00s (159 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_005.cnf    15 minutes       no solution          15.29s (142 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_007.cnf    15 minutes       no solution          16.61s (149 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_008.cnf    15 minutes       no solution           0.96s (18 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_014.cnf    15 minutes       no solution          15.32s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_016.cnf    15 minutes       no solution          15.05s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_018.cnf    15 minutes       no solution           4.66s (57 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_020.cnf    15 minutes       no solution          17.20s (166 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_024.cnf    15 minutes       no solution           6.22s (62 steps*)
```

Comparison with YalSat, winner of the SAT Competition 2017 Random Track [https://github.com/arminbiere/yalsat], a single threaded solver.

```
PROBLEM INSTANCE                                                  Max.Walltime     YalSat 1.0.1         Dynex 1 Core
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf    15 minutes       14.1s (20M flips)    10.67s (98 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_002.cnf    15 minutes       no solution          13.90s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_004.cnf    15 minutes       no solution          17.00s (159 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_005.cnf    15 minutes       no solution          15.29s (142 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_007.cnf    15 minutes       no solution          16.61s (149 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_008.cnf    15 minutes       no solution           0.96s (18 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_014.cnf    15 minutes       no solution          15.32s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_016.cnf    15 minutes       no solution          15.05s (138 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_018.cnf    15 minutes       no solution           4.66s (57 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_020.cnf    15 minutes       62.36s (104M flips)  17.20s (166 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_024.cnf    15 minutes       no solution           6.22s (62 steps*)
```

Comparison with PalSat [https://github.com/arminbiere/yalsat], a multii threaded solver, on 3 CPU cores (command line option -t 3) and Dynex on 3 CPU cores (command line option -w 3)

```
PROBLEM INSTANCE                                                  Max.Walltime     PalSat 3 Cores      Dynex 3 Cores
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf    15 minutes       50.43s (16M flips)   14.83s (98 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_002.cnf    15 minutes       no solution           1.08s (14 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_004.cnf    15 minutes       no solution           5.83s (52 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_005.cnf    15 minutes       no solution           7.72s (63 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_007.cnf    15 minutes       no solution           6.09s (57 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_008.cnf    15 minutes       28.56s (28M flips)    1.31s (18 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_014.cnf    15 minutes       no solution           6.08s (56 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_016.cnf    15 minutes       no solution           5.11s (48 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_018.cnf    15 minutes       no solution           6.07s (57 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_020.cnf    15 minutes       21.35s (20M flips)    6.60s (69 steps*)
transformed_barthel_n_100000_r_8.000_p0_0.080_instance_024.cnf    15 minutes       no solution           8.36s (62 steps*)
```

PalSat is a very efficient parallel implementation and has been designed to achieve high computational speed. It requires PalSat 16-28 million flips (=computing operations) for finding solutions. 

*In comparison, the Dynex reference implementation, which has not been optimised for speed but to showcase the calculations in a clear and understandable way in the source code, requires only between 15 - 183 integration steps to find a solution. TTS for Dynex is greatly improved in the DNX mining software implementations.

### Notes and Remarks

- The Dynex Chip is not a physical produced chip, DynexSolve simulates them by integrating its equations of motion (ODE integration)
Extrapolating the equations of motions from the chip design requires a number of variables / parameters to be introduced (as explained in the DynexSolve paper)

- In the case of SAT problems, there are in total 6 parameters to be defined. These parameters depend on the underlying structure for a problem type. The reference implementation used in our example is applying tuned parameters for the problem type stemming from Barthel instances and may not be suitable / applicable out of the box for other problem types

- Identifying the optimal parameters is an iterative tuning process which happens on the Dynex platform where a large number of GPUs are tuning in until they have identified the best values

- With parameters tuned, TTS is optimised and applicable for every similar problem types, meaning problems of similar type can be solved even more efficient on the platform.

## CUDA Reference Implementation (GPU)

File "dynex.cu" is a reference implementation of the neuromorphic circuit. It integrates the system of equations (ODE) derived from the circuit model on a GPU. The reference design solves the Boolean satisfiability problem (sometimes called propositional satisfiability problem and abbreviated SATISFIABILITY, SAT or B-SAT). A propositional logic formula, also called Boolean expression, is built from variables, operators AND (conjunction, also denoted by ∧), OR (disjunction, ∨), NOT (negation, ¬), and parentheses. A formula is said to be satisfiable if it can be made TRUE by assigning appropriate logical values (i.e. TRUE, FALSE) to its variables. The Boolean satisfiability problem (SAT) is, given a formula, to check whether it is satisfiable. This decision problem is of central importance in many areas of computer science, including theoretical computer science, complexity theory, algorithmics, cryptography and artificial intelligence. Please note that this reference implementation proves the mathematical model.

### Requirements:
Please note that it is required to have CUDA as well as the [Boost library](https://www.boost.org) (Version 1.74.0 or better) installed: 

```
sudo apt-get install libboost-all-dev (Ubuntu Linux)
brew install boost (MacOS)
```

### Build from source:

```
g++ nvcc dynex.cu -o dynex_gpu -std=c++17 -O4
```

### Run the neuromorphic sat solver:

```
./dynex_gpu -i cnf/transformed_barthel_n_100_r_8.000_p0_0.080_instance_001.cnf
./dynex_gpu -i cnf/transformed_barthel_n_200_r_8.000_p0_0.080_instance_001.cnf
./dynex_gpu -i cnf/transformed_barthel_n_500_r_8.000_p0_0.080_instance_001.cnf
./dynex_gpu -i cnf/transformed_barthel_n_1000_r_8.000_p0_0.080_instance_001.cnf
./dynex_gpu -i cnf/transformed_barthel_n_10000_r_8.000_p0_0.080_instance_001.cnf
./dynex_gpu -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf
```





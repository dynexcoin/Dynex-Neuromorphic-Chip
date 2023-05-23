# Dynex-Neuromorphic-Chip

Reference implementations and Verilog code for FPGAs of the Dynex Neuromorphic Chip as described [on our website](https://dynexcoin.org/dynex-neuromorphic-chip/).

## C++ Reference Implementation (CPU)

File "dynex.cc" is a reference implementation of the neuromorphic circuit. It integrates the system of equations (ODE) derived from the circuit model on a CPU. The reference design solves the Boolean satisfiability problem (sometimes called propositional satisfiability problem and abbreviated SATISFIABILITY, SAT or B-SAT). A propositional logic formula, also called Boolean expression, is built from variables, operators AND (conjunction, also denoted by ∧), OR (disjunction, ∨), NOT (negation, ¬), and parentheses. A formula is said to be satisfiable if it can be made TRUE by assigning appropriate logical values (i.e. TRUE, FALSE) to its variables. The Boolean satisfiability problem (SAT) is, given a formula, to check whether it is satisfiable. This decision problem is of central importance in many areas of computer science, including theoretical computer science, complexity theory, algorithmics, cryptography and artificial intelligence. Please note that this reference implementation proves the mathematical model and is not as efficient as an implementation on a FPGA chip.

To illustrate the superior performance of neuromorphic computing, the following example showcases an implementation of a constraint satisfaction problem, where a problem formulation with complexity O(n^100,000) is being solved using the Dynex Neuromorphic Chip. The problem consists of 100,000 unique variables. Existing methodologies based on current and Quantum technology (reducing the complexity with Shor’s algorithm to O(n^50,000) cannot solve this problem today as it would require longer than the existence of the universe to find a solution. The Dynex Neuromorphic Chip solves the problem in 2.23s because of its inherent parallelization, it’s long-range order and its capability to utilize instantons.

```
Current             Quantum           Dynex
Method              Method            Neuromorphic Chip
-------------------------------------------------------
O(n^100,000)         O(n^50,000)        2.23s
*Longer than        *Longer than
the universe exits  the universe exists
```

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

## CUDA Reference Implementation (GPU)

File "dynex.cu" is a reference implementation of the neuromorphic circuit. It integrates the system of equations (ODE) derived from the circuit model on a GPU. The reference design solves the Boolean satisfiability problem (sometimes called propositional satisfiability problem and abbreviated SATISFIABILITY, SAT or B-SAT). A propositional logic formula, also called Boolean expression, is built from variables, operators AND (conjunction, also denoted by ∧), OR (disjunction, ∨), NOT (negation, ¬), and parentheses. A formula is said to be satisfiable if it can be made TRUE by assigning appropriate logical values (i.e. TRUE, FALSE) to its variables. The Boolean satisfiability problem (SAT) is, given a formula, to check whether it is satisfiable. This decision problem is of central importance in many areas of computer science, including theoretical computer science, complexity theory, algorithmics, cryptography and artificial intelligence. Please note that this reference implementation proves the mathematical model and is not as efficient as an implementation on a FPGA chip.

To illustrate the superior performance of neuromorphic computing, the following example showcases an implementation of a constraint satisfaction problem, where a problem formulation with complexity O(n^100,000) is being solved using the Dynex Neuromorphic Chip. The problem consists of 100,000 unique variables. Existing methodologies based on current and Quantum technology (reducing the complexity with Shor’s algorithm to O(n^50,000) cannot solve this problem today as it would require longer than the existence of the universe to find a solution. The Dynex Neuromorphic Chip solves the problem in 2.23s because of its inherent parallelization, it’s long-range order and its capability to utilize instantons.

```
Current             Quantum           Dynex
Method              Method            Neuromorphic Chip
-------------------------------------------------------
O(n^100,000)         O(n^50,000)        2.23s
*Longer than        *Longer than
the universe exits  the universe exists
```

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





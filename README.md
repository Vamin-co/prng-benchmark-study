### PRNG Project: Evaluation of Pseudo-Random Number Generators

#### Overview
This project evaluates the performance of two pseudo-random number generators (PRNGs): Linear Congruential Generator (LCG) and XORShift. Using ten standard benchmark functions, the study assesses their effectiveness for computational tasks requiring randomness.

#### Features
- Comparison of LCG and XORShift.
- Benchmarking with functions like Schwefel, Rosenbrock, Rastrigin, and Ackley.
- Statistical analysis of fitness, range, median, and runtime for each PRNG.

#### Requirements
- C Compiler (e.g., GCC)
- Make utility

#### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/Vamin-co/prng-benchmark-study
   cd prng-benchmark-study
   ```
2. Compile the code:
   ```bash
   make
   ```

#### Usage
Run the executable to evaluate PRNGs using benchmark functions:
```bash
./prng_benchmark
```

#### Results
Performance metrics, including average fitness, standard deviation, and runtime, are displayed during execution. Benchmark functions are evaluated for both PRNGs across multiple iterations.

#### Benchmark Functions
- Schwefel
- De Jong (Sphere)
- Rosenbrock
- Rastrigin
- Griewangk
- Sine Envelope Sine Wave
- Stretched V Sine Wave
- Ackley’s One
- Ackley’s Two
- Egg Holder


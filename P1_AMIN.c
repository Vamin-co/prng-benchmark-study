/*
* Name: Vandan Amin
* Student ID: 44006979
* Date: Oct 25 2024 
* Honor Code: I pledge that this submission is solely my work.
*
* CS 465 - Scientific Computing
* Department of Computer Science, Central Washington University
* Dr. Donald Davendra
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <float.h>

#define POPULATION_SIZE 30
#define ELEMENTS 30
#define ITERATIONS 30
#define ILS_ITERATIONS 10
#define LS_STEP_SIZE 0.01
#define LS_TOLERANCE 1e-6


// Linear Congruential Generator (LCG) parameters
#define A 1664525
#define C 1013904223
#define M UINT_MAX

unsigned int seed = 42;

// XORShift parameters
unsigned int x = 123456789;
unsigned int y = 362436069;
unsigned int z = 521288629;
unsigned int w = 88675123;

// Linear Congruential Generator (LCG)
unsigned int lcg() {
    seed = (A * seed + C) % M;
    return seed;
}

// XORShift Generator
unsigned int xorshift() {
    unsigned int t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
}

// Benchmark functions
// Schwefel's function
double schwefel(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS; i++) {
        sum += -x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * ELEMENTS - sum;
}

// De Jong's (Sphere) function
double dejong(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

// Rosenbrock's function
double rosenbrock(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += 100 * pow((x[i + 1] - x[i] * x[i]), 2) + pow((1 - x[i]), 2);
    }
    return sum;
}

// Rastrigin's function
double rastrigin(double *x) {
    double sum = 10.0 * ELEMENTS;
    for (int i = 0; i < ELEMENTS; i++) {
        sum += x[i] * x[i] - 10.0 * cos(2 * M_PI * x[i]);
    }
    return sum;
}

// Griewangk's function
double griewangk(double *x) {
    double sum = 0.0;
    double product = 1.0;
    for (int i = 0; i < ELEMENTS; i++) {
        sum += x[i] * x[i] / 4000.0;
        product *= cos(x[i] / sqrt(i + 1));
    }
    return sum - product + 1;
}

// Sine Envelope Sine Wave function
double sine_envelope(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += 0.5 + (sin(sqrt(x[i] * x[i] + x[i + 1] * x[i + 1])) - 0.5) / pow(1 + 0.001 * (x[i] * x[i] + x[i + 1] * x[i + 1]), 2);
    }
    return -sum;
}

// Stretched V Sine Wave function
double stretched_v(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += pow((x[i] * x[i] + x[i + 1] * x[i + 1]), 0.25) * (sin(50 * pow((x[i] * x[i] + x[i + 1] * x[i + 1]), 0.1)) + 1);
    }
    return sum;
}

// Ackley's One function
double ackley_one(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += exp(-0.2) * sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]) + 3 * (cos(2 * x[i]) + sin(2 * x[i + 1]));
    }
    return sum;
}

// Ackley's Two function
double ackley_two(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += 20 + exp(1) - 20 * exp(-0.2 * sqrt((x[i] * x[i] + x[i + 1] * x[i + 1]) / 2)) - exp(0.5 * (cos(2 * M_PI * x[i]) + cos(2 * M_PI * x[i + 1])));
    }
    return sum;
}

// Egg Holder function
double egg_holder(double *x) {
    double sum = 0.0;
    for (int i = 0; i < ELEMENTS - 1; i++) {
        sum += -x[i] * sin(sqrt(fabs(x[i] - x[i + 1] - 47))) - (x[i + 1] + 47) * sin(sqrt(fabs(x[i + 1] + 47 + x[i] / 2)));
    }
    return sum;
}

// Initialize population
void initialize_population(double population[POPULATION_SIZE][ELEMENTS], double lower_bound, double upper_bound, unsigned int (*prng)()) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        for (int j = 0; j < ELEMENTS; j++) {
            double rand_val = (double)prng() / M;
            population[i][j] = lower_bound + rand_val * (upper_bound - lower_bound);
        }
    }
}

// Evaluate population
void evaluate_population(double population[POPULATION_SIZE][ELEMENTS], double (*func)(double *), double *results) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        results[i] = func(population[i]);
    }
}

// Comparison function for qsort
static int compare(const void *a, const void *b) {
    double diff = (*(double *)a - *(double *)b);
    return (diff > 0) - (diff < 0);
}

// Calculate statistical analysis
void statistical_analysis(double *results, int size, double *avg, double *std_dev, double *range, double *median) {
    double sum = 0.0;
    double min_val = DBL_MAX;
    double max_val = -DBL_MAX;
    for (int i = 0; i < size; i++) {
        sum += results[i];
        if (results[i] < min_val) {
            min_val = results[i];
        }
        if (results[i] > max_val) {
            max_val = results[i];
        }
    }
    *avg = sum / size;
    *range = max_val - min_val;

    // Calculate standard deviation
    double variance_sum = 0.0;
    for (int i = 0; i < size; i++) {
        variance_sum += pow(results[i] - *avg, 2);
    }
    *std_dev = sqrt(variance_sum / size);

    // Calculate median
    qsort(results, size, sizeof(double), &compare);
    if (size % 2 == 0) {
        *median = (results[size / 2 - 1] + results[size / 2]) / 2.0;
    } else {
        *median = results[size / 2];
    }
}

// Local Search (Gradient Descent)
void local_search(double *best_solution, double (*func)(double *)) {
    double gradient[ELEMENTS];
    double current_solution[ELEMENTS];
    double best_fitness = func(best_solution);
    double new_fitness;

    for (int i = 0; i < ELEMENTS; i++) {
        current_solution[i] = best_solution[i];
    }

    int iterations = 0;
    while (1) {
        // Calculate gradient (numerical approximation)
        for (int i = 0; i < ELEMENTS; i++) {
            double original_value = current_solution[i];
            current_solution[i] += LS_STEP_SIZE;
            double forward_value = func(current_solution);
            current_solution[i] = original_value - LS_STEP_SIZE;
            double backward_value = func(current_solution);
            gradient[i] = (forward_value - backward_value) / (2 * LS_STEP_SIZE);
            current_solution[i] = original_value;
        }

        // Update solution
        for (int i = 0; i < ELEMENTS; i++) {
            current_solution[i] -= LS_STEP_SIZE * gradient[i];
        }

        // Evaluate new solution
        new_fitness = func(current_solution);
        if (fabs(new_fitness - best_fitness) < LS_TOLERANCE) {
            break;
        }

        best_fitness = new_fitness;
        for (int i = 0; i < ELEMENTS; i++) {
            best_solution[i] = current_solution[i];
        }

        iterations++;
    }
}

// Iterated Local Search (ILS)
double iterated_local_search(unsigned int (*prng)(), double (*func)(double *), double lower_bound, double upper_bound) {
    double population[POPULATION_SIZE][ELEMENTS];
    double results[POPULATION_SIZE];
    double best_solution[ELEMENTS];
    double best_fitness = DBL_MAX;

    for (int ils_iter = 0; ils_iter < ILS_ITERATIONS; ils_iter++) {
        // Blind Search (BS)
        initialize_population(population, lower_bound, upper_bound, prng);
        evaluate_population(population, func, results);

        // Find the best individual in the population
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (results[i] < best_fitness) {
                best_fitness = results[i];
                for (int j = 0; j < ELEMENTS; j++) {
                    best_solution[j] = population[i][j];
                }
            }
        }

        // Local Search (LS)
        local_search(best_solution, func);
        double ls_fitness = func(best_solution);

        if (ls_fitness < best_fitness) {
            best_fitness = ls_fitness;
        }
    }

    return best_fitness;
}


int main() {
    double population[POPULATION_SIZE][ELEMENTS];
    double results[POPULATION_SIZE];
    double avg, std_dev, range, median;
    clock_t start, end;
    double time_spent;

    // Blind Search for all 10 functions over 30 iterations
    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Schwefel's function with LCG
        start = clock();
        initialize_population(population, -512, 512, lcg);
        evaluate_population(population, schwefel, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Schwefel's function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Schwefel's function with XORShift
        start = clock();
        initialize_population(population, -512, 512, xorshift);
        evaluate_population(population, schwefel, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Schwefel's function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // De Jong's function with LCG
        start = clock();
        initialize_population(population, -100, 100, lcg);
        evaluate_population(population, dejong, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - De Jong's function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // De Jong's function with XORShift
        start = clock();
        initialize_population(population, -100, 100, xorshift);
        evaluate_population(population, dejong, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - De Jong's function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Rastrigin's function with LCG
        start = clock();
        initialize_population(population, -30, 30, lcg);
        evaluate_population(population, rastrigin, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Rastrigin's function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Rastrigin's function with XORShift
        start = clock();
        initialize_population(population, -30, 30, xorshift);
        evaluate_population(population, rastrigin, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Rastrigin's function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Griewangk's function with LCG
        start = clock();
        initialize_population(population, -500, 500, lcg);
        evaluate_population(population, griewangk, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Griewangk's function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Griewangk's function with XORShift
        start = clock();
        initialize_population(population, -500, 500, xorshift);
        evaluate_population(population, griewangk, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Griewangk's function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Sine Envelope Sine Wave function with LCG
        start = clock();
        initialize_population(population, -30, 30, lcg);
        evaluate_population(population, sine_envelope, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Sine Envelope Sine Wave function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Sine Envelope Sine Wave function with XORShift
        start = clock();
        initialize_population(population, -30, 30, xorshift);
        evaluate_population(population, sine_envelope, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Sine Envelope Sine Wave function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Stretched V Sine Wave function with LCG
        start = clock();
        initialize_population(population, -30, 30, lcg);
        evaluate_population(population, stretched_v, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Stretched V Sine Wave function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Stretched V Sine Wave function with XORShift
        start = clock();
        initialize_population(population, -30, 30, xorshift);
        evaluate_population(population, stretched_v, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Stretched V Sine Wave function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Ackley's One function with LCG
        start = clock();
        initialize_population(population, -32, 32, lcg);
        evaluate_population(population, ackley_one, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Ackley's One function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Ackley's One function with XORShift
        start = clock();
        initialize_population(population, -32, 32, xorshift);
        evaluate_population(population, ackley_one, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Ackley's One function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Ackley's Two function with LCG
        start = clock();
        initialize_population(population, -32, 32, lcg);
        evaluate_population(population, ackley_two, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Ackley's Two function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Ackley's Two function with XORShift
        start = clock();
        initialize_population(population, -32, 32, xorshift);
        evaluate_population(population, ackley_two, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Ackley's Two function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Egg Holder function with LCG
        start = clock();
        initialize_population(population, -500, 500, lcg);
        evaluate_population(population, egg_holder, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Egg Holder function (LCG):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {
        // Egg Holder function with XORShift
        start = clock();
        initialize_population(population, -500, 500, xorshift);
        evaluate_population(population, egg_holder, results);
        end = clock();
        time_spent = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;

        statistical_analysis(results, POPULATION_SIZE, &avg, &std_dev, &range, &median);
        printf("Iteration %d - Egg Holder function (XORShift):\n", iter + 1);
        printf("  Average: %f\n  Standard Deviation: %f\n  Range: %f\n  Median: %f\n  Time (ms): %f\n\n", avg, std_dev, range, median, time_spent);
    }


    return 0;
}

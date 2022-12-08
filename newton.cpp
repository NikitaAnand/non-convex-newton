#include <iostream>
#include <cmath>
#include <random>
#include <getopt.h>
#include <time.h>

using namespace std;

#define PI 3.14175

void usage(const char* p) {
    printf("Usage: %s [options]\n", p);
    printf("Program Options:\n");
    printf("  -N  Number of Points\n");
    printf("  -f  Function to optimize\n");
}

int checkInputArguments(int n, int f)
{
    if (n == -1) {
        printf("Number of points was not specified\n");
        //return -1;
    }

    if (f == -1) {
        printf("Function to optimize was not specified\n");
        return -1;
    }

    return 0;
}


// Define the non-convex function
double function(double x, double y, int f) {
    switch (f)
    {
        case 0:
            return pow(x + 1, 2) + pow(y, 2);
        //Polynomial curve with 2 minima
        case 1:
            return pow(x - 2, 4) + pow(x - 2 * y, 2);
        //Rastrigin Function 2D
        case 2:
            return 20 + pow(x, 2) - (10 * cos(2 * PI * x)) + pow(y, 2) - (10 * cos(2 * PI * y));
        default:
            return pow(x, 4) - 3 * pow(x, 2) + y;
    }
    
}

// Define the gradient of the non-convex function
pair<double, double> function_gradient(double x, double y, int f) {
    switch(f)
    {
        case 0:
            return make_pair(2 * x + 2, 2 * y);
        //Derivative of polynomial curve with 2 minima
        case 1:
            return make_pair(2 * (2 * pow(x - 2, 3) + x - 2 * y), 8 * y - 4 * x);
        //Derivative of Rastrigin Function
        case 2:
            return make_pair((2 * x + 20 * PI * sin(2 * PI * x)), (2 * y + 20 * PI * sin(2 * PI * y)));
        default:
            return make_pair(4 * pow(x, 3) - 6 * x, 1);
    }
    
}

// Define the Hessian matrix of the non-convex function
pair<pair<double, double>, pair<double, double>> hessian_f(double x, double y, int f) {
    switch(f)
    {
        case 0:
            return make_pair(make_pair(2, 0),
                make_pair(0, 2));
        //Hessian of Polynomial curve with 2 minima
        case 1:
            return make_pair(make_pair(12*pow(x-2,2) + 2, -4),
                make_pair(-4, 8));
        //Hessian of Rastrigin function
        case 2:
            return make_pair(make_pair((40 * pow(PI, 2) * cos(2 * PI * x) + 2), 0),
                make_pair(0, (40 * pow(PI, 2) * cos(2 * PI * y) + 2)));
        default:
            return make_pair(make_pair(12 * pow(x, 2) - 6, 0),
                make_pair(0, 0));
    }
    
}

// Perform a single iteration of Newton's Method at a given point (x, y)
pair<double, double> newton_iteration(double x, double y, int f) {
    // Compute the gradient and Hessian at the current point
    auto gradient = function_gradient(x, y, f);
    auto hessian = hessian_f(x, y, f);

    // Compute the update using the gradient and Hessian
    double dx = (gradient.first * hessian.second.second - gradient.second * hessian.first.second)
        / (hessian.first.first * hessian.second.second - hessian.first.second * hessian.second.first);
    double dy = (gradient.second * hessian.first.first - gradient.first * hessian.second.first)
        / (hessian.first.first * hessian.second.second - hessian.first.second * hessian.second.first);

    // Return the updated point
    return make_pair(x - dx, y - dy);
}

// Perform Newton's Method with a multi-start approach to find the global minimum of the non-convex function
pair<double, double> cpu_newton_global_minimum(int N, int f) {
    // Set the number of random starting points and the number of iterations per starting point
    const int num_iterations = 100;

    // Set the lower and upper bounds for the random starting points
    const double lower_bound = -500;
    const double upper_bound = 500;

    // Create a random number generator to generate starting points
    mt19937 generator(random_device{}());
    uniform_real_distribution<double> distribution(lower_bound, upper_bound);

    // Set the initial minimum to the maximum possible value
    double min_f = numeric_limits<double>::max();
    pair<double, double> min_x = make_pair(0, 0);

    // Generate random starting points and perform Newton's Method from each starting point
    for (int i = 0; i < N; ++i) {
        // Generate a random starting point
        double x = distribution(generator);
            double y = distribution(generator);

            // Perform Newton's Method for the specified number of iterations
            for (int j = 0; j < num_iterations; ++j) {
                auto new_xy = newton_iteration(x, y, f);
                x = new_xy.first;
                y = new_xy.second;
            }

        // Update the minimum if the current point has a lower value of the function
        double f_x = function(x, y, f);
        if (f_x < min_f) {
            min_f = f_x;
            min_x = make_pair(x, y);
        }
    }

    // Return the global minimum
    return min_x;
}

int main(int argc, char* argv[]) {

    clock_t start, end;
    double time_taken;

    int opt;
    int f = -1;
    int N = -1;
    while ((opt = getopt(argc, argv, "N:f:")) != EOF) {
        switch (opt) {
        case 'N':
            N = atoi(optarg);
            break;
        case 'f':
            f = atoi(optarg);
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (checkInputArguments(N, f) == -1) {
        usage(argv[0]);
        return 1;
    }

    //Sequential Version

    printf("Newton's Method sequential version:\n");

    start = clock();
    auto cpu_global_minimum = cpu_newton_global_minimum(N, f);
    end = clock();

    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("Time taken = %lf\n", time_taken);
    cout << "The global minimum is at (" << cpu_global_minimum.first << ", " << cpu_global_minimum.second << ")" << endl;

    return 0;
}

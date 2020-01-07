#define EIGEN_DONT_PARALLELIZE

#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <string>
#include <fstream>
#include <Eigen/Dense>

void initial_row(double* const&, double* const* const&, Eigen::MatrixXd&, const int&, const int&);
void initial_col(double* const&, double* const* const&, Eigen::MatrixXd&, const int&, const int&);
void operate_row(double* const&, double* const* const&, Eigen::MatrixXd&, const int&, const int&);
void operate_col(double* const&, double* const* const&, Eigen::MatrixXd&, const int&, const int&);

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "No. of input: " << argc << std::endl;
        puts("./a.out <dim> <nrep>");
        exit(0);
    }
    
    int n = std::stoi(argv[1]);
    int nrep = std::stoi(argv[2]);

    double* x = new double[n * n];
    double** xx = new double*[n];
    for (size_t i = 0; i < n; i++) {
        xx[i] = new double[n];
    }
    Eigen::MatrixXd X(n, n);
    
    // Initialize Matrix
    initial_row(x, xx, X, n, nrep);
    initial_col(x, xx, X, n, nrep);
    
    // Matrix Elementwise Operation
    operate_row(x, xx, X, n, nrep);
    operate_col(x, xx, X, n, nrep);

    for (size_t i = 0; i < n; i++) {
        delete[] xx[i];
    }
    delete[] xx;
    delete[] x;
 
    return 0;
}

void initial_row(double* const& x, double* const* const& xx, Eigen::MatrixXd& X, const int& n, const int& nrep) {
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed;
    double t = 0;
    
    // Vectorized Matrix
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                x[i * n + j] = i < j ? i : j;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "x* row time: " << t / nrep << std::endl;

    // Double array Matrix
    t = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                xx[i][j] = i < j ? i : j;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "xx row time: " << t / nrep << std::endl;

    // Eigen Matrix
    t = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                X(i, j) = i < j ? i : j;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "e* row time: " << t / nrep << std::endl;
}

void initial_col(double* const& x, double* const* const& xx, Eigen::MatrixXd& X, const int& n, const int& nrep) {
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed;
    double t = 0;
    
    // Vectorized Matrix
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                x[j * n + i] = j < i ? j : i;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "x* col time: " << t / nrep << std::endl;

    // Double array Matrix
    t = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                xx[j][i] = j < i ? j : i;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "xx col time: " << t / nrep << std::endl;

    // Eigen Matrix
    t = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                X(j, i) = j < i ? j : i;
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "e* col time: " << t / nrep << std::endl;
}

void operate_row(double* const& x, double* const* const& xx, Eigen::MatrixXd& X, const int& n, const int& nrep) {
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed;
    double t = 0;
    double s = 0.0;
    
    // Vectorized Matrix
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                s += std::sqrt(x[i * n + j]);
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " x* row time: " << t / nrep << std::endl;

    // Double array Matrix
    t = 0.0;
    s = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                s += std::sqrt(xx[i][j]);
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " xx row time: " << t / nrep << std::endl;

    // Eigen Matrix
    t = 0.0;
    s = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
               s += std::sqrt(X(i, j));
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " e* row time: " << t / nrep << std::endl;

    // Eigen Matrix Function
    t = 0.0;
    s = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        s += X.array().sqrt().sum();
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " ee all time: " << t / nrep << std::endl;
}

void operate_col(double* const& x, double* const* const& xx, Eigen::MatrixXd& X, const int& n, const int& nrep) {
    auto t0 = std::chrono::system_clock::now();
    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed;
    double t = 0;
    double s = 0.0;
    
    // Vectorized Matrix
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                s += std::sqrt(x[j * n + i]);
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " x* col time: " << t / nrep << std::endl;

    // Double array Matrix
    t = 0.0;
    s = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                s += std::sqrt(xx[j][i]);
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " xx col time: " << t / nrep << std::endl;

    // Eigen Matrix
    t = 0.0;
    s = 0.0;
    for (size_t z = 0; z < nrep; z++) {
        t0 = std::chrono::system_clock::now();
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                s += std::sqrt(X(j, i));
            }
        }
        t1 = std::chrono::system_clock::now();
        elapsed = t1 - t0;
        t += static_cast<double>(elapsed.count());
    }
    std::cout << "s: " << s / nrep << " e* col time: " << t / nrep << std::endl;
}
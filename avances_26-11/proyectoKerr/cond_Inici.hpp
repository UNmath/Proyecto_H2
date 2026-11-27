
#pragma once

#include "/usr/include/eigen3/Eigen/Dense"
#include <array>

namespace coin {

using Vector5d = Eigen::Matrix<double, 5, 1>;

double A(double r, double theta, double delta, double a);
double B(double A1, double theta, double sigma);
double alpha(double A2, double sigma, double delta);

std::array<double, 6> condini(const Vector5d& estado, double x, double y,
                              double sigma, double delta, double a);

} 

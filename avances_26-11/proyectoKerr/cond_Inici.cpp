#include "cond_Inici.hpp"
#include <cmath>

using Vector5d = Eigen::Matrix<double, 5, 1>;
namespace coin {

double A(double r, double theta, double delta, double a){
    return (r*r + a*a)*(r*r + a*a) - a*a*delta*std::sin(theta)*std::sin(theta);
}
double B(double A1, double theta, double sigma){
    return std::sqrt(A1/sigma)*std::sin(theta);
}
double alpha(double A2, double sigma, double delta){
    return std::sqrt(sigma*delta/A2);
}

std::array<double, 6> condini(const Vector5d& estado, double x, double y, double sigma, double delta, double a){
    double A3 = A(estado[0], estado[1], delta, a);
    double B3 = B(A3, estado[1], sigma);
    double alpha3 = alpha(A3, sigma, delta);
    double w3 = (2.0 * a * estado[0]) / A3;

    double nr = -std::cos(y);
    double nth = std::sin(y) * std::cos(x);
    double nph = std::sin(y) * std::sin(x);

    double t_dot = 1.0 / alpha3;
    double r_dot = std::sqrt(delta / sigma) * nr;
    double th_dot = nth / std::sqrt(sigma);
    double ph_dot = (w3 / alpha3) + (nph / B3);

    double pr = std::sqrt(sigma / delta) * nr;
    double pth = std::sqrt(sigma) * nth;

    return std::array<double, 6>{t_dot, r_dot, th_dot, ph_dot, pr, pth};

}

}

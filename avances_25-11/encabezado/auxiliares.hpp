#ifndef AUXILIARES_HPP
#define AUXILIARES_HPP

#include "/usr/include/eigen3/Eigen/Dense"

namespace aux {

using Vector8d = Eigen::Matrix<double, 8, 1>;

double sigma(const Vector8d& estado, double a);
double delta(const Vector8d& estado, double a);

double calcular_Lz(const Vector8d& estado, double a, double E);
double calcular_phi_dot(const Vector8d& estado, double a, double E, double Lz);
double calcular_Q(const Vector8d& estado, double a, double E);
double calcular_kappa(const Vector8d& estado, double Lz, double a, double E);
double calcularTDot(const Vector8d& estado, double a, double E, double Lz);
double calcularE2(const Vector8d& estado, double a);
double calcular_pr_dot(const Vector8d& estado, double a, double E, double kappa);
double calcular_ptheta_pun(const Vector8d& estado, double Lzi, double Ei);

Vector8d rk4(double a, double h, double Ei, double Lzi, const Vector8d& estado);

} 

#endif 
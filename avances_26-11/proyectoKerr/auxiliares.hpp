#pragma once
#include "/usr/include/eigen3/Eigen/Dense"

// Convención del vector de estado:
// [0]: r, [1]: theta, [2]: phi, [3]: pr, [4]: ptheta

using Vector5d = Eigen::Matrix<double, 5, 1>;

namespace aux {

double sigma(const Vector5d& estado, double a);
double delta(const Vector5d& estado, double a);

double calcular_r_dot(const Vector5d& estado, double sigma, double delta);
double calcular_theta_dot(const Vector5d& estado, double sigma);
double calcular_phi_dot(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta);
double calcularTDot(const Vector5d& estado, double a, double E, double Lz, double sigma, double delta);

double calcular_Lz(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta);
double calcularE2(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta);

double calcular_Q(const Vector5d& estado, double a, double E, double Lz);
double calcular_kappa(const Vector5d& estado, double Lz, double a, double E);

double calcular_ptheta_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma);
double calcular_pr_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma, double delta, double kappa);

double calcularEi2(const Vector5d& estado, std::array<double, 6> condi, double a, double sigma, double delta);

double calcular_Lzi(const Vector5d& estado, std::array<double, 6> condi, double a, double Ei, double sigma, double delta);

} // namespace aux

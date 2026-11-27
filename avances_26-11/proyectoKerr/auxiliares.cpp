#include "auxiliares.hpp"
#include <cmath>

// el vector de estado tendrá 
// [0]: r, [1]: theta, [2]: phi, [3]: pr, [4]: ptheta
namespace aux {

double sigma(const Vector5d& estado, double a) {
    double c = std::cos(estado(1));
    return estado(0) * estado(0) + a * a * c * c;
}

double delta(const Vector5d& estado, double a) {
    return estado(0) * estado(0) + a * a - 2.0 * estado(0);
}

double calcular_r_dot(const Vector5d& estado, double sigma, double delta) {
    return estado(3) * delta / sigma;
}

double calcular_theta_dot(const Vector5d& estado, double sigma) {
    return estado(4) / sigma;
}

double calcular_Lz(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta) {
    double r = estado(0), theta = estado(1);
    double phi_dot = calcular_phi_dot(estado, a, Ei, Lzi, sigma, delta); // ajusta si Lz es dato previo
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;

    double numerador = (sigma * delta * phi_dot - 2.0 * a * r * Ei) * sin2_theta;
    double denominador = sigma - 2.0 * r;
    return numerador / denominador;
}

double calcular_Lzi(const Vector5d& estado, std::array<double, 6> condi, double a, double Ei, double sigma, double delta) {
    double r = estado(0), theta = estado(1);
    double phi_dot = condi[3]; // ajusta si Lz es dato previo
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;

    double numerador = (sigma * delta * phi_dot - 2.0 * a * r * Ei) * sin2_theta;
    double denominador = sigma - 2.0 * r;
    return numerador / denominador;
}



double calcularE2(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta) {
    double r = estado(0), theta = estado(1);
    double r_dot = calcular_r_dot(estado, sigma, delta);
    double theta_dot = calcular_theta_dot(estado, sigma);
    double phi_dot = calcular_phi_dot(estado, a, Ei, Lzi, sigma, delta);

    double sin_theta = std::sin(theta);

    return ((sigma - 2.0 * r) / delta) * ((r_dot * r_dot) + delta * (theta_dot * theta_dot))
         + delta * (phi_dot * phi_dot) * (sin_theta * sin_theta);
}

double calcularEi2(const Vector5d& estado, std::array<double, 6> condi, double a, double sigma, double delta) {
    double r = estado(0), theta = estado(1);
    double r_dot = condi[1];
    double theta_dot = condi[2];
    double phi_dot = condi[3];

    double sin_theta = std::sin(theta);

    return ((sigma - 2.0 * r) / delta) * ((r_dot * r_dot) + delta * (theta_dot * theta_dot))
         + delta * (phi_dot * phi_dot) * (sin_theta * sin_theta);
}


double calcular_phi_dot(const Vector5d& estado, double a, double Ei, double Lzi, double sigma, double delta) {
    double r = estado(0), theta = estado(1);
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;

    double numerador = 2.0 * a * r * Ei * sin2_theta + (sigma - 2.0 * r) * Lzi;
    double denominador = sigma * delta * sin2_theta;
    return numerador / denominador;
}

double calcularTDot(const Vector5d& estado, double a, double E, double Lz, double sigma, double delta) {
    double r = estado(0);
    double numerador = 2.0 * r * ((r * r) + (a * a)) * E - 2.0 * a * r * Lz;
    return E + numerador / (sigma * delta);
}

double calcular_Q(const Vector5d& estado, double a, double E, double Lz) {
    double theta = estado(1), p_theta = estado(4);

    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin2_theta = sin_theta * sin_theta;
    double cos2_theta = cos_theta * cos_theta;

    return (p_theta * p_theta) + ((Lz * Lz) / sin2_theta - (a * a) * (E * E)) * cos2_theta;
}

double calcular_kappa(const Vector5d& estado, double Lz, double a, double E) {
    double Q = calcular_Q(estado, a, E, Lz);
    return Q + (Lz * Lz) + (a * a) * (E * E);
}

double calcular_ptheta_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma) {
    double theta = estado(1);
    double s = std::sin(theta);
    double s2 = s * s;
    return (std::sin(2.0 * theta) / (2.0 * sigma)) * ((Lzi * Lzi) / (s2 * s2) - (a * a) * (Ei * Ei));
}

double calcular_pr_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma, double delta, double kappa) {
    double r = estado(0), pr = estado(3);
    double num1 = 2.0 * r * (r * r + a * a) * (Ei * Ei) - kappa * (r - 1.0) - 2.0 * a * Ei * Lzi;
    return (num1 / (sigma * delta)) - ((2.0 * pr * pr * (r - 1.0)) / sigma);
}

} // namespace aux


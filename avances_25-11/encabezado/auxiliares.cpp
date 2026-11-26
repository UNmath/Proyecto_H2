#include "auxiliares.hpp"
#include <cmath>


namespace aux {

double sigma(const Vector8d& estado, double a) {
    double c = std::cos(estado(2));
    return estado(0) * estado(0) + a * a * c * c;
}

double delta(const Vector8d& estado, double a) {
    return estado(0) * estado(0) + a * a - 2.0 * estado(0);
}

double calcular_Lz(const Vector8d& estado, double a, double E) {
    double r = estado(0), theta = estado(2), phi_dot = estado(5);   
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;
    double sig = sigma(estado, a), del = delta(estado, a); 

    double numerador = (sig * del * phi_dot - 2.0 * a * r * E) * sin2_theta;
    double denominador = sig - 2.0 * r;
    return numerador / denominador;
}

double calcularE2(const Vector8d& estado, double a) {
    double r = estado(0), r_dot = estado(1);
    double theta = estado(2), theta_dot = estado(3);
    double phi_dot = estado(5);
    double sig = sigma(estado, a), del = delta(estado, a); 

    double sin_theta = std::sin(theta);

    return ((sig - 2.0 * r) / del) * ((r_dot * r_dot) + del * (theta_dot * theta_dot))
         + del * (phi_dot * phi_dot) * (sin_theta * sin_theta); 
}

double calcular_phi_dot(const Vector8d& estado, double a, double E, double Lz) {
    double r = estado(0), theta = estado(2);
    double sin_theta = std::sin(theta);
    double sin2_theta = sin_theta * sin_theta;
    double sig = sigma(estado, a), del = delta(estado, a);

    double numerador = 2.0 * a * r * E + (sig - 2.0 * r) * Lz;
    double denominador = sig * del * sin2_theta;
    return numerador / denominador;
}

double calcularTDot(const Vector8d& estado, double a, double E, double Lz) {
    double r = estado(0);
    double numerador = 2.0 * r * ((r * r) + (a * a)) * E - 2.0 * a * r * Lz; 
    return E + numerador / (sigma(estado, a) * delta(estado, a)); 
}

double calcular_Q(const Vector8d& estado, double a, double E) {
    double theta = estado(2), p_theta = estado(7);
    double Lz = calcular_Lz(estado, a, E);

    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin2_theta = sin_theta * sin_theta;
    double cos2_theta = cos_theta * cos_theta;

    return (p_theta * p_theta) + ((Lz * Lz) / sin2_theta - (a * a) * (E * E)) * cos2_theta;
}

double calcular_kappa(const Vector8d& estado, double Lz, double a, double E) {
    double Q = calcular_Q(estado, a, E);
    return Q + (Lz * Lz) + (a * a) * (E * E);
}

double calcular_pr_dot(const Vector8d& estado, double a, double E, double kappa) {
    double r = estado(0);
    double pr = estado(6);
    
    double sig = sigma(estado, a);
    double del = delta(estado, a);
    double Lz = calcular_Lz(estado, a, E);

    double termino1 = (2.0 * r * (r * r + a * a) * (E * E) - kappa * (r - 1.0) - 2.0 * a * E * Lz) / (sig * del);
    double termino2 = (2.0 * (pr * pr) * (r - 1.0)) / sig;

    return termino1 - termino2;
}

double calcular_ptheta_pun(const Vector8d& estado, double Lzi, double Ei) {
    return 0.0;
}

Vector8d rk4(double a, double h, double Ei, double Lzi, Vector8d& estado) {
    double k1theta = 0.0, k2theta = 0.0, k3theta = 0.0, k4theta = 0.0;
    double k1r=0.0,k2r=0.0, k3r=0.0, k4r=0.0;
    double k1pr=0.0,k2pr=0.0,k3pr=0.0,k4pr=0.0;

    //k1theta = h * calcular_phi_dot(estado);
    //k1r = h * calcularprpunto(estado);
    //k1pr = h; 
    return estado;
}
}
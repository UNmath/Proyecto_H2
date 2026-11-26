#include <iostream>
#include <fstream>
#include "/usr/include/eigen3/Eigen/Dense"

// Alias para mayor claridad y facilidad al escribir el tipo de vector de 5 elementos
using Vector5d = Eigen::Matrix<double, 5, 1>;

// Orden del vector estado:
// [0]: r, [1]: theta, [2]: phi, [3]: pr, [4]: ptheta

double sigma(Vector5d& estado, double a) {return estado(0) * estado(0) + a * a * std::cos(estado(2)) * std::cos(estado(2));}
double delta(Vector5d& estado, double a) {return estado(0) * estado(0) + a * a - 2.0 * estado(0);}

// Prototipos (ahora reciben sig y del como parámetros)
double calcular_Lz(Vector5d& estado, double a, double E, double sig, double del);
double calcular_phi_dot(Vector5d& estado, double a, double E, double Lz, double sig, double del);
double calcular_Q(Vector5d& estado, double a, double E, double sig, double del);
double calcular_kappa(Vector5d& estado, double Lz, double a, double E, double sig, double del);
double calcularTDot(Vector5d& estado, double a, double E, double Lz, double sig, double del);
double calcularE2(Vector5d& estado, double a, double sig, double del);
double calcular_pr_dot(Vector5d& estado, double a, double E, doubl, double sig, double del);
double calcular_ptheta_dot(Vector5d& estado, double Lzi, double Ei);
double calcular_pr(Vector5d& estado, double del, double sig);
double calcular_ptheta(Vector5d& estado, double sig); 
Vector5d rkf45(double a, double h, double Ei, double Lzi, Vector5d& estado);

Vector5d estado = Vector5d::Zero();

int main() {
    double a = 0.9;
    
    // Ejemplo de cálculo único de sig y del para un estado dado:
    double r = state(0);
    double cos_theta = std::cos(state(1));
    double sig = sigma(estado(0), a,);
    double del = delta(estado(0), a);

    return 0; 
}

// funciones auxiliares
double sigma(const Vector5d& estado, double a);
double delta(const Vector5d& estado, double a);

double calcular_r_dot(const Vector5d& estado, double sigma, double delta);
double calcular_theta_dot(const Vector5d& estado, double sigma);
double calcular_phi_dot(const Vector5d& estado, double a, double E, double Lz, double sigma, double delta);
double calcularTDot(const Vector5d& estado, double a, double E, double Lz, double sigma, double delta);

double calcular_Lz(const Vector5d& estado, double a, double E, double sigma, double delta);
double calcularE2(const Vector5d& estado, double a, double sigma, double delta);
double calcular_Q(const Vector5d& estado, double a, double E, double Lz);
double calcular_kappa(const Vector5d& estado, double Lz, double a, double E);

double calcular_ptheta_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma);
double calcular_pr_dot(const Vector5d& estado, double Lzi, double Ei, double a, double sigma, double delta);

// pa que se guíe

Vector5d rkf45(Vector5d& estado, double del, double sig, double a, double h, double Ei, double Lzi, double kappa) {
    Vector5d k1 = Vector5d::Zero();
    Vector5d k2 = Vector5d::Zero();
    Vector5d k3 = Vector5d::Zero();
    Vector5d k4 = Vector5d::Zero();
    Vector5d k5 = Vector5d::Zero();
    Vector5d k6 = Vector5d::Zero();
    
    k1 << h * calcular_pr(estado, del, sig), h * calcular_ptheta(estado, sig),h * calcular_phi_dot(estado, a, Ei, Lzi, sig, del), h * calcular_pr_dot(estado, a, Ei, sig, del, kappa), h * calcular_ptheta_dot(estado, Lzi, Ei); 
    Vector5d estado2 = estado + 0.25 * k1;
    double del2= delta(estado2, a), sig2 = sigma(estado2, a);
    k2 << h * calcular_pr(estado2, del2, sig2), h * calcular_ptheta(estado2, sig2),h * calcular_phi_dot(estado2, a, Ei, Lzi, sig2, del2), h * calcular_pr_dot(estado2, a, Ei, sig2, del2, kappa), h * calcular_ptheta_dot(estado2, Lzi, Ei); 
    Vector5d estado3 = estado + (3.0/32.0) * k1 + (9.0/32.0) * k2;
    double del3 = delta(estado3, a), sig3 = sigma(estado3, a);
    k3 << h * calcular_pr(estado3, del3, sig3), h * calcular_ptheta(estado3, sig3),h * calcular_phi_dot(estado3, a, Ei, Lzi, sig3, del3), h * calcular_pr_dot(estado3, a, Ei, sig3, del3, kappa), h * calcular_ptheta_dot(estado3, Lzi, Ei); 
    Vector5d estado4 = estado + (1932.0/2197.0) * k1 - (7200.0/2197.0) * k2 + (7296.0/2197.0) * k3;
    double del4 = delta(estado4, a), sig4 = sigma(estado4, a);
    k4 << h * calcular_pr(estado4, del4, sig4), h * calcular_ptheta(estado4, sig4),h * calcular_phi_dot(estado4, a, Ei, Lzi, sig4, del4), h * calcular_pr_dot(estado4, a, Ei, sig4, del4, kappa), h * calcular_ptheta_dot(estado4, Lzi, Ei); 
    Vector5d estado5 = estado + (439.0/216.0) * k1 - 8.0 * k2 + (3680.0/513.0) * k3 - (845.0/4104.0) * k4;
    double del5 = delta(estado5, a), sig5 = sigma(estado5, a);
    k5 << h * calcular_pr(estado5, del5, sig5), h * calcular_ptheta(estado5, sig5),h * calcular_phi_dot(estado5, a, Ei, Lzi, sig5, del5), h * calcular_pr_dot(estado5, a, Ei, sig5, del5, kappa), h * calcular_ptheta_dot(estado5, Lzi, Ei); 
    Vector5d estado6 = estado - (8.0/27.0) * k1 + 2 * k2 - (3544.0/2565.0) * k3 + (1859.0/4104.0) * k4 - (11.0/40.0) * k5;
    double del6 = delta(estado6, a), sig6 = sigma(estado6, a);
    k6 << h * calcular_pr(estado6, del6, sig6), h * calcular_ptheta(estado6, sig6),h * calcular_phi_dot(estado6, a, Ei, Lzi, sig6, del6), h * calcular_pr_dot(estado6, a, Ei, sig6, del6, kappa), h * calcular_ptheta_dot(estado6, Lzi, Ei); 


    Vector5d y4 = estado + (25.0/216.0) * k1 + (1408.0/2565.0) * k3 + (2197.0/4104.0) * k4 - (1.0/5.0) * k5;
    Vector5d y5 = estado + (16.0/135.0) * k1 + (6656.0/12825.0) * k3 + (28561.0/56430.0) * k4 - (9.0/50.0) * k5 + (2.0/55.0) * k6;

    return y5;
}


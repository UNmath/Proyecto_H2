#include "auxiliares.hpp"
#include "cond_Inici.hpp"
#include <iostream>
#include <fstream>
#include "/usr/include/eigen3/Eigen/Dense"

// [0]: r, [1]: theta, [2]: phi, [3]: pr, [4]: ptheta
Vector5d rkf45(Vector5d& estado, double del, double sig, double a, double& h, double Ei, double Lzi, double kappa, double tol, double h_max, double h_min);

struct Particula{
    int id;
    Vector5d estado;
    double Ei;            
    double Lzi;           
    double kappa;         
    double h;             
    bool activa;
}

Particula crear_particula(int id, Vector5d estado, double Ei, double Lzi, double kappa, double h) {
    Particula p;
    p.id = id;
    p.estado = estado;
    p.Ei = Ei;
    p.Lzi = Lzi;
    p.kappa = kappa;
    p.h = h;
    p.activa = true; // Inicialmente activa
    return p;
}

int main(){
    
    std::cout<<"Muchas gracias"<<"\n";
    
    double N=2* std::pow(10,6)

    for (int i = 0; i < N; ++i) {
        
        Vector5d state(12.0, M_PI / 3.0, 0.0, 0.0, 0.0);
    double a = 0.9;

    double delta = aux::delta(state, a);
    double sigma = aux::sigma(state, a);
    double xp = 0.6;
    double yp = 0.25;

    std::array<double, 6> condi = coin::condini(state, xp, yp, sigma, delta, a);
    state(3) = condi[4]; // pr inicial
    state(4) = condi[5]; // ptheta inicial

    double Ei = std::sqrt(aux::calcularEi2(state, condi, a, sigma, delta));
    double Lzi = aux::calcular_Lzi(state, condi, a, Ei, sigma, delta);
    double kappa = aux::calcular_kappa(state, Lzi, a, Ei);

    double h = 0.05; // Paso inicial
    double h_max = 1.0;
    double h_min = 1e-6;
    double r_horizonte = 1.0 + std::sqrt(1.0 - a * a); // Horizonte de eventos de Kerr
    double max_pasos = 10000;
    double tol = 1e-8;
        
        
    }

    std::ofstream archivo1("datosKerr.txt");
        for (int paso = 0; paso < max_pasos; ++paso) {
            delta = aux::delta(state, a);
            sigma = aux::sigma(state, a);
            state = rkf45(state, delta, sigma, a, h, Ei, Lzi, kappa, tol, h_max, h_min);
            
            double r     = state(0);
            double theta = state(1);
            double phi   = state(2);

            // Transformación de Boyer-Lindquist a Cartesianas (x, y, z)
            double radio_kerr = std::sqrt(r * r + a * a);
            double x = radio_kerr * std::sin(theta) * std::cos(phi);
            double y = radio_kerr * std::sin(theta) * std::sin(phi);
            double z = r * std::cos(theta);

            // Guardar en el archivo separados por comas
            archivo1 << paso << "," << r << "," << theta << "," << phi << ","
                    << x << "," << y << "," << z << "\n";

            if (state(0) <= r_horizonte * 1.01 || state(0) > 100.0) {
                break;
            }
        }
    archivo1.close();
    return 0; 
}

Vector5d rkf45(Vector5d& estado, double del, double sig, double a, double& h, double Ei, double Lzi, double kappa, double tol, double h_max, double h_min) {
    Vector5d k1 = Vector5d::Zero();
    Vector5d k2 = Vector5d::Zero();
    Vector5d k3 = Vector5d::Zero();
    Vector5d k4 = Vector5d::Zero();
    Vector5d k5 = Vector5d::Zero();
    Vector5d k6 = Vector5d::Zero();
    
    while(true) {
        k1 << h * aux::calcular_r_dot(estado, sig, del), h * aux::calcular_theta_dot(estado, sig), h * aux::calcular_phi_dot(estado, a, Ei, Lzi, sig, del), h * aux::calcular_pr_dot(estado, Lzi, Ei, a, sig, del, kappa), h * aux::calcular_ptheta_dot(estado, Lzi, Ei, a, sig); 
        Vector5d estado2 = estado + 0.25 * k1;
        double del2= aux::delta(estado2, a), sig2 = aux::sigma(estado2, a);
        k2 << h * aux::calcular_r_dot(estado2, sig2, del2), h * aux::calcular_theta_dot(estado2, sig2), h * aux::calcular_phi_dot(estado2, a, Ei, Lzi, sig2, del2), h * aux::calcular_pr_dot(estado2, Lzi, Ei, a, sig2, del2, kappa), h * aux::calcular_ptheta_dot(estado2, Lzi, Ei, a, sig2); 
        Vector5d estado3 = estado + (3.0/32.0) * k1 + (9.0/32.0) * k2;
        double del3 = aux::delta(estado3, a), sig3 = aux::sigma(estado3, a);
        k3 << h * aux::calcular_r_dot(estado3, sig3, del3), h * aux::calcular_theta_dot(estado3, sig3), h * aux::calcular_phi_dot(estado3, a, Ei, Lzi, sig3, del3), h * aux::calcular_pr_dot(estado3, Lzi, Ei, a, sig3, del3, kappa), h * aux::calcular_ptheta_dot(estado3, Lzi, Ei, a, sig3); 
        Vector5d estado4 = estado + (1932.0/2197.0) * k1 - (7200.0/2197.0) * k2 + (7296.0/2197.0) * k3;
        double del4 = aux::delta(estado4, a), sig4 = aux::sigma(estado4, a);
        k4 << h * aux::calcular_r_dot(estado4, sig4, del4), h * aux::calcular_theta_dot(estado4, sig4), h * aux::calcular_phi_dot(estado4, a, Ei, Lzi, sig4, del4), h * aux::calcular_pr_dot(estado4, Lzi, Ei, a, sig4, del4, kappa), h * aux::calcular_ptheta_dot(estado4, Lzi, Ei, a, sig4); 
        Vector5d estado5 = estado + (439.0/216.0) * k1 - 8.0 * k2 + (3680.0/513.0) * k3 - (845.0/4104.0) * k4;
        double del5 = aux::delta(estado5, a), sig5 = aux::sigma(estado5, a);
        k5 << h * aux::calcular_r_dot(estado5, sig5, del5), h * aux::calcular_theta_dot(estado5, sig5), h * aux::calcular_phi_dot(estado5, a, Ei, Lzi, sig5, del5), h * aux::calcular_pr_dot(estado5, Lzi, Ei, a, sig5, del5, kappa), h * aux::calcular_ptheta_dot(estado5, Lzi, Ei, a, sig5); 
        Vector5d estado6 = estado - (8.0/27.0) * k1 + 2 * k2 - (3544.0/2565.0) * k3 + (1859.0/4104.0) * k4 - (11.0/40.0) * k5;
        double del6 = aux::delta(estado6, a), sig6 = aux::sigma(estado6, a);
        k6 << h * aux::calcular_r_dot(estado6, sig6, del6), h * aux::calcular_theta_dot(estado6, sig6), h * aux::calcular_phi_dot(estado6, a, Ei, Lzi, sig6, del6), h * aux::calcular_pr_dot(estado6, Lzi, Ei, a, sig6, del6, kappa), h * aux::calcular_ptheta_dot(estado6, Lzi, Ei, a, sig6); 

        Vector5d y4 = estado + (25.0/216.0) * k1 + (1408.0/2565.0) * k3 + (2197.0/4104.0) * k4 - (1.0/5.0) * k5;
        Vector5d y5 = estado + (16.0/135.0) * k1 + (6656.0/12825.0) * k3 + (28561.0/56430.0) * k4 - (9.0/50.0) * k5 + (2.0/55.0) * k6;

        double error = (y5 - y4).lpNorm<Eigen::Infinity>();

        double s;
        if (error == 0.0) {
            s = 2.0; // Evita división por cero si el error es nulo
        } else {
            s = 0.84 * std::pow(tol / error, 0.25);
            s = std::max(0.1, std::min(s, 4.0));
        }

        // 4. Aceptar o rechazar el paso
        if (error <= tol || h <= h_min) {
            h = std::max(h_min, std::min(h * s, h_max)); // Nuevo h para el siguiente paso
            return y5;
        } else {
            h = std::max(h_min, h * s); // Reducir h y repetir el bucle
        }
    }

}
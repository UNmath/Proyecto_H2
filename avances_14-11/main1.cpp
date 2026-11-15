#include <cmath>
#include <iostream>
#include <array>
#include <fstream>

// se debe ajustar el N o el theta final para que no de tantas vueltas y se balla a la verguchi
const int N=4000;


double coeficiente(double M, double G, double c, double r);

void movimiento(double M, double h, double t0, double tf, std::array<double, N>& u_prim, std::array<double, N>& u_seg);

int main(){
    double G = 1, c=1, M=1, r0=10, L=12, E=2;
    double b = L/E;
    double u0=1/r0;
    double v0 = std::sqrt((1/(b*b) - u0*u0 * coeficiente(M, G, c, 1/u0)));
    double a=0;

    std::array<double, N> u_p;
    std::array<double, N> u_s;
    u_p[0] = u0;
    u_s[0] = v0;
    double t0 = 0;
    double tf = 5;
    double h = 0.001;

    movimiento(M, h, t0, tf, u_p, u_s);

    std::ofstream archivo1("datos.txt");
    for (int i = 0; i < N; i++) {
        double phi = i * h;
        double r = 1.0 / u_p[i];
        // Guardamos phi, u, y también (x, y) para graficar fácilmente
        double x = r * std::cos(phi);
        double y = r * std::sin(phi);
        archivo1 << phi << "," << u_p[i] << "," << x << "," << y << "\n";
    }
    archivo1.close();

    return 0;
}

double coeficiente(double M, double G, double c, double r){
    return (1 - (2*M*G/(c*c*r)));
}

double velInicial()

void movimiento(double M, double h, double t0, double tf, std::array<double, N>& u_prim, std::array<double, N>& u_seg){
    double k1x=0.0;
    double k2x=0.0;
    double k3x=0.0;
    double k4x=0.0;

    double k1v=0.0;
    double k2v=0.0;
    double k3v=0.0;
    double k4v=0.0;

    for (int i = 0; i<(N-1); i++ ){
        k1x = h * u_seg[i];
        k1v = h * (-u_prim[i] + 3*M*u_prim[i]*u_prim[i]);

        k2x = h * (u_seg[i] + k1v/2); 
        k2v = h * (-(u_prim[i] + k1x/2) + 3*M*(u_prim[i] + k1x/2)*(u_prim[i] + k1x/2));

        k3x = h * (u_seg[i] + k2v/2);
        k3v = h * (-(u_prim[i] + k2x/2) + 3*M*(u_prim[i] + k2x/2)*(u_prim[i] + k2x/2));

        k4x = h * (u_seg[i] + k3v);
        k4v = h * (-(u_prim[i] + k3x) + 3*M*(u_prim[i] + k3x)*(u_prim[i] + k3x));

        u_prim[i+1] = u_prim[i] + (k1x + 2*k2x + 2*k3x + k4x) / 6.0; 
        u_seg[i+1] = u_seg[i] + (k1v + 2*k2v + 2*k3v + k4v) / 6.0; 
    }
}

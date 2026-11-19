#include <cmath>
#include <fstream>
#include <vector>
#include <array>

const double a = -1.0;
const double b = 7.0;

using State = std::array<double, 2>; // {r, v = dr/dphi}

double Delta(double r) { return r*r + a*a - 2.0*r; }
double Dn(double r)    { return b*(r - 2.0) + 2.0*a; }
double Rr(double r)    { return r*r + 2.0/r*(a-b)*(a-b) + a*a - b*b; }

// G(r) = (dr/dphi)^2
double G(double r) {
    double d = Dn(r), De = Delta(r);
    return De*De*Rr(r)/(d*d);
}

// dG/dr analítica
double dG(double r) {
    double De = Delta(r), dDe = 2.0*r - 2.0;
    double d = Dn(r), R = Rr(r);
    double dR = 2.0*r - 2.0*(a-b)*(a-b)/(r*r);
    return (2.0*De*dDe*R + De*De*dR)/(d*d) - 2.0*De*De*R*b/(d*d*d);
}

State f(const State& y) {
    return { y[1], 0.5*dG(y[0]) };
}

State rk4_step(const State& y, double h) {
    auto add = [](const State& p, const State& q, double s) {
        return State{ p[0] + s*q[0], p[1] + s*q[1] };
    };
    State k1 = f(y);
    State k2 = f(add(y, k1, h/2));
    State k3 = f(add(y, k2, h/2));
    State k4 = f(add(y, k3, h));
    return { y[0] + h*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0,
             y[1] + h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0 };
}

int main() {
    const double PI = 3.14159265358979323846;
    const double r_start = 12.5;
    const double r_plus  = 1.0 + std::sqrt(1.0 - a*a); // horizonte
    const double h = 1e-3;                             // paso en phi
    const int Nmax = 20000;

    State y = { r_start, -std::sqrt(G(r_start)) };     // entrante
    double phi = 0.0;

    std::ofstream out("datos.txt");
    for (int i = 0; i < Nmax; i++) {
        double rho = std::sqrt(y[0]*y[0] + a*a);  // plano ecuatorial
        out << phi << "," << y[0] << "," << rho*std::cos(phi)
            << "," << rho*std::sin(phi) << "\n";

        y = rk4_step(y, h);
        phi += h;

        if (y[0] < r_plus + 0.01) break;               // capturado
        if (y[0] > r_start && y[1] > 0) break;         // ya escapó
    }
    return 0;
}
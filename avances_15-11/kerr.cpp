#include <cmath>
#include <iostream>
#include <array>
#include <fstream>
#include <vector>
#include <random>
#include <functional>

const int N = 8000;
const double a=0.95;
const double b=7;


double funcion(double r);

template <std::size_t M>
void metodorungekutte(double ti, double tf, double h, int N, std::array<double, M>& r1, std::function<double(double)> func);

int main(){

    double PI = 3.14159265358979323846;
    double pi = PI;
    double pf = 2.0 * PI +PI;
    double h1 = (pf - pi) / N;

    std::array<double, N> r_sol;
    r_sol[0] = 12.0;

    metodorungekutte(pi, pf, h1, N, r_sol, funcion);
    

    std::ofstream archivo1("datos.txt");
    for (int i = 0; i < N; i++) {
        double phi = pi + i * h1;
        double x = std::sqrt((r_sol[i]*r_sol[i] + a*a)) * std::cos(phi);
        double y = std::sqrt((r_sol[i]*r_sol[i] + a*a)) * std::sin(phi);
        archivo1 << phi << "," << r_sol[i] << "," << x << "," << y << "\n";
    }
    archivo1.close();

    return 0; 
}

double funcion(double r){
    double dis = r*r + (2/r)*(a-b)*(a-b) + a*a - b*b;
    if (dis < 0.0){
        return 0.0;
    }
   return ((r*r + a*a-2*r)/(b*(r - 2) + 2*a))* std::sqrt(dis);
}

template <std::size_t M>
void metodorungekutte(double ti, double tf, double h, int N, std::array<double, M>& r1, std::function<double(double)> func){
    double k1=0.0;
    double k2=0.0;
    double k3=0.0;
    double k4=0.0;
    double sig = -1.0;
    for (int i = 0; i<(N-1); i++ ){
        k1 = h * func(r1[i]);
        k2 = h * (func(r1[i] + k1/2)); 
        k3 = h * (func(r1[i] + k2/2));
        k4 = h * (func(r1[i] + k3));
        r1[i+1] = r1[i] + sig*(k1 + 2*k2 + 2*k3 + k4)/6.0; 

        if (sig < 0.0 && r1[i+1] < 1e-6){
            sig = 1.0;
        }
    }
}

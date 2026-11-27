#include "auxiliares.hpp"
#include "cond_Inici.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include "/usr/include/eigen3/Eigen/Dense"

// [0]: r, [1]: theta, [2]: phi, [3]: pr, [4]: ptheta
Vector5d rkf45(Vector5d& estado, double del, double sig, double a, double& h, double Ei, double Lzi, double kappa, double tol, double h_max, double h_min);

// Modelo analítico rápido de emisividad + efecto Doppler en el disco de acreción ecuatorial
double evaluar_disco(double r, double a, double Ei, double Lzi, double r_isco, double r_disk_max) {
    if (r < r_isco || r > r_disk_max) return 0.0;

    // Velocidad angular kepleriana relativista en Kerr (disco prógrado)
    double omega = 1.0 / (std::pow(r, 1.5) + a);
    
    // Factor de corrimiento al rojo (gravitacional + Doppler) g = u^t_em^-1 / (1 - omega * Lz / E)
    double g_tt = -(1.0 - 2.0 / r);
    double g_tphi = -2.0 * a / r;
    double g_phiphi = r * r + a * a + (2.0 * a * a) / r;
    
    double norm = -(g_tt + 2.0 * omega * g_tphi + omega * omega * g_phiphi);
    if (norm <= 0.0) return 0.0;
    
    double ut = 1.0 / std::sqrt(norm);
    double b = Lzi / Ei; // Parámetro de impacto azimutal
    double g = 1.0 / (ut * (1.0 - omega * b)); // Factor g (corrimiento azul/rojo)

    // Perfil radial de temperatura tipo Novikov-Thorne / Shackura-Sunyaev simplificado
    double perfil_radial = std::pow(r_isco / r, 3.0) * (1.0 - std::sqrt(r_isco / r));
    if (perfil_radial < 0.0) perfil_radial = 0.0;

    // Intensidad observada escalada por g^4 (invariante relativista I_nu / nu^3)
    return perfil_radial * std::pow(g, 4.0) * 150.0;
}

int main() {
    // --- CONFIGURACIÓN DE RENDER SEMI-REALISTA ---
    const int width = 1280;
    const int height = 720;
    const int samples_per_pixel = 4; // 2x2 Antialiasing (Total: 3,686,400 fotones)
    const long long total_fotones = static_cast<long long>(width) * height * samples_per_pixel;

    // Parámetros físicos de Kerr
    const double a = 0.9;
    const double r_obs = 50.0; // Observador más lejano (más realista)
    const double theta_obs = 82.0 * M_PI / 180.0; // Inclinación casi de canto (82°) para ver el arco gravitacional
    const double r_horizonte = 1.0 + std::sqrt(1.0 - a * a);
    const double limite_horizonte = r_horizonte * 1.01;
    const double r_escape = 120.0;

    // Radio ISCO (Órbita Circular Estable Más Interna) aproximado para a = 0.9
    const double r_isco = 2.32; 
    const double r_disk_max = 22.0;

    // Parámetros del integrador RKF45
    const double h_inicial = 0.1;
    const double h_max = 2.0;
    const double h_min = 1e-6;
    const int max_pasos = 15000;
    const double tol = 1e-8;

    // Campo de visión (FOV) de la cámara
    const double aspect_ratio = static_cast<double>(width) / height;
    const double escala_fov = 0.45; // Ajusta el zoom sobre el agujero negro

    // Buffer de imagen RGB en memoria
    std::vector<double> framebuffer(width * height * 3, 0.0);

    long long pasos_totales = 0;
    long long fotones_capturados = 0;
    long long cruces_disco = 0;

    std::cout << "Renderizando " << width << "x" << height << " (" << samples_per_pixel << "x MSAA = " 
              << total_fotones << " fotones)...\n";

    auto inicio = std::chrono::high_resolution_clock::now();

    // Si compilas con -fopenmp, esta línea reparte los píxeles entre todos los núcleos de tu CPU
    #pragma omp parallel for schedule(dynamic, 4) reduction(+:pasos_totales, fotones_capturados, cruces_disco)
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            double r_acum = 0.0, g_acum = 0.0, b_acum = 0.0;

            // Sub-muestreo 2x2 por píxel (Antialiasing)
            for (int sy = 0; sy < 2; ++sy) {
                for (int sx = 0; sx < 2; ++sx) {
                    // Coordenadas normalizadas de pantalla con desplazamiento sub-píxel
                    double u = ((px + (sx + 0.5) * 0.5) / width) * 2.0 - 1.0;
                    double v = 1.0 - ((py + (sy + 0.5) * 0.5) / height) * 2.0;

                    double xp = u * aspect_ratio * escala_fov;
                    double yp = v * escala_fov;

                    Vector5d state(r_obs, theta_obs, 0.0, 0.0, 0.0);
                    double h = h_inicial;

                    double delta = aux::delta(state, a);
                    double sigma = aux::sigma(state, a);

                    std::array<double, 6> condi = coin::condini(state, xp, yp, sigma, delta, a);
                    state(3) = condi[4];
                    state(4) = condi[5];

                    double Ei = std::sqrt(aux::calcularEi2(state, condi, a, sigma, delta));
                    double Lzi = aux::calcular_Lzi(state, condi, a, Ei, sigma, delta);
                    double kappa = aux::calcular_kappa(state, Lzi, a, Ei);

                    double intensidad_rayo = 0.0;
                    double transmitancia = 1.0; // Para anillos secundarios semi-transparentes

                    for (int paso = 0; paso < max_pasos; ++paso) {
                        double theta_anterior = state(1);
                        double r_anterior = state(0);

                        delta = aux::delta(state, a);
                        sigma = aux::sigma(state, a);
                        state = rkf45(state, delta, sigma, a, h, Ei, Lzi, kappa, tol, h_max, h_min);
                        pasos_totales++;

                        double r_actual = state(0);
                        double theta_actual = state(1);

                        // 1. Verificar si cruzó el plano ecuatorial (theta = PI/2)
                        if ((theta_anterior - M_PI_2) * (theta_actual - M_PI_2) <= 0.0) {
                            // Interpolación lineal del radio de cruce
                            double t_cruce = std::abs(theta_anterior - M_PI_2) / 
                                             (std::abs(theta_anterior - M_PI_2) + std::abs(theta_actual - M_PI_2) + 1e-15);
                            double r_cruce = r_anterior + t_cruce * (r_actual - r_anterior);

                            double brillo = evaluar_disco(r_cruce, a, Ei, Lzi, r_isco, r_disk_max);
                            if (brillo > 0.0) {
                                intensidad_rayo += transmitancia * brillo;
                                transmitancia *= 0.45; // Atenuación al atravesar el plasma del disco
                                cruces_disco++;
                            }
                        }

                        // 2. Verificar captura en el horizonte de eventos
                        if (r_actual <= limite_horizonte) {
                            fotones_capturados++;
                            break;
                        }

                        // 3. Verificar escape a la esfera celeste
                        if (r_actual > r_escape) {
                            // Fondo estelar procedural simple cuando el rayo escapa al infinito
                            double hash_estrella = std::sin(state(1) * 127.1 + state(2) * 311.7) * 43758.5453;
                            hash_estrella -= std::floor(hash_estrella);
                            if (hash_estrella > 0.9985) {
                                intensidad_rayo += transmitancia * 0.8;
                            }
                            break;
                        }
                    }

                    // Mapeo de color térmico simple (negro -> naranja -> blanco azulado)
                    r_acum += std::min(1.0, intensidad_rayo * 1.4);
                    g_acum += std::min(1.0, intensidad_rayo * 0.85);
                    b_acum += std::min(1.0, intensidad_rayo * 0.45);
                }
            }

            int idx = (py * width + px) * 3;
            framebuffer[idx]     = r_acum / samples_per_pixel;
            framebuffer[idx + 1] = g_acum / samples_per_pixel;
            framebuffer[idx + 2] = b_acum / samples_per_pixel;
        }
    }

    auto fin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracion = fin - inicio;

    // --- VOLCADO DE IMAGEN A DISCO (.PPM) ---
    auto inicio_io = std::chrono::high_resolution_clock::now();
    std::ofstream img("kerr_render.ppm", std::ios::binary);
    img << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height * 3; ++i) {
        // Corrección gamma (2.2) para visualización realista
        double val = std::pow(std::max(0.0, std::min(1.0, framebuffer[i])), 1.0 / 2.2);
        unsigned char byte = static_cast<unsigned char>(val * 255.0);
        img.write(reinterpret_cast<char*>(&byte), 1);
    }
    img.close();
    auto fin_io = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracion_io = fin_io - inicio_io;

    std::cout << "==================================================\n";
    std::cout << "Tiempo de trazado (CPU):  " << duracion.count() << " s\n";
    std::cout << "Tiempo guardar imagen:    " << duracion_io.count() << " s\n";
    std::cout << "Rendimiento promedio:     " << (total_fotones / duracion.count()) / 1000.0 << " kFotones/s\n";
    std::cout << "Pasos RKF45 evaluados:    " << pasos_totales << "\n";
    std::cout << "Impactos en el disco:     " << cruces_disco << "\n";
    std::cout << "Fotones en horizonte:     " << fotones_capturados << "\n";
    std::cout << "Imagen guardada en:       kerr_render.ppm\n";

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

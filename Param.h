#include <vector>
#include <cmath>

#ifndef PARAM_H
#define PARAM_H

/**
 * \brief Structure contenant les paramètres de la simulation
 */
struct Param{

    // Paramètres physiques

    // Taille de domaine
    double xmin=-0.5; 
    double xmax=0.5;
    double ymin=-0.5; 
    double ymax=0.5; 

    double kappa = 1e-6;    // coefficient de diffusion constant
    double Tmax = 1350;     // amplitude
    double sigma = 0.3;     // ecart-type de la solution gaussienne
    double t_final = 1e4;   // temps auquel est calculé T


    // Paramètres numériques
    std::vector<int> sizes = {64, 128, 256, 512};   // tailles des matrices à tester
    double dx, dy, dt;                              // discrétisation de l'esapce et du temps
    int Nt;                                         // nombre d'itérations temporelles

    /**
     * \brief Calcul des paramètres numériques selon la taille du système
     * 
     * \param size taille du carré
     */
    void update(int size){
        dx = (xmax-xmin)/size;
        dy = (ymax-ymin)/size;
        dt = 1./(4 * kappa * (1/pow(dx,2) + 1/pow(dy,2)));  // dt vaut 1/2 fois la condition CFL
        Nt = (int)(t_final/dt);                             // nombre d'itérations pour les différences finies
    }
};

#endif
#include <vector>
#include <cmath>
#include <nlohmann/json.hpp>
#include <fstream>

#ifndef PARAM_H
#define PARAM_H

using json = nlohmann::json;

/**
 * \brief Structure containing the simulation parameters
 */
struct Param{

    // Physical parameters

    // domain size
    double xmin=-0.5;   ///< minimum x coordinate
    double xmax=0.5;    ///< maximum x coordinate
    double ymin=-0.5;   ///< minimum y coordinate
    double ymax=0.5;    ///< maximum y coordinate

    double kappa;       ///< constant diffusion coefficient
    double Tmax;        ///< amplitude
    double sigma;       ///< standard deviation of the Gaussian solution
    double t_final;     ///< time at which T is calculated


    // Numerical parameters

    std::vector<int> sizes; ///< sizes of matrices to be tested
    double dx, dy, dt;      ///< spatial and temporal steps
    int Nt;                 ///< number of iterations for finite differences

    /**
     * @brief   Load parameters from a JSON file
     * 
     * @param filename 
     */
    void load(const std::string& filename);

    /**
     * @brief   Update numerical parameters according to system size
     * 
     * @param size 
     */
    void update(int size);
};

#endif
#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] is use energy equation
     * @param[in] initial temperature
     * @param[in] Prantl number
     */
    Fields(const Grid &grid, double _nu, double _dt, double _tau, int imax, int jmax, double UI, double VI, double PI,
           double Gx, double Gy, bool energy_on, double TI, double PR, double beta);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */
    void calculate_fluxes(Grid &grid);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid);

    /**
     * @brief Temperature calculation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_temperature(Grid &grid);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressre index based access and modify
    double &p(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);

    /// temperature index based access and modify
    double &T(int i, int j);

    /// get timestep size
    double dt() const;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

    /// get _energy_on
    bool energy_on();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// x-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;

    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;

    /// temperature matrix
    Matrix<double> _T;

    /// kinematic viscosity
    double _nu;
    /// gravitional accelearation in x direction
    double _gx{0.0};
    /// gravitional accelearation in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    /// Prandtl number
    double _Pr;
    /// thermal expansion coefficient
    double _beta;
    /// indicate whether apply temperature or not
    bool _energy_on;
};

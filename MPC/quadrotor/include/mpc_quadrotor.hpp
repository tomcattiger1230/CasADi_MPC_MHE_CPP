/*
 * @Author: Wei Luo
 * @Date: 2023-01-06 19:47:07
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-13 23:27:16
 * @Note: Note
 */
#ifndef __MPC_QUADROTOR__
#define __MPC_QUADROTOR__
#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <chrono>
#include <iostream>
#include <math.h>

namespace ca = casadi;


class MPCQuadrotor {
public:
  MPCQuadrotor(const double dt, const double N,
               const double gravity_acceleration, const double roll_tau,
               const double roll_gain, const double pitch_tau,
               const double pitch_gain);
  virtual ~MPCQuadrotor();
  void initialization_formulation();
  void set_boundary(const std::vector<double> u_min, const std::vector<double> u_max, const std::vector<double> x_min, const std::vector<double> x_max);
  void get_results(std::vector<double> init_value,
                   std::vector<double> desired_trajectory);
  ca::Function get_system_dynamics()
  {
    return system_dynamics_;
  };

  int num_controls;
  int num_states;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  double g_acceleration_;
  int prediction_horizon_;
  double dt_;
  ca::Slice slice_all = ca::Slice();
  ca::Slice slice_state = ca::Slice(0, 6);

  // control gain
  double roll_tau_;
  double roll_gain_;
  double pitch_gain_;
  double pitch_tau_;

  // system dynamics
  ca::Function system_dynamics_;

  ca::Function solver_;

  // boundary conditions
  std::vector<double> u_min_;
  std::vector<double> u_max_;
  std::vector<double> x_min_;
  std::vector<double> x_max_;

  std::map<std::string, casadi::DM> opt_results_;

  template <typename T1, typename T2>
  T1 RK4_function(T1 state, T1 control, T2 dynamics)
  {
    std::vector<T1> input;
    input.push_back(state);
    input.push_back(control);
    T1 k1 = dynamics(input).at(0);
    input[0] = state + dt_ / 2.0 * k1;
    T1 k2 = dynamics(input).at(0);
    input[0] = state + dt_ / 2.0 * k2;
    T1 k3 = dynamics(input).at(0);
    input[0] = state + dt_ * k3;
    T1 k4 = dynamics(input).at(0);
    return state + dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
};
#endif
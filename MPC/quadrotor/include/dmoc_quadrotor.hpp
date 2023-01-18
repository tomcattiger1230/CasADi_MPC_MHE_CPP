/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:33:44
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-18 17:49:06
 * @Note: Note
 */

#ifndef __DMOC_QUADROTOR__
#define __DMOC_QUADROTOR__
#include "derivation_uav.hpp"
#include <chrono>
#include <iostream>

class DMOCUAV : public DerivationUAV {
public:
  DMOCUAV(const double dt, const double N, const double mass_quadrotor,
          const std::vector<double> inertia_moment, const double frame_size,
          const double motor_torque_const, const double g);

  void initialization_formulation();
  void
  set_boundary(const std::vector<double> u_min, const std::vector<double> u_max,
               const std::vector<double> x_min,
               const std::vector<double> x_max);
  void get_results(std::vector<double> init_value,
                   std::vector<double> desired_trajectory,
                   Eigen::MatrixXd &result_x_matrix,
                   Eigen::MatrixXd &result_u_matrix);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int prediction_horizon_;
  double dt_;
  ca::Slice slice_state_ = ca::Slice(0, 6);
  ca::Slice slice_all = ca::Slice();
  ca::Function solver_;

  // boundary conditions
  std::vector<double> u_min_;
  std::vector<double> u_max_;
  std::vector<double> x_min_;
  std::vector<double> x_max_;

  std::map<std::string, casadi::DM> opt_results_;
};

#endif
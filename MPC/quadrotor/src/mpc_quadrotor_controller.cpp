/*
 * @Author: Wei Luo
 * @Date: 2023-01-06 21:14:04
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-14 21:21:26
 * @Note: Note
 */

#include "mpc_quadrotor.hpp"
#include <Eigen/Dense>

int main() {
  constexpr double quadrotor_mass = 1.659;
  constexpr double g_acceleration = 9.8066;
  double roll_tau = 0.257;
  double roll_gain = 0.75;
  double pitch_tau = 0.259;
  double pitch_gain = 0.78;
  auto mpc_quadrotor_handle = std::make_shared<MPCQuadrotor>(MPCQuadrotor(
      0.33, 20, g_acceleration, roll_tau, roll_gain, pitch_tau, pitch_gain));
  mpc_quadrotor_handle->initialization_formulation();

  Eigen::MatrixXd init_trajectory(9, 20);
  init_trajectory.block<9, 1>(0, 0) << 0, 0, 0.2, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 1) << 0, 0, 0.4, 0, 0, 0.01, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 2) << 0, 0, 0.5, 0, 0, 0.0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 3) << 0.1, 0.0, 0.6, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 4) << 0.1, 0, 0.67, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 5) << 0.1, 0, 0.69, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 6) << 0.2, 0.0, 0.73, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 7) << 0.2, 0.0, 0.76, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 8) << 0.2, 0.0, 0.8, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 9) << 0.2, 0.0, 0.83, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 10) << 0.2, 0.0, 0.85, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 11) << 0.3, 0.0, 0.88, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 12) << 0.3, 0.2, 0.91, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 13) << 0.4, 0.2, 0.93, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 14) << 0.5, 0.2, 0.95, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 15) << 0.5, 0.2, 0.97, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 16) << 0.7, 0.2, 0.99, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 17) << 0.8, 0.2, 1.0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 18) << 0.9, 0.2, 1.0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<9, 1>(0, 19) << 1.0, 0.2, 1.0, 0, 0, 0, 0, 0, 0;

  std::vector<double> control_lb = {-M_PI / 4.0, -M_PI / 4.0,
                                    g_acceleration * 0.5};
  std::vector<double> control_ub = {M_PI / 4.0, M_PI / 4.0,
                                    g_acceleration * 1.5};
  std::vector<double> state_lb = {-ca::inf, -ca::inf, -ca::inf,
                                  -ca::inf, -ca::inf, -ca::inf,
                                  -ca::inf, -ca::inf, -ca::inf};
  std::vector<double> state_ub = {ca::inf, ca::inf, ca::inf, ca::inf, ca::inf,
                                  ca::inf, ca::inf, ca::inf, ca::inf};

  mpc_quadrotor_handle->set_boundary(control_lb, control_ub, state_lb,
                                     state_ub);
  int iter_index = 0;
  const int max_iter = 1;
  Eigen::MatrixXd opt_x = Eigen::MatrixXd::Zero(9, 20);
  Eigen::MatrixXd opt_u = Eigen::MatrixXd::Zero(3, 19);
  Eigen::MatrixXd ref_trajectory = init_trajectory;

  while (iter_index < max_iter)
  {
    std::vector<double> init_state_std(opt_x.data(),
                                       opt_x.size() + opt_x.data());
    std::vector<double> init_control_std(opt_u.data(),
                                         opt_u.size() + opt_u.data());
    std::vector<double> opt_init;
    opt_init.insert(opt_init.end(), init_state_std.begin(), init_state_std.end());
    opt_init.insert(opt_init.end(), init_control_std.begin(), init_control_std.end());
    std::vector<double> ref_trajectory_std(
        ref_trajectory.data(), ref_trajectory.data() + ref_trajectory.size());
    mpc_quadrotor_handle->get_results(opt_init, ref_trajectory_std, opt_x, opt_u);
    std::cout << opt_x << std::endl;
    iter_index += 1;
  }

  return 0;
}
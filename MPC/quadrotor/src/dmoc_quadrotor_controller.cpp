/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:35:40
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-02-09 19:12:09
 * @Note: Note
 */

#include "dmoc_quadrotor.hpp"

void circle_trajectory(Eigen::VectorXd current_state,
                       Eigen::MatrixXd &trajectory, int iter,
                       int trajectory_horizon = 20) {
  if (iter <= 30) {
    if (current_state[2] >= 0.2) {
      trajectory(Eigen::all, Eigen::seqN(0, trajectory.cols() - 1)) =
          trajectory(Eigen::all, Eigen::seqN(1, trajectory.cols()));
    }
  } else {
    Eigen::VectorXd idx = Eigen::VectorXd::LinSpaced(19, 0.0, 18.0);
    trajectory(Eigen::all, Eigen::seqN(0, trajectory.cols() - 1)) =
        trajectory(Eigen::all, Eigen::seqN(1, trajectory.cols()));
    for (int i = 1; i < trajectory_horizon; i++) {
      double temp = ((double)i - 1.0 + (double)iter - 30.0) / 360.0 * M_PI;
      trajectory.block<2, 1>(0, i) << std::cos(temp), std::sin(temp);
    }
  }
}

int main() {
  const double quadrotor_mass = 1.659;
  const double g_acceleration = 9.8066;
  const double dt = 0.33;
  const int horizon = 20;
  std::vector<double> inertia = {0.01576, 0.01540, 0.01948};
  const double frame_size = 0.33;
  const double motor_torque_const = 0.013;

  auto dmoc_uav_handle = std::make_shared<DMOCUAV>(
      DMOCUAV(dt, horizon, quadrotor_mass, inertia, frame_size,
              motor_torque_const, g_acceleration));

  dmoc_uav_handle->initialization_formulation();

  Eigen::VectorXd current_state(12);
  current_state << 0.0, 0.0, 0.1, 0, 0, 0., 0.0, 0.0, 0, 0, 0, 0;
  Eigen::MatrixXd init_trajectory(12, 20);
  init_trajectory.block<12, 1>(0, 0) << 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<12, 1>(0, 1) << 0, 0, 0.2, 0, 0, 0.0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<12, 1>(0, 2) << 0, 0, 0.3, 0, 0, 0.0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<12, 1>(0, 3) << 0.1, 0.0, 0.4, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 4) << 0.1, 0, 0.57, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<12, 1>(0, 5) << 0.1, 0, 0.69, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  init_trajectory.block<12, 1>(0, 6) << 0.2, 0.0, 0.73, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 7) << 0.2, 0.0, 0.76, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 8) << 0.2, 0.0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 9) << 0.2, 0.0, 0.83, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 10) << 0.2, 0.0, 0.85, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 11) << 0.3, 0.0, 0.88, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 12) << 0.3, 0.2, 0.91, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 13) << 0.4, 0.2, 0.93, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 14) << 0.5, 0.2, 0.95, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 15) << 0.5, 0.2, 0.97, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 16) << 0.7, 0.2, 0.99, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 17) << 0.8, 0.2, 1.0, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 18) << 0.9, 0.2, 1.0, 0, 0, 0, 0, 0, 0, 0, 0,
      0;
  init_trajectory.block<12, 1>(0, 19) << 1.0, 0.2, 1.0, 0, 0, 0, 0, 0, 0, 0, 0,
      0;

  const double hover_force = g_acceleration * quadrotor_mass * 0.25;
  std::vector<double> control_lb = {hover_force * 0.5, hover_force * 0.5,
                                    hover_force * 0.5, hover_force * 0.5};
  std::vector<double> control_ub = {hover_force * 1.5, hover_force * 1.5,
                                    hover_force * 1.5, hover_force * 1.5};
  std::vector<double> state_lb = {-ca::inf, -ca::inf, 0.0, -M_PI, -M_PI, -M_PI};
  std::vector<double> state_ub = {ca::inf, ca::inf, ca::inf, M_PI, M_PI, M_PI};
  std::vector<double> v_limit{ca::inf, ca::inf, ca::inf};
  std::vector<double> d_rpy{ca::inf, ca::inf, ca::inf};

  dmoc_uav_handle->set_boundary(control_lb, control_ub, state_lb, state_ub,
                                v_limit, d_rpy);

  int iter_index = 0;
  const int max_iter = 100;
  Eigen::MatrixXd opt_x = Eigen::MatrixXd::Zero(6, 20);
  Eigen::MatrixXd opt_u = Eigen::MatrixXd::Zero(4, 20);
  opt_u = Eigen::MatrixXd::Constant(4, 20, hover_force);
  Eigen::MatrixXd ref_trajectory = init_trajectory;
  std::vector<Eigen::VectorXd> real_trajectory_std;
  real_trajectory_std.push_back(current_state);
  std::vector<double> parameter_std;

  while (iter_index < max_iter) {
    ref_trajectory.block<12, 1>(0, 0) = current_state;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<double> init_state_std(opt_x.data(),
                                       opt_x.size() + opt_x.data());
    std::vector<double> init_control_std(opt_u.data(),
                                         opt_u.size() + opt_u.data());
    std::vector<double> opt_init;
    opt_init.insert(opt_init.end(), init_state_std.begin(),
                    init_state_std.end());
    opt_init.insert(opt_init.end(), init_control_std.begin(),
                    init_control_std.end());

    // std::vector<double> current_state_std(
    //     current_state.data(), current_state.data() + current_state.size());
    std::vector<double> ref_trajectory_std(
        ref_trajectory.data(), ref_trajectory.data() + ref_trajectory.size());
    // std::cout << "reference trajectory" << ref_trajectory.transpose() <<
    // std::endl; std::merge(current_state.begin(), current_state.end(),
    //            ref_trajectory_std.begin(),
    //            ref_trajectory_std.end(),
    //            std::back_inserter(parameter_std));
    dmoc_uav_handle->get_results(opt_init, ref_trajectory_std, opt_x, opt_u);
    // std::cout << opt_x.transpose() << std::endl;
    // std::cout << opt_u.transpose() << std::endl;
    dmoc_uav_handle->model_based_movement(current_state, opt_u.col(0), opt_x,
                                          opt_u);
    // std::cout << "robot's current state: " << current_state.transpose() <<
    // std::endl;
    auto stop_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        stop_time - start_time);
    std::cout << "average calculation time for each iteration [s]: "
              << duration.count() / 1e6 << std::endl;
    real_trajectory_std.push_back(current_state);
    circle_trajectory(current_state, ref_trajectory, iter_index);
    // next iteration
    iter_index += 1;
  }
  for (Eigen::VectorXd i : real_trajectory_std) {
    std::cout << i.transpose() << std::endl;
  }
}
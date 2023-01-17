/*
 * @Author: Wei Luo
 * @Date: 2022-12-27 22:50:51
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-17 16:52:44
 * @Note: Note
 */

#include <dmcc_uav_manipulator.hpp>

casadi::Function fixed_ref_position_function(const std::vector<double> x) {
  auto dt = casadi::MX::sym("dt");
  auto index = casadi::MX::sym("index");
  casadi::DM pos{x};
  return casadi::Function("fixed_ref_function", {dt, index}, {x});
}

casadi::Function fixed_ref_velocity_function() {
  auto dt = casadi::MX::sym("dt");
  auto index = casadi::MX::sym("index");
  auto v = casadi::DM::zeros(3);
  return casadi::Function("fixed_ref_function", {dt, index}, {v});
}

int main() {
  const double quadrotor_mass = 1.659;
  const double manipulator_mass = 0.36113;
  const double g = 9.8066;
  auto dmcc_uavm_handle =
      std::make_shared<DMCCUAVManipulator>(DMCCUAVManipulator(
          quadrotor_mass, manipulator_mass, {0.01576, 0.01540, 0.01948},
          {6.1324e-5, 0.00163814, 0.00162087}, 0.34, true, {0.0, 0.0, -0.107},
          0.33, 0.013, g));

  // an example
  Eigen::VectorXd init_pose(14);
  init_pose << -0.8, 0.0, 0.7, 0.0, 0.0, 0.0, M_PI / 2.0, 0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0;
  Eigen::MatrixXd designed_waypoints(14, 2);
  designed_waypoints.block<14, 1>(0, 0) << 0.0, 0.5, 0.54, 0.0, 0.0, 0.0,
      M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

  designed_waypoints.block<14, 1>(0, 1) << 1.0, 0.0, 0.7, 0.0, 0.0, 0.0,
      M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  std::cout << "designed waypoints :" << std::endl
            << designed_waypoints << std::endl;
  const std::vector<double> position_of_target = {0.0, 0.5, 0.28};
  dmcc_uavm_handle->get_target_position_function(
      fixed_ref_position_function(position_of_target));
  dmcc_uavm_handle->get_target_velocity_function(fixed_ref_velocity_function());

  dmcc_uavm_handle->initialization_formulation(
      designed_waypoints.cols(), 15); // num_def_points, num_pred_points
  dmcc_uavm_handle->get_path_waypoints(6, init_pose, designed_waypoints, 0.3);

  const double hover_force = g * (quadrotor_mass + manipulator_mass) * 0.25;
    std::vector<double> v_limit{0.18, 0.18, 0.12};
    std::vector<double> d_rpy_alpha{0.2, 0.2, 0.1, M_PI / 5.0};
    std::vector<double> min_x {-0.5, -1.0, 0.0, -M_PI / 4.0, -M_PI / 4.0,
    -M_PI, M_PI / 3.0};
    std::vector<double> max_x {1.8, 1.0, 1.5, M_PI / 4.0, M_PI / 4.0, M_PI,
                              2.0 * M_PI / 3.0};
    std::vector<double> min_u {0.5 * hover_force, 0.5 * hover_force,
                              0.5 * hover_force, 0.5 * hover_force, -2.0};
    std::vector<double> max_u {1.5 * hover_force, 1.5 * hover_force,
                              1.5 * hover_force, 1.5 * hover_force, 2.0};
    const double max_accept_distance = 0.02;
    const double max_accept_velocity = 0.01 * 0.01;
    const double max_accept_heading = 0.01 * 0.01;
    const double high_interation_threshold =
        -0.01;

    dmcc_uavm_handle->set_constraints(
        min_x, max_x, min_u, max_u, max_accept_distance, v_limit,
        d_rpy_alpha, max_accept_velocity, max_accept_heading,
        high_interation_threshold);
      dmcc_uavm_handle->get_results();
    return 0;
}

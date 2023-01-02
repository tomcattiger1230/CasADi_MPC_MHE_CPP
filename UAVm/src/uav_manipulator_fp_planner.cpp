/*
 * @Author: Wei Luo
 * @Date: 2022-12-27 22:50:51
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-02 23:01:46
 * @Note: Note
 */

#include <dmcc_uav_manipulator.hpp>

casadi::Function fixed_ref_position_function() {
  auto dt = casadi::MX::sym("dt");
  auto index = casadi::MX::sym("index");
  auto x = casadi::MX::sym("x", 3);
  return casadi::Function("fixed_ref_function", {dt, index}, {x});
}

casadi::Function fixed_ref_velocity_function() {
  auto dt = casadi::MX::sym("dt");
  auto index = casadi::MX::sym("index");
  auto v = casadi::DM::zeros(3);
  return casadi::Function("fixed_ref_function", {dt, index}, {v});
}

int main() {
  auto dmcc_uavm_handle = std::make_shared<DMCCUAVManipulator>(
      DMCCUAVManipulator(1.659, 0.36113, {0.01576, 0.01540, 0.01948},
                         {6.1324e-5, 0.00163814, 0.00162087}, 0.34, true,
                         {0.0, 0.0, -0.107}, 0.33, 0.013, 9.8066));

  // an example
  Eigen::VectorXd init_pose(14);
  init_pose << 0.0, 0.0, 0.65, 0.0, 0.0, 0.0, M_PI / 2.0, 0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0;
  Eigen::MatrixXd designed_waypoints(2, 14);
  designed_waypoints.block<1, 14>(0, 0) << 1.0, 0.0, 0.54, 0.0, 0.0, 0.0,
      M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

  designed_waypoints.block<1, 14>(1, 0) << 2.5, 0.0, 0.65, 0.0, 0.0, 0.0,
      M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  std::cout << "designed waypoints :" << std::endl
            << designed_waypoints << std::endl;

  dmcc_uavm_handle->get_target_position_function(fixed_ref_position_function());
  dmcc_uavm_handle->get_target_velocity_function(fixed_ref_velocity_function());

  dmcc_uavm_handle->initialization_formulation(
      designed_waypoints.rows(), 15); // num_def_points, num_pred_points
  dmcc_uavm_handle->get_path_waypoints(6, init_pose, designed_waypoints, 0.3);

  //   double hover_force = 9.8066 * (1.659 + 0.36113) * 0.25;
  //   std::vector<double> v_limit{0.18, 0.18, 0.12};
  //   std::vector<double> d_rpy_alpha{0.2, 0.2, 0.1, M_PI / 5.0};
  //   std::vector<double> min_x {-0.5, -1.0, 0.0, -M_PI / 4.0, -M_PI / 4.0,
  //   -M_PI,
  //                             M_PI / 3.0};
  //   std::vector<double> max_x {1.8, 1.0, 1.5, M_PI / 4.0, M_PI / 4.0, M_PI,
  //                             2.0 * M_PI / 3.0};
  //   std::vector<double> min_u {0.5 * hover_force, 0.5 * hover_force,
  //                             0.5 * hover_force, 0.5 * hover_force, -2.0};
  //   std::vector<double> max_u {1.5 * hover_force, 1.5 * hover_force,
  //                             1.5 * hover_force, 1.5 * hover_force, 2.0};
  //   std::vector<double> lambda_range {0.0, 1.0};

  //   std::vector<double> nu_range {0.0, 0.05};

  //   dmcc_uavm_handle->set_constraints(v_limit, d_rpy_alpha, min_x, max_x,
  //   min_u,
  //                                     max_u, lambda_range, nu_range);
  //   dmcc_uavm_handle->get_results();
  return 0;
}

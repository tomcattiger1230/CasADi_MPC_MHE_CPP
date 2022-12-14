/*
 * @Author: Wei Luo
 * @Date: 2022-11-30 15:48:33
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-04 17:41:31
 * @Note: Note
 */

#include <dmcc_uav_manipulator.hpp>

int main() {
    const double quadrotor_mass = 1.659;
    const double manipulator_mass = 0.36113;
    const double g_acceleration = 9.8066;
    auto dmcc_uavm_handle =
        std::make_shared<DMCCUAVManipulator>(DMCCUAVManipulator(
            quadrotor_mass, manipulator_mass, {0.01576, 0.01540, 0.01948},
            {6.1324e-5, 0.00163814, 0.00162087}, 0.34, false,
            {0.0, 0.0, -0.107}, 0.33, 0.013, g_acceleration));

    // an example
    Eigen::VectorXd init_pose(14);
    init_pose << 0.0, 0.0, 0.55, 0.0, 0.0, 0.0, M_PI / 2.0, 0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0;
    Eigen::MatrixXd designed_waypoints(14, 3);
    //   designed_waypoints.block<1, 14>(0, 0) << 0.0, 0.5, 0.5623, 0.0, 0.0, 0.0,
    //       M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    designed_waypoints.block<14, 1>(0, 0) << 0.55, 0.2, 0.4, 0.0, 0.0, 0.0,
        M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    designed_waypoints.block<14, 1>(0, 1) << 0.95, -0.1, 0.5, 0.0, 0.0, 0.0,
        M_PI / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    designed_waypoints.block<14, 1>(0, 2) << 1.5, 0, 0.55, 0.0, 0.0, 0.0,
        M_PI / 2.0, .0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    std::cout << designed_waypoints << std::endl;
    dmcc_uavm_handle->initialization_formulation(designed_waypoints.cols(),
                                                 15); // num_def_points, num_pred_points
    dmcc_uavm_handle->get_path_waypoints(init_pose, designed_waypoints, 0.3);

    const double hover_force =
        g_acceleration * (quadrotor_mass + manipulator_mass) * 0.25;
      std:: vector<double> v_limit{0.25, 0.25, 0.12};
      std::vector<double> d_rpy_alpha{0.2, 0.2, 0.1, M_PI / 5.0};
      std::vector<double> min_x {-0.5, -1.0, 0.0, -M_PI / 4.0, -M_PI / 4.0, -M_PI,
                                M_PI / 3.0};
      std::vector<double> max_x {1.8, 1.0, 1.5, M_PI / 4.0, M_PI / 4.0, M_PI,
                                2.0 * M_PI / 3.0};
      std::vector<double> min_u {0.5 * hover_force, 0.5 * hover_force,
                                0.5 * hover_force, 0.5 * hover_force, -2.0};
      std::vector<double> max_u {1.5 * hover_force, 1.5 * hover_force,
                                1.5 * hover_force, 1.5 * hover_force, 2.0};
      std::vector<double> lambda_range {0.0, 1.0};

      std::vector<double> nu_range {0.0, 0.05};

      dmcc_uavm_handle->set_constraints(v_limit, d_rpy_alpha, min_x, max_x, min_u,
                                        max_u, lambda_range, nu_range);
      dmcc_uavm_handle->get_results();
    return 0;
}

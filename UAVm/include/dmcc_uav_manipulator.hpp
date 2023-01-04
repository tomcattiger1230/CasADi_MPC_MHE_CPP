/*
 * @Author: Wei Luo
 * @Date: 2022-12-12 17:54:20
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-04 19:44:48
 * @Note: Note
 */

#ifndef __DMCC_UAV_MANIPULATOR_HPP__
#define __DMCC_UAV_MANIPULATOR_HPP__

#include <chrono>
#include <derivation_uav_manipulator_pos_force.hpp>
#include <iostream>

class DMCCUAVManipulator : public DerivationUAVm {
public:
  DMCCUAVManipulator(const double mass_quadrotor, const double mass_manipulator,
                     const std::vector<double> inertia_moment,
                     const std::vector<double> manipulator_inertia_moment,
                     const double manipulator_length,
                     const bool has_contact_target,
                     const std::vector<double> montage_offset_b,
                     const double frame_size, const double motor_torque_const,
                     const double g = 9.8066);
  ~DMCCUAVManipulator();


  void get_path_waypoints(const Eigen::VectorXd current_pose,
                          const Eigen::MatrixXd path_waypoints,
                          const double guessed_velocity);
  void get_path_waypoints(const int kappa0, const Eigen::VectorXd current_pose,
                          const Eigen::MatrixXd path_waypoints,
                          const double guessed_velocity);
  void initialization_formulation(const int num_def_waypoints,
                                  const int num_pred_waypoints);
  void set_constraints(const std::vector<double> v,
                       const std::vector<double>
                           d_rpy,
                       const std::vector<double>
                           min_state,
                       const std::vector<double>
                           max_state,
                       const std::vector<double>
                           min_control_input,
                       const std::vector<double>
                           max_control_input,
                       const std::vector<double>
                           lambda,
                       std::vector<double> nu);
  void set_constraints(const std::vector<double> min_state,
                       const std::vector<double> max_state,
                       const std::vector<double> min_control_input,
                       const std::vector<double> max_control_input,
                       const double acceptable_distance,
                       const std::vector<double> v,
                       const std::vector<double> d_rpy,
                       const double acceptable_velocity_difference, const double acceptable_heading_difference,
                       const double high_interation);

  void get_target_position_function(ca::Function function){
    target_position_function_ = function;
  }
  void get_target_velocity_function(ca::Function function) {
    target_velocity_function_ = function;
  }
  void get_results();
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      private : bool has_contact_target_;
  int kappa0_;
  int num_pred_waypoints_;
  int num_def_waypoints_;
  std::vector<int> i_switch_;
  double vel_guess_;
  double time_guess_;
  std::vector<double> init_values_;

  Eigen::VectorXd init_pose_;
  Eigen::MatrixXd path_waypoints_;

  ca::MX dt_;
  ca::MX obj_function_ = 0; // cost function
  ca::MX constraint_vector_; // constraint vector
  std::vector<ca::MX> constraint_vector_std_; // constraint
  void init_state_guess(const bool with_additional_relaxation = false);
  ca::Slice all = ca::Slice();
  ca::Function solver_;

  int number_equality_constraint_;

  // constraints
  std::vector<double> lbg_;
  std::vector<double> ubg_;
  std::vector<double> lbx_;
  std::vector<double> ubx_;
  std::vector<double> optimization_param_;
  std::map<std::string, casadi::DM> opt_results_;

  // without the contact target
  ca::MX lambda_param_;
  ca::MX tolerance_param_;

  // with the contact target
  ca::Function target_position_function_;
  ca::Function target_velocity_function_;
  ca::MX epsilon_param_;
  ca::MX kappa_param_;
  ca::MX contact_relax_param_;
};

#endif
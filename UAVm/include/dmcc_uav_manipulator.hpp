/*
 * @Author: Wei Luo
 * @Date: 2022-12-12 17:54:20
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-12-27 16:17:06
 * @Note: Note
 */

#ifndef __DMCC_UAV_MANIPULATOR_HPP__
#define __DMCC_UAV_MANIPULATOR_HPP__

#include <chrono>
#include <derivation_uav_manipulator_pos_force.hpp>

class DMCCUAVManipulator : public DerivationUAV {
public:
  DMCCUAVManipulator(
      const double manipulator_length, const double mass_quadrotor,
      const double mass_manipulator, const std::vector<double> inertia_moment,
      const std::vector<double> manipulator_inertia_moment,
      const bool has_contact_target, const bool has_manipulator,
      const std::vector<double> montage_offset_b, const double frame_size,
      const double motor_torque_const, const double g = 9.8066);
  ~DMCCUAVManipulator();

  void initialization_formulation(const int kappa0, const int num_def_waypoints,
                                  const int num_pred_waypoints);
  void get_path_waypoints(const Eigen::VectorXd current_pose,
                          const Eigen::MatrixXd path_waypoints,
                          const double guessed_velocity);
  void set_constraints(const std::vector<double> v,
                       const std::vector<double> d_rpy,
                       const std::vector<double> min_state,
                       const std::vector<double> max_state,
                       const std::vector<double> min_control_input,
                       const std::vector<double> max_control_input,
                       const std::vector<double> lambda,
                       std::vector<double> nu);
  void set_constraints(
      const std::vector<double> v, const std::vector<double> d_rpy,
      const std::vector<double> tolerance_param,
      const std::vector<double> high_interation);
  void get_results();
  ca::Function solver_;

private:
  bool has_contact_target_;
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
  void init_state_guess(const bool with_additional_relaxation = false);
  ca::Slice all;

  ca::MX lambda_param_;
  ca::MX tolerance_param_;
  int number_equal_constraint;

  // constraints
  std::vector<double> lbg_;
  std::vector<double> ubg_;
  std::vector<double> lbx_;
  std::vector<double> ubx_;
  std::vector<double> optimization_param_;
  std::map<std::string, casadi::DM> opt_results_;
};

#endif
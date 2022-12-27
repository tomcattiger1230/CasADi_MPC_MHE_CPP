/*
 * @Author: Wei Luo
 * @Date: 2022-12-14 11:00:34
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-12-27 16:50:32
 * @Note: Note
 */

#include <dmcc_uav_manipulator.hpp>

DMCCUAVManipulator::DMCCUAVManipulator(
    const double manipulator_length, const double mass_quadrotor,
    const double mass_manipulator, const std::vector<double> inertia_moment,
    const std::vector<double> manipulator_inertia_moment,
    const bool has_contact_target, const bool has_manipulator,
    const std::vector<double> montage_offset_b, const double frame_size,
    const double motor_torque_const, const double g)
    : DerivationUAV(g, manipulator_length, mass_manipulator, mass_quadrotor,
                    inertia_moment, manipulator_inertia_moment, frame_size,
                    motor_torque_const, montage_offset_b) {

  // only a trajectory planning or considering a specific mission
  has_contact_target_ = has_contact_target;
};
DMCCUAVManipulator::~DMCCUAVManipulator(){};

void DMCCUAVManipulator::get_path_waypoints(
    const Eigen::VectorXd current_pose, const Eigen::MatrixXd path_waypoints,
    const double guessed_velocity) {
  init_pose_ = current_pose;
  path_waypoints_ = path_waypoints;
  if (num_def_waypoints_ == path_waypoints.rows())
    std::cout << "num of defined path waypoints: " << num_def_waypoints_
              << std::endl;
  else {
    std::cout << "num of defined path waypoints is not equal to that is "
                 "initalized !!!"
              << std::endl;
  }

  std::vector<double> dist;
  dist.push_back((init_pose_.block<3, 1>(0, 0).transpose() -
                  path_waypoints_.block<1, 3>(0, 0))
                     .norm());

  if (dist[0] < 0.01) {
    dist[0] += 0.01; // incase of 0 or almost zero distance
  }

  for (int i = 1; i < num_def_waypoints_; i++) {
    dist.push_back(dist[i - 1] + (path_waypoints_.block<1, 3>(i, 0) -
                                  path_waypoints_.block<1, 3>(i - 1, 0))
                                     .norm());
  }

  // for (int i = 0; i<dist.size(); i++) {
  //   i_switch_.push_back(num_pred_waypoints_*dist[i]/dist[-1]);
  // }
  for (auto it : dist) {
    i_switch_.push_back(num_pred_waypoints_ * it / dist.back());
  }
  // if (num_pred_waypoints_ <= i_switch_.back()) {
  //   num_pred_waypoints_ += 1;
  // }

  vel_guess_ = guessed_velocity;
  time_guess_ = dist.back() / vel_guess_;

  std::cout << "distance between each pair of nodes: " << dist << std::endl;
  for (auto it : i_switch_)
    std::cout << "switch index: " << it << std::endl;
  std::cout << "we will generate " << num_pred_waypoints_ << " waypoints."
            << std::endl;
  std::cout << "the guessed travel time is : " << time_guess_ << "[s]"
            << std::endl;

  init_state_guess(false);
}

void DMCCUAVManipulator::init_state_guess(
    const bool with_additional_relaxation) {

  // Eigen::MatrixXd u0(num_pred_waypoints_, num_dofs_ - 2);
  std::vector<double> u0(num_pred_waypoints_ * num_controls_,
                         (mass_quadrotor_ + mass_manipulator_) * g_ * 0.25);

  if (num_controls_ == 7) {
    for (int i = 0; i < num_pred_waypoints_; i++) {
      u0[4 * (i + 1)] = 0.0;
    }
  }

  int i_wp = 0;
  Eigen::Vector3d wp_last;
  Eigen::Vector3d wp_next;
  double interp;
  Eigen::MatrixXd x_guess_matrix =
      Eigen::MatrixXd::Zero(num_pred_waypoints_, num_states_);

  if (num_dofs_ == 6) {
    x_guess_matrix.block<1, 12>(0, 0) = init_pose_;
  } else if (num_dofs_ == 7) {
    x_guess_matrix.block<1, 14>(0, 0) = init_pose_;
  }

  std::vector<double> x0(init_pose_.data(), init_pose_.data() + num_dofs_);
  std::vector<int> lambda0;

  if (has_contact_target_) {
    //
  } else {
    for (int i = 0; i < num_pred_waypoints_; ++i) {
      if (i > i_switch_[i_wp]) {
        i_wp += 1;
      }
      if (i_wp == 0) {
        wp_last = init_pose_.block<3, 1>(0, 0).transpose();
      } else {
        wp_last = path_waypoints_.block<1, 3>(i_wp - 1, 0);
      }
      wp_next = path_waypoints_.block<1, 3>(i_wp, 0);
      if (i_wp > 0) {
        interp = ((double)i - i_switch_[i_wp - 1]) /
                 (i_switch_[i_wp] - i_switch_[i_wp - 1]);
      } else {
        interp = (double)i / (i_switch_[i_wp] + 1e-6);
      }
      Eigen::Vector3d pos_guess = (1.0 - interp) * wp_last + interp * wp_next;
      Eigen::Vector3d vel_guess =
          vel_guess_ * (wp_next - wp_last) / (wp_next - wp_last).norm();

      if (i > 0) {
        if (num_dofs_ == 6) {
          x_guess_matrix.block<1, 3>(i, 0) = pos_guess;
          x_guess_matrix.block<1, 3>(i, 6) = vel_guess;
          std::vector<double> pos_std(pos_guess.data(),
                                      pos_guess.data() +
                                          pos_guess.rows() * pos_guess.cols());
          // std::vector<double> vel_std(vel_guess.data(),
          //                             vel_guess.data() +
          //                                 vel_guess.rows() *
          //                                 vel_guess.cols());
          x0.insert(x0.end(), pos_std.begin(), pos_std.end());
          x0.push_back(0.0);
          x0.push_back(0.0);
          x0.push_back(0.0);
          // x0.insert(x0.end(), vel_std.begin(), vel_std.end());
          // x0.push_back(0.0);
          // x0.push_back(0.0);
          // x0.push_back(0.0);
        } else if (num_dofs_ == 7) {
          x_guess_matrix.block<1, 3>(i, 0) = pos_guess;
          x_guess_matrix(i, 6) = M_PI;
          x_guess_matrix.block<1, 3>(i, 7) = vel_guess;

          std::vector<double> pos_std(pos_guess.data(),
                                      pos_guess.data() +
                                          pos_guess.rows() * pos_guess.cols());
          // std::vector<double> vel_std(vel_guess.data(),
          //                             vel_guess.data() +
          //                                 vel_guess.rows() *
          //                                 vel_guess.cols());
          x0.insert(x0.end(), pos_std.begin(), pos_std.end());
          x0.push_back(0.0);
          x0.push_back(0.0);
          x0.push_back(0.0);
          x0.push_back(M_PI / 2.0);
          // x0.insert(x0.end(), vel_std.begin(), vel_std.end());
          // x0.push_back(0.0);
          // x0.push_back(0.0);
          // x0.push_back(0.0);
          // x0.push_back(0.0);
        } else {
          std::cout << "Error: undefined degrees of freedom!" << std::endl;
        }
      }

      std::vector<int> lambda_temp(num_def_waypoints_, 0.0);
      if (i != num_pred_waypoints_ - 1) {
        for (int j = 0; j < lambda_temp.size(); j++) {
          if (j >= i_wp)
            lambda_temp[j] = 1;
        }
      }

      lambda0.insert(lambda0.end(), lambda_temp.begin(), lambda_temp.end());
    }

    // fill in the guessed speed of each waypoint
    for (int i = 0; i < i_switch_.size(); i++) {
      path_waypoints_.block<1, 3>(i, num_dofs_) =
          x_guess_matrix.block<1, 3>(i, 6);
    }

    std::vector<double> tolerance_param0(
        (num_pred_waypoints_ - 1) * num_def_waypoints_, 0.0);

    // guessed optimization states
    init_values_.clear();
    init_values_.push_back(time_guess_);
    init_values_.insert(init_values_.end(), x0.begin(), x0.end());
    init_values_.insert(init_values_.end(), u0.begin(), u0.end());
    init_values_.insert(init_values_.end(), lambda0.begin(), lambda0.end());
    init_values_.insert(init_values_.end(), tolerance_param0.begin(),
                        tolerance_param0.end());

  } // no contact target
}

void DMCCUAVManipulator::initialization_formulation(
    const int kappa0, const int num_def_waypoints,
    const int num_pred_waypoints) {
  kappa0_ = kappa0;
  num_def_waypoints_ = num_def_waypoints;
  num_pred_waypoints_ = num_pred_waypoints;
  // CasADi formulation
  auto Tn = ca::MX::sym("Tn");
  dt_ = Tn / (num_pred_waypoints_ - 1);

  auto U = ca::MX::sym("U", num_controls_, num_pred_waypoints_);
  auto X = ca::MX::sym("X", num_dofs_, num_pred_waypoints_);
  auto init_pose = ca::MX::sym("init_pose", num_states_);
  auto waypoint_reference =
      ca::MX::sym("W_ref", num_states_, num_def_waypoints_);

  if (has_contact_target_) {

  } // has contact target
  else {
    lambda_param_ =
        ca::MX::sym("lambda_param", num_def_waypoints_, num_pred_waypoints_);
    tolerance_param_ = ca::MX::sym("tolerance_param", num_def_waypoints_,
                                   num_pred_waypoints_ - 1);
  }
  obj_function_ = Tn;
  ca::DM Q = ca::DM::eye(5) * 10.0;
  ca::DM u0 = ca::DM({1.0, 1.0, 1.0, 1.0, 0.0}) *
              (mass_quadrotor_ + mass_manipulator_) * g_;
  for (int i = 0; i < num_pred_waypoints_; i++) {
    obj_function_ += 0.0003 * dt_ *
                     ca::MX::mtimes({(U(all, i) - u0).T(), Q, U(all, i) - u0});
  }

  // ca::DM init_pose_dm{std::vector<double>(
  //     init_pose_.data(), init_pose_.size() + init_pose_.data())};
  // auto end_pose = path_waypoints_.block<1, 14>(num_def_waypoints_ - 1, 0);
  // ca::DM end_pose_dm{
  //     std::vector<double>(end_pose.data(), end_pose.size() +
  //     end_pose.data())};

  // constraint_vector_ = X(0, all).T() - init_pose_dm(ca::Slice(0, num_dofs_));
  // using parameters
  constraint_vector_ = X(all, 0) - init_pose(ca::Slice(0, num_dofs_));
  constraint_vector_ = ca::MX::reshape(constraint_vector_, -1, 1);

  constraint_vector_ = ca::MX::vertcat(
      {constraint_vector_,
       X(all, -1) - waypoint_reference(ca::Slice(0, num_dofs_), -1)});
  constraint_vector_ = ca::MX::vertcat({constraint_vector_, U(all, 0) - u0});
  constraint_vector_ = ca::MX::vertcat({constraint_vector_, U(all, -1) - u0});

  // discrete Lagrange equations
  auto q_nm1 = ca::MX::sym("q_nm1", num_dofs_);
  auto q_n = ca::MX::sym("q_n", num_dofs_);
  auto q_np1 = ca::MX::sym("q_np1", num_dofs_);

  auto D2L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_nm1, q_n, get_lagrangian_function()),
      q_n);

  auto D1L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_n, q_np1, get_lagrangian_function()),
      q_n);

  ca::Function derivative_EulerLagrange_function =
      ca::Function("dEL", {q_nm1, q_n, q_np1}, {D2L_d + D1L_d});

  auto q_b = ca::MX::sym("q_b", num_dofs_);
  auto q_b_dot = ca::MX::sym("q_b_dot", num_dofs_);
  std::vector<ca::MX> q_b_vec(2);
  q_b_vec[0] = q_b;
  q_b_vec[1] = q_b_dot;
  auto l_dot_function = get_derivative_lagrangian_function();
  auto D2L = l_dot_function(q_b_vec).at(0);
  auto d_EulerLagrange_init_function =
      ca::Function("dEl_init", {q_b, q_b_dot, q_n, q_np1}, {D2L + D1L_d});
  auto d_EulerLagrange_end_function =
      ca::Function("dEl_end", {q_b, q_b_dot, q_nm1, q_n}, {-D2L + D2L_d});

  // waypoint connection
  for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
    auto f_d_nm1 = discrete_forces(dt_, get_system_dynamics_function(),
                                   (ca::MX)X(all, i - 1), (ca::MX)U(all, i - 1),
                                   (ca::MX)U(all, i));
    auto f_d_n =
        discrete_forces(dt_, get_system_dynamics_function(), (ca::MX)X(all, i),
                        (ca::MX)U(all, i), (ca::MX)U(all, i + 1));
    std::vector<ca::MX> input(3);
    input[0] = X(all, i - 1);
    input[1] = X(all, i);
    input[2] = X(all, i + 1) + f_d_nm1 + f_d_n;
    auto sum = derivative_EulerLagrange_function(input).at(0);
    constraint_vector_ = ca::MX::vertcat({constraint_vector_, sum});
  }

  // boundary conditions
  auto f_0 =
      discrete_forces(dt_, get_system_dynamics_function(), (ca::MX)X(all, 0),
                      (ca::MX)U(all, 0), (ca::MX)U(all, 1));
  std::vector<ca::MX> input(4);
  input[0] = init_pose(ca::Slice(0, num_dofs_));
  input[1] = init_pose(ca::Slice(num_dofs_, 2 * num_dofs_));
  input[2] = X(all, 0);
  input[3] = X(all, 1) + f_0;
  auto init_condition = d_EulerLagrange_init_function(input).at(0);
  constraint_vector_ = ca::MX::vertcat({constraint_vector_, init_condition});

  auto f_N_1 = discrete_forces(dt_, get_system_dynamics_function(),
                               (ca::MX)X(all, num_pred_waypoints_ - 2),
                               (ca::MX)U(all, num_pred_waypoints_ - 2),
                               (ca::MX)U(all, num_pred_waypoints_ - 1));

  // input[0] = end_pose_dm(ca::Slice(0, num_dofs_));
  // input[1] = end_pose_dm(ca::Slice(num_dofs_, 2 * num_dofs_));
  input[0] = waypoint_reference(ca::Slice(0, num_dofs_), -1);
  input[1] = waypoint_reference(ca::Slice(num_dofs_, 2 * num_dofs_), -1);
  input[2] = X(all, num_pred_waypoints_ - 2);
  input[3] = X(all, num_pred_waypoints_ - 1) + f_N_1;

  auto end_condition = d_EulerLagrange_end_function(input).at(0);
  constraint_vector_ = ca::MX::vertcat({constraint_vector_, end_condition});

  if (has_contact_target_) {
  } else {
    for (int i = 0; i < num_def_waypoints_; ++i) {
      constraint_vector_ =
          ca::MX::vertcat({constraint_vector_, lambda_param_(i, 0) - 1.0});
      constraint_vector_ =
          ca::MX::vertcat({constraint_vector_, lambda_param_(i, -1)});
    }
    for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
      for (int j = 0; j < num_def_waypoints_; ++j) {
        auto mu = lambda_param_(j, i) - lambda_param_(j, i + 1);
        auto cost = ca::MX::norm_2(X(ca::Slice(0, 3), i) -
                                   waypoint_reference(ca::Slice(0, 3), j)) -
                    tolerance_param_(j, i);
        constraint_vector_ = ca::MX::vertcat({constraint_vector_, mu * cost});

      }
    }
  } // no contact target

  number_equal_constraint = constraint_vector_.size1();

  std::cout << "number of equal constraints: " << number_equal_constraint
            << std::endl;

  // velocity constraints
  for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
    constraint_vector_ = ca::MX::vertcat(
        {constraint_vector_,
         average_velocity(dt_, (ca::MX)X(all, i), (ca::MX)X(all, i + 1), 7)});

  }

  if (has_contact_target_) {
  } else {
    ca::MXDict nlp = {
        {"x", ca::MX::vertcat({Tn, ca::MX::reshape(X, -1, 1),
                               ca::MX::reshape(U, -1, 1),
                               ca::MX::reshape(lambda_param_, -1, 1),
                               ca::MX::reshape(tolerance_param_, -1, 1)})},
        {"f", obj_function_},
        {"g", constraint_vector_},
        {"p", ca::MX::vertcat({ca::MX::reshape(init_pose, -1, 1),
                               ca::MX::reshape(waypoint_reference, -1, 1)})}};
    std::string solver_name = "ipopt";
    casadi::Dict solver_opts;
    solver_opts["expand"] = true;
    solver_opts["ipopt.max_iter"] = 1000;
    solver_opts["ipopt.print_level"] = 3;
    solver_opts["print_time"] = 0;
    solver_opts["ipopt.acceptable_tol"] = 1e-8;
    solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
    solver_ = casadi::nlpsol("nlpsol", solver_name, nlp, solver_opts);
  }
}

void DMCCUAVManipulator::set_constraints(
    const std::vector<double> v, const std::vector<double> d_rpy,
    const std::vector<double> tolerance_param,
    const std::vector<double> high_interation) {
  for (int i = 0; i < number_equal_constraint; i++) {
    lbg_.push_back(0.0);
    ubg_.push_back(0.0);
  }
}

void DMCCUAVManipulator::set_constraints(
    const std::vector<double> v, const std::vector<double> d_rpy,
    const std::vector<double> min_state, const std::vector<double> max_state,
    const std::vector<double> min_control_input,
    const std::vector<double> max_control_input,
    const std::vector<double> lambda, std::vector<double> nu) {
  for (int i = 0; i < number_equal_constraint; i++) {
    lbg_.push_back(0.0);
    ubg_.push_back(0.0);
  }

  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    for (int j = 0; j < v.size(); j++) {
      lbg_.push_back(-v[j]);
      ubg_.push_back(v[j]);
    }
    for (int j = 0; j < d_rpy.size(); j++) {
      lbg_.push_back(-d_rpy[j]);
      ubg_.push_back(d_rpy[j]);
    }
  }

  // Tn
  lbx_.push_back(0.01);
  ubx_.push_back(ca::inf);
  // X
  for (int i = 0; i < num_pred_waypoints_; i++) {
    lbx_.insert(lbx_.end(), min_state.begin(), min_state.end());
    ubx_.insert(ubx_.end(), max_state.begin(), max_state.end());
  }
  // U
  for (int i = 0; i < num_pred_waypoints_; i++) {
    lbx_.insert(lbx_.end(), min_control_input.begin(), min_control_input.end());
    ubx_.insert(ubx_.end(), max_control_input.begin(), max_control_input.end());
  }
  // lambda
  for (int i = 0; i < num_pred_waypoints_; i++) {
    for (int j = 0; j < num_def_waypoints_; j++) {
      lbx_.push_back(lambda[0]);
      ubx_.push_back(lambda[1]);
    }
  }
  // nu
  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    for (int j = 0; j < num_def_waypoints_; j++) {
      lbx_.push_back(nu[0]);
      ubx_.push_back(nu[1]);
    }
  }
}

void DMCCUAVManipulator::get_results() {
  optimization_param_.clear();
  std::vector<double> init_pose_std(init_pose_.data(),
                                    init_pose_.size() + init_pose_.data());
  std::vector<double> reference_waypoint_std(
      path_waypoints_.data(), path_waypoints_.size() + path_waypoints_.data());
  std::merge(init_pose_std.begin(), init_pose_std.end(),
             reference_waypoint_std.begin(), reference_waypoint_std.end(),
             std::back_inserter(optimization_param_));

  std::cout << init_values_ << std::endl;
  ca::DMDict arg = {{"lbx", lbx_},        {"ubx", ubx_},
                    {"lbg", lbg_},        {"ubg", ubg_},
                    {"x0", init_values_}, {"p", optimization_param_}};

  auto start_time = std::chrono::high_resolution_clock::now();

  opt_results_ = solver_(arg);

  std::vector<double> result_all(opt_results_.at("x"));
  std::vector<double> result_x, result_u;
  double opt_time = result_all[0];
  result_x.assign(result_all.begin() + 1,
                  result_all.begin() + num_pred_waypoints_ * num_dofs_ + 1);
  Eigen::MatrixXd result_x_matrix =
      Eigen::MatrixXd::Map(result_x.data(), num_dofs_, num_pred_waypoints_);

  result_u.assign(
      result_all.begin() + num_pred_waypoints_ * num_dofs_ + 1,
      result_all.begin() +
          num_pred_waypoints_ * (num_controls_ + num_dofs_) + 1);
  Eigen::MatrixXd result_u_matrix =
      Eigen::MatrixXd::Map(result_u.data(), num_controls_, num_pred_waypoints_);

  Eigen::MatrixXd result_lambda_matrix;

  if (has_contact_target_)
  {}
  else{
    std::vector<double> result_lambda, result_nu;
    result_lambda.assign(
        result_all.begin() + num_pred_waypoints_ * (num_controls_ + num_dofs_) +
            1,
        result_all.begin() + num_pred_waypoints_ * (num_controls_ + num_dofs_) +
            1 + num_pred_waypoints_ * num_def_waypoints_);
    result_lambda_matrix = Eigen::MatrixXd::Map(
        result_lambda.data(), num_def_waypoints_, num_pred_waypoints_);
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
      stop_time - start_time);
  std::cout << "estimated travel time: " << opt_time << std::endl;
  std::cout << "average calculation time for each iteration [s]: "
            << duration.count() / 1e6 << std::endl;
            // std::cout << result_x << std::endl;
            // std::cout << result_u<< std::endl;
  std::cout << result_x_matrix.transpose() << std::endl;
  std::cout << "======" << std::endl;
  std::cout << result_u_matrix.transpose() << std::endl;

  if (has_contact_target_){

  }
  else{
    std::cout << result_lambda_matrix << std::endl;
  }

  std::cout << "full results: " << result_all << std::endl;
}

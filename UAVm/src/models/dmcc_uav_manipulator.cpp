/*
 * @Author: Wei Luo
 * @Date: 2022-12-14 11:00:34
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-19 14:03:33
 * @Note: Note
 */

#include <dmcc_uav_manipulator.hpp>

DMCCUAVManipulator::DMCCUAVManipulator(
    const double mass_quadrotor, const double mass_manipulator,
    const std::vector<double> inertia_moment,
    const std::vector<double> manipulator_inertia_moment,
    const double manipulator_length, const bool has_contact_target,
    const std::vector<double> montage_offset_b, const double frame_size,
    const double motor_torque_const, const double g)
    : DerivationUAVm(g, manipulator_length, mass_manipulator, mass_quadrotor,
                     inertia_moment, manipulator_inertia_moment, frame_size,
                     motor_torque_const, montage_offset_b) {

  // only a trajectory planning or considering a specific mission
  has_contact_target_ = has_contact_target;
};
DMCCUAVManipulator::~DMCCUAVManipulator(){};

void DMCCUAVManipulator::initialization_formulation(
    const int num_def_waypoints, const int num_pred_waypoints) {

  num_def_waypoints_ = num_def_waypoints;
  num_pred_waypoints_ = num_pred_waypoints;

  std::cout << "number of waypoint prediction :" << num_pred_waypoints_
            << std::endl;

  // CasADi formulation
  auto Tn = ca::MX::sym("Tn");
  dt_ = Tn / (num_pred_waypoints_ - 1);
  auto U = ca::MX::sym("U", num_controls_, num_pred_waypoints_);
  auto X = ca::MX::sym("X", num_dofs_, num_pred_waypoints_);
  auto init_pose = ca::MX::sym("init_pose", num_states_);
  auto waypoint_reference =
      ca::MX::sym("W_ref", num_states_, num_def_waypoints_);

  if (has_contact_target_) {
    epsilon_param_ = ca::MX::sym("epsilon_param", num_pred_waypoints_ - 1);
    kappa_param_ = ca::MX::sym("kappa_param", num_pred_waypoints_);
    contact_relax_param_ =
        ca::MX::sym("contact_relax_param", num_pred_waypoints_ - 1);
  } // has contact target
  else {
    // lambda_param_ =
    // ca::MX::sym("lambda_param", num_def_waypoints_, num_pred_waypoints_);
    lambda_param_ =
        ca::MX::sym("lambda_param", num_pred_waypoints_, num_def_waypoints_);
    tolerance_param_ = ca::MX::sym("tolerance_param", num_def_waypoints_,
                                   num_pred_waypoints_ - 1);
  }

  obj_function_ = Tn;
  ca::DM u0 = ca::DM({1.0, 1.0, 1.0, 1.0, 0.0}) *
              (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ / 4.0;
  std::cout << mass_quadrotor_ << std::endl;
  std::cout << mass_manipulator_ << std::endl;
  std::cout << (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ / 4.0
            << std::endl;
  ca::DM Q = ca::DM::eye(5) * 10.0;

  for (int i = 0; i < num_pred_waypoints_; i++) {
    obj_function_ += 0.0003 * dt_ *
                     ca::MX::mtimes(ca::MX::mtimes((U(all, i) - u0).T(), Q), U(all, i) - u0
  );
    // obj_function_ +=
    //     0.0003 * dt_ * ca::MX::dot(Q * (U(all, i) - u0).T(), (U(all, i) -
    //     u0).T());
  }

  // using parameters
  constraint_vector_std_.push_back(X(all, 0) -
                                   init_pose(ca::Slice(0, num_dofs_)));
  constraint_vector_std_.push_back(
      X(all, num_pred_waypoints_ - 1) -
      waypoint_reference(ca::Slice(0, num_dofs_), num_def_waypoints_ - 1));

  constraint_vector_std_.push_back(U(all, 0) - u0);
  constraint_vector_std_.push_back(U(all, -1) - u0);

  // discrete Lagrange equations
  ca::MX q_nm1 = ca::MX::sym("q_nm1", num_dofs_);
  ca::MX q_n = ca::MX::sym("q_n", num_dofs_);
  ca::MX q_np1 = ca::MX::sym("q_np1", num_dofs_);

  ca::MX D2L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_nm1, q_n, get_lagrangian_function()),
      q_n);

  ca::MX D1L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_n, q_np1, get_lagrangian_function()),
      q_n);

  ca::Function derivative_EulerLagrange_function =
      ca::Function("dEL", {q_nm1, q_n, q_np1}, {D2L_d + D1L_d});

  ca::MX q_b = ca::MX::sym("q_b", num_dofs_);
  ca::MX q_b_dot = ca::MX::sym("q_b_dot", num_dofs_);
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
  for (int i = 1; i < num_pred_waypoints_ - 1; ++i) {
    auto f_d_nm1 = discrete_forces(dt_, get_external_force_function(),
                                   (ca::MX)X(all, i - 1), (ca::MX)U(all, i - 1),
                                   (ca::MX)U(all, i));
    auto f_d_n =
        discrete_forces(dt_, get_external_force_function(), (ca::MX)X(all, i),
                        (ca::MX)U(all, i), (ca::MX)U(all, i + 1));
    std::vector<ca::MX> input(3);
    input[0] = X(all, i - 1);
    input[1] = X(all, i);
    input[2] = X(all, i + 1) ;
    auto sum = derivative_EulerLagrange_function(input).at(0) + f_d_nm1 + f_d_n;
    constraint_vector_std_.push_back(sum);
  }

  // boundary conditions
  auto f_0 =
      discrete_forces(dt_, get_external_force_function(), (ca::MX)X(all, 0),
                      (ca::MX)U(all, 0), (ca::MX)U(all, 1));

  std::cout << " f0 size is " << f_0.size() << std::endl;
  std::vector<ca::MX> input(4);
  input[0] = init_pose(ca::Slice(0, num_dofs_));
  input[1] = init_pose(ca::Slice(num_dofs_, 2 * num_dofs_));
  input[2] = X(all, 0);
  input[3] = X(all, 1);
  auto init_condition = d_EulerLagrange_init_function(input).at(0) + f_0;
  constraint_vector_std_.push_back(init_condition);

  auto f_N_1 = discrete_forces(dt_, get_external_force_function(),
                               (ca::MX)X(all, num_pred_waypoints_ - 2),
                               (ca::MX)U(all, num_pred_waypoints_ - 2),
                               (ca::MX)U(all, num_pred_waypoints_ - 1));

  // input[0] = end_pose_dm(ca::Slice(0, num_dofs_));
  // input[1] = end_pose_dm(ca::Slice(num_dofs_, 2 * num_dofs_));
  input[0] = waypoint_reference(ca::Slice(0, num_dofs_), -1);
  input[1] = waypoint_reference(ca::Slice(num_dofs_, 2 * num_dofs_), -1);
  input[2] = X(all, num_pred_waypoints_ - 2);
  input[3] = X(all, num_pred_waypoints_ - 1);

  auto end_condition = d_EulerLagrange_end_function(input).at(0) + f_N_1;
  constraint_vector_std_.push_back(end_condition);

  if (has_contact_target_) {
    // # initial value of kappa_param = init kappa
    // endpoint value of kappa_param = 0
    constraint_vector_std_.push_back(kappa_param_(0) - kappa0_);
    constraint_vector_std_.push_back(kappa_param_(-1));
    // kappa_param[i+1] = kappa_param[i]-epsilon[i]
    for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
      constraint_vector_std_.push_back(kappa_param_(i + 1) - kappa_param_(i) +
                                       epsilon_param_(i));
    }
    // epsilon_param[i]*(|End_Manip-pickup_point|-relax)==0
    ca::Function end_effector_pos_function =
        get_end_effector_position_function();
    for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
      ca::MX dist_temp =
          ca::MX::norm_2(end_effector_pos_function(X(all, i)).at(0) -
                         target_position_function_({dt_, (double)i}).at(0));
      constraint_vector_std_.push_back(epsilon_param_(i) *
                                       (dist_temp - contact_relax_param_(i)));
    }
  } else {
    for (int i = 0; i < num_def_waypoints_; ++i) {
      constraint_vector_std_.push_back(lambda_param_(0, i) - 1.0);
      constraint_vector_std_.push_back(lambda_param_(-1, i));
    }
    for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
      for (int j = 0; j < num_def_waypoints_; ++j) {
        // auto mu = lambda_param_(j, i) - lambda_param_(j, i + 1);
        ca::MX mu = lambda_param_(i, j) - lambda_param_(i + 1, j);
        ca::MX cost = ca::MX::norm_2(X(ca::Slice(0, 3), i) -
                                     waypoint_reference(ca::Slice(0, 3), j)) -
                      tolerance_param_(j, i); // ca::Slice(0, 3)
        constraint_vector_std_.push_back(mu * cost);
      }
    }
  } // no contact target

  number_equality_constraint_ =
      ca::MX::vertcat(constraint_vector_std_).size().first;

  std::cout << "number of equal constraints: " << number_equality_constraint_
            << std::endl;

  // velocity constraints
  for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
    constraint_vector_std_.push_back(
        average_velocity(dt_, (ca::MX)X(all, i), (ca::MX)X(all, i + 1), 7));
  }

  // auto a =
  //     ca::MX::vertcat({Tn, ca::MX::reshape(X, -1, 1), ca::MX::reshape(U, -1,
  //     1),
  //                      ca::MX::reshape(lambda_param_, -1, 1),
  //                      ca::MX::reshape(tolerance_param_, -1, 1)});
  // std::cout << a.size() << std::endl;

  if (has_contact_target_) {
    //
    ca::Function end_effector_pos_function =
        get_end_effector_position_function();
    ca::Function end_effector_vel_function =
        get_end_effector_velocity_function();
    std::vector<ca::MX> input(2);
    for (int i = 0; i < num_pred_waypoints_ - 1; ++i) {
      input.at(0) = (ca::MX)X(all, i);
      input.at(1) =
          average_velocity(dt_, (ca::MX)X(all, i), (ca::MX)X(all, i + 1), 7);
      ca::MX v1 = end_effector_vel_function(input).at(0);
      ca::MX v2 = target_velocity_function_({dt_, (double)i}).at(0);
      constraint_vector_std_.push_back(epsilon_param_(i) *
                                       ca::MX::norm_2(v1 - v2));

      ca::DM unit_e = ca::DM::zeros(3);
      unit_e(0) = 1.0;
      ca::MX angle = X(ca::Slice(3, 6), i);
      ca::MX unit_heading_uav =
          ca::MX::mtimes({rotation_matrix(angle), unit_e});
      ca::MX temp_result =
          unit_heading_uav(0) * v2(1) - v2(0) * unit_heading_uav(1);
      constraint_vector_std_.push_back(epsilon_param_(i) * temp_result *
                                       temp_result);
    }

    for (int i = 0; i < num_pred_waypoints_; ++i) {
      ca::MX end_effector_pos = end_effector_pos_function(X(all, i)).at(0);
      ca::MX target_position =
          target_position_function_({dt_, (double)i}).at(0);
      constraint_vector_std_.push_back(end_effector_pos(2) -
                                       target_position(2));
    }

    ca::MXDict nlp = {
        {"x", ca::MX::vertcat({Tn, ca::MX::reshape(X, -1, 1),
                               ca::MX::reshape(U, -1, 1),
                               ca::MX::reshape(epsilon_param_, -1, 1),
                               ca::MX::reshape(kappa_param_, -1, 1),
                               ca::MX::reshape(contact_relax_param_, -1, 1)})},
        {"f", obj_function_},
        {"g", ca::MX::vertcat(constraint_vector_std_)},
        {"p", ca::MX::vertcat({ca::MX::reshape(init_pose, -1, 1),
                               ca::MX::reshape(waypoint_reference, -1, 1)})}};
    std::string solver_name = "ipopt";
    casadi::Dict solver_opts;
    solver_opts["expand"] = true;
    solver_opts["ipopt.max_iter"] = 1000;
    solver_opts["ipopt.print_level"] = 3;
    solver_opts["print_time"] = 0;
    solver_opts["ipopt.acceptable_tol"] = 1e-5;
    solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
    solver_ = casadi::nlpsol("nlpsol", solver_name, nlp, solver_opts);

  } else {
    for (int i = 0; i < num_pred_waypoints_; ++i) {
      for (int j = 0; j < num_def_waypoints_ - 1; ++j) {
        constraint_vector_std_.push_back(lambda_param_(i, j) -
                                         lambda_param_(i, j + 1));
      }
    }
    ca::MXDict nlp = {
        {"x", ca::MX::vertcat({Tn, ca::MX::reshape(X, -1, 1),
                               ca::MX::reshape(U, -1, 1),
                               ca::MX::reshape(lambda_param_, -1, 1),
                               ca::MX::reshape(tolerance_param_, -1, 1)})},
        {"f", obj_function_},
        {"g", ca::MX::vertcat(constraint_vector_std_)},
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

/* * get path waypoints without contact target * */
void DMCCUAVManipulator::get_path_waypoints(
    const Eigen::VectorXd current_pose, const Eigen::MatrixXd path_waypoints,
    const double guessed_velocity) {
  init_pose_ = current_pose;
  path_waypoints_ = path_waypoints;
  std::cout << "path waypoints_ are " << path_waypoints_ << std::endl;
  if (num_def_waypoints_ == path_waypoints.cols())
    std::cout << "num of defined path waypoints: " << num_def_waypoints_
              << std::endl;
  else {
    std::cout << "num of defined path waypoints is not equal to that is "
                 "initalized !!!"
              << std::endl;
  }

  std::vector<double> dist;
  dist.push_back(
      (init_pose_.block<3, 1>(0, 0) - path_waypoints_.block<3, 1>(0, 0))
          .norm());

  if (dist[0] < 0.01) {
    dist[0] += 0.01; // incase of 0 or almost zero distance
  }

  for (int i = 1; i < num_def_waypoints_; i++) {
    dist.push_back(dist[i - 1] + (path_waypoints_.block<3, 1>(0, i) -
                                  path_waypoints_.block<3, 1>(0, i - 1))
                                     .norm());
  }

  for (auto it : dist) {
    i_switch_.push_back(num_pred_waypoints_ * it / dist.back());
  }

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

/* * get path waypoints without contact target * */
void DMCCUAVManipulator::get_path_waypoints(
    const int kappa0, const Eigen::VectorXd current_pose,
    const Eigen::MatrixXd path_waypoints, const double guessed_velocity) {
  kappa0_ = kappa0;
  init_pose_ = current_pose;
  path_waypoints_ = path_waypoints;
  std::cout << "path waypoints_ are " << path_waypoints_ << std::endl;
  if (num_def_waypoints_ == path_waypoints.cols())
    std::cout << "num of defined path waypoints: " << num_def_waypoints_
              << std::endl;
  else {
    std::cout << "num of defined path waypoints is not equal to that is "
                 "initalized !!!"
              << std::endl;
  }

  std::vector<double> dist;
  dist.push_back(
      (init_pose_.block<3, 1>(0, 0) - path_waypoints_.block<3, 1>(0, 0))
          .norm());

  if (dist.at(0) < 0.01) {
    dist.at(0) += 0.01; // incase of 0 or almost zero distance
  }

  for (int i = 1; i < num_def_waypoints_; i++) {
    dist.push_back(dist[i - 1] + (path_waypoints_.block<3, 1>(0, i) -
                                  path_waypoints_.block<3, 1>(0, i - 1))
                                     .norm());
  }

  for (auto it : dist) {
    i_switch_.push_back(num_pred_waypoints_ * it / dist.back());
  }

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

  Eigen::MatrixXd u0_eigen =
      Eigen::MatrixXd::Zero(num_controls_, num_pred_waypoints_);
  for (int i = 0; i < num_pred_waypoints_; i++) {
    u0_eigen.col(i) << (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ *
                           0.25,
        (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ * 0.25,
        (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ * 0.25,
        (mass_quadrotor_ + mass_manipulator_) * g_acceleration_ * 0.25, 0.0;
  }

  std::cout << u0_eigen << std::endl;

  int i_wp = 0;
  Eigen::Vector3d wp_last;
  Eigen::Vector3d wp_next;
  double interp;
  Eigen::MatrixXd x_guess_matrix =
      Eigen::MatrixXd::Zero(num_states_, num_pred_waypoints_);

  if (num_dofs_ == 6) {
    x_guess_matrix.block(0, 0, 12, 1) = init_pose_;
  } else if (num_dofs_ == 7) {
    x_guess_matrix.block(0, 0, 14, 1) = init_pose_;
    x_guess_matrix.row(num_dofs_ - 1) =
        Eigen::MatrixXd::Constant(1, num_pred_waypoints_, M_PI / 2.0);
  }

  // without contact target
  // std::vector<int> lambda0_std;
  // std::vector<double> lambda0_std;
  Eigen::MatrixXd lambda0_eigen =
      Eigen::MatrixXd::Ones(num_pred_waypoints_, num_def_waypoints_);
  // with contact target
  std::vector<int> kappa0_std;
  std::vector<int> epsilon0_std(num_pred_waypoints_ - 1, 0);
  if (has_contact_target_) {
    for (int i = 0; i < num_pred_waypoints_; ++i) {
      if (i > i_switch_.at(i_wp)) {
        i_wp += 1;
      }
      if (i_wp == 0) {
        wp_last = init_pose_.block<3, 1>(0, 0);
      } else {
        wp_last = path_waypoints_.block(0, i_wp - 1, 3, 1);
      }
      wp_next = path_waypoints_.block(0, i_wp, 3, 1);
      if (i_wp > 0) {
        interp = (double)(i - i_switch_.at(i_wp - 1)) /
                 (double)(i_switch_.at(i_wp) - i_switch_.at(i_wp - 1) + 1e-6);
      } else {
        interp = (double)i / (i_switch_.at(i_wp) + 1e-6);
      }
      Eigen::Vector3d pos_guess = (1.0 - interp) * wp_last + interp * wp_next;
      Eigen::Vector3d vel_guess =
          vel_guess_ * (wp_next - wp_last) / (wp_next - wp_last).norm();

      if (i > 0) {
        x_guess_matrix.block(0, i, 3, 1) = pos_guess;
        x_guess_matrix.block(num_dofs_, i, 3, 1) = vel_guess;
        std::vector<double> pos_std(pos_guess.data(),
                                    pos_guess.data() +
                                        pos_guess.rows() * pos_guess.cols());
      }
    }

    for (int i = 0; i < num_pred_waypoints_; i++) {
      kappa0_std.push_back(kappa0_);
    }
    for (int i = 0; i < num_pred_waypoints_; i++) {
      if (i == i_switch_.at(0)) {
        kappa0_std.at(i) = kappa0_ - 1;
        epsilon0_std.at(i) = 1;
        int temp_index = i;
        while (kappa0_std.at(temp_index) != 0) {
          kappa0_std.at(temp_index + 1) = kappa0_std.at(temp_index) - 1;
          temp_index += 1;
          epsilon0_std.at(temp_index) = 1;
        }
        for (int j = i + kappa0_ - 1; j < num_pred_waypoints_; j++) {
          kappa0_std.at(j) = 0;
        }
      }
    }

    // guessed optimization states
    init_values_.clear();
    init_values_.push_back(time_guess_);
    auto init_x0 = x_guess_matrix.block(0, 0, num_dofs_, num_pred_waypoints_);
    std::vector<double> l0_std(num_pred_waypoints_ - 1, 0.0);
    std::vector<double> init_x0_std(init_x0.data(),
                                    init_x0.size() + init_x0.data());
    init_values_.insert(init_values_.end(), init_x0_std.begin(),
                        init_x0_std.end());
    std::vector<double> init_u0_std(u0_eigen.data(),
                                    u0_eigen.size() + u0_eigen.data());
    init_values_.insert(init_values_.end(), init_u0_std.begin(),
                        init_u0_std.end());
    init_values_.insert(init_values_.end(), epsilon0_std.begin(),
                        epsilon0_std.end());
    init_values_.insert(init_values_.end(), kappa0_std.begin(),
                        kappa0_std.end());
    init_values_.insert(init_values_.end(), l0_std.begin(), l0_std.end());

    if (with_additional_relaxation) {
      std::vector<double> v0_std(num_pred_waypoints_ - 1, 0.0);
      std::vector<double> gamma0_std(num_pred_waypoints_ - 1, 0.0);
      init_values_.insert(init_values_.end(), v0_std.begin(), v0_std.end());
      init_values_.insert(init_values_.end(), gamma0_std.begin(),
                          gamma0_std.end());
    }
  } else {
    for (int i = 0; i < num_pred_waypoints_; ++i) {
      if (i > i_switch_.at(i_wp)) {
        i_wp += 1;
      }
      if (i_wp == 0) {
        wp_last = init_pose_.block<3, 1>(0, 0);
      } else {
        wp_last = path_waypoints_.block(0, i_wp - 1, 3, 1);
      }
      wp_next = path_waypoints_.block(0, i_wp, 3, 1);
      if (i_wp > 0) {
        interp = (double)(i - i_switch_.at(i_wp - 1)) /
                 (double)(i_switch_.at(i_wp) - i_switch_.at(i_wp - 1) + 1e-6);
      } else {
        interp = (double)i / (i_switch_.at(i_wp) + 1e-6);
      }
      Eigen::Vector3d pos_guess = (1.0 - interp) * wp_last + interp * wp_next;
      Eigen::Vector3d vel_guess =
          vel_guess_ * (wp_next - wp_last) / (wp_next - wp_last).norm();

      if (i > 0) {
        x_guess_matrix.block(0, i, 3, 1) = pos_guess;
        x_guess_matrix.block(num_dofs_, i, 3, 1) = vel_guess;
        std::vector<double> pos_std(pos_guess.data(),
                                    pos_guess.data() +
                                        pos_guess.rows() * pos_guess.cols());
      }
      if (i != num_pred_waypoints_ - 1) {
        // std::vector<int> lambda_temp(num_def_waypoints_, 1);
        lambda0_eigen.row(i) =
            Eigen::MatrixXd::Constant(1, num_def_waypoints_, 1.0);
        for (int j = 0; j < i_wp; j++) {
          // lambda_temp[j] = 0;
          lambda0_eigen(i, j) = 0.0;
        }
        // std::cout << lambda_temp << std::endl;
        // lambda0_std.insert(lambda0_std.end(), lambda_temp.begin(),
        //                    lambda_temp.end());
      } else {
        lambda0_eigen.row(i) =
            Eigen::MatrixXd::Constant(1, num_def_waypoints_, 0.0);
        // std::vector<int> lambda_temp(num_def_waypoints_, 0);
        // std::cout << lambda_temp << std::endl;
        // lambda0_std.insert(lambda0_std.end(), lambda_temp.begin(),
        //                    lambda_temp.end());
      }
    }
    std::cout << lambda0_eigen << std::endl;
    std::vector<double> lambda0_std(
        lambda0_eigen.data(), lambda0_eigen.size() + lambda0_eigen.data());
    // fill in the guessed speed of each waypoint
    for (int i = 0; i < i_switch_.size(); i++) {
      path_waypoints_.block<3, 1>(num_dofs_, i) =
          x_guess_matrix.block<3, 1>(6, i);
    }

    std::vector<double> tolerance_param0(
        (num_pred_waypoints_ - 1) * num_def_waypoints_, 0.0);

    // guessed optimization states
    init_values_.clear();
    init_values_.push_back(time_guess_);
    auto init_x0 = x_guess_matrix.block(0, 0, num_dofs_, num_pred_waypoints_);
    std::vector<double> init_x0_std(init_x0.data(),
                                    init_x0.size() + init_x0.data());
    init_values_.insert(init_values_.end(), init_x0_std.begin(),
                        init_x0_std.end());
    std::vector<double> init_u0_std(u0_eigen.data(),
                                    u0_eigen.size() + u0_eigen.data());
    init_values_.insert(init_values_.end(), init_u0_std.begin(),
                        init_u0_std.end());

    init_values_.insert(init_values_.end(), lambda0_std.begin(),
                        lambda0_std.end());
    init_values_.insert(init_values_.end(), tolerance_param0.begin(),
                        tolerance_param0.end());
  } // no contact target

  std::cout << "Initial values: " << init_values_ << std::endl;
  std::cout << x_guess_matrix << std::endl;
}

void DMCCUAVManipulator::set_constraints(
    const std::vector<double> min_state, const std::vector<double> max_state,
    const std::vector<double> min_control_input,
    const std::vector<double> max_control_input,
    const double acceptable_distance, const std::vector<double> v,
    const std::vector<double> d_rpy,
    const double acceptable_velocity_difference,
    const double acceptable_heading_difference, const double high_interation) {
  lbg_.clear();
  ubg_.clear();
  lbx_.clear();
  ubx_.clear();
  // equality constraints
  for (int i = 0; i < number_equality_constraint_; i++) {
    lbg_.push_back(0.0);
    ubg_.push_back(0.0);
  }

  assert(v != 3);
  assert(d_rpy != 4);
  // velocity constraints
  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    for (int j = 0; j < v.size(); j++) {
      lbg_.push_back(-v.at(j));
      ubg_.push_back(v.at(j));
    }
    for (int j = 0; j < d_rpy.size(); j++) {
      lbg_.push_back(-d_rpy.at(j));
      ubg_.push_back(d_rpy.at(j));
    }
  }

  // heading and velocity constraints
  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    lbg_.push_back(0.0);
    lbg_.push_back(0.0);
    ubg_.push_back(acceptable_velocity_difference);
    ubg_.push_back(acceptable_heading_difference);
  }

  // altitude constraints
  for (int i = 0; i < num_pred_waypoints_; i++) {
    lbg_.push_back(high_interation);
    ubg_.push_back(ca::inf);
  }

  // states
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
  // epsilon
  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    lbx_.push_back(0.0);
    ubx_.push_back(1.0);
  }
  // kappa
  for (int i = 0; i < num_pred_waypoints_ - 1; i++) {
    lbx_.push_back(0.0);
    ubx_.push_back(kappa0_);
  }
  // contact distance
  for (int i = 0; i < num_pred_waypoints_; i++) {
    lbx_.push_back(0.0);
    ubx_.push_back(acceptable_distance);
  }
}

void DMCCUAVManipulator::set_constraints(
    const std::vector<double> v, const std::vector<double> d_rpy,
    const std::vector<double> min_state, const std::vector<double> max_state,
    const std::vector<double> min_control_input,
    const std::vector<double> max_control_input,
    const std::vector<double> lambda, std::vector<double> nu) {
  for (int i = 0; i < number_equality_constraint_; i++) {
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

  for (int i = 0; i < num_pred_waypoints_; i++) {
    for (int j = 0; j < num_def_waypoints_ - 1; j++) {
      lbg_.push_back(-ca::inf);
      ubg_.push_back(0.0);
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
  optimization_param_.insert(optimization_param_.end(), init_pose_std.begin(),
                             init_pose_std.end());
  optimization_param_.insert(optimization_param_.end(),
                             reference_waypoint_std.begin(),
                             reference_waypoint_std.end());
  std::cout << init_pose_std << std::endl;
  std::cout << reference_waypoint_std << std::endl;
  std::cout << optimization_param_ << std::endl;
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

  result_u.assign(result_all.begin() + num_pred_waypoints_ * num_dofs_ + 1,
                  result_all.begin() +
                      num_pred_waypoints_ * (num_controls_ + num_dofs_) + 1);
  Eigen::MatrixXd result_u_matrix =
      Eigen::MatrixXd::Map(result_u.data(), num_controls_, num_pred_waypoints_);

  Eigen::MatrixXd result_lambda_matrix;

  if (has_contact_target_) {
  } else {
    std::vector<double> result_lambda, result_nu;
    result_lambda.assign(
        result_all.begin() + num_pred_waypoints_ * (num_controls_ + num_dofs_) +
            1,
        result_all.begin() + num_pred_waypoints_ * (num_controls_ + num_dofs_) +
            1 + num_pred_waypoints_ * num_def_waypoints_);
    result_lambda_matrix = Eigen::MatrixXd::Map(
        result_lambda.data(), num_pred_waypoints_, num_def_waypoints_);
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

  if (has_contact_target_) {

  } else {
    std::cout << result_lambda_matrix << std::endl;
  }

  std::cout << "full results: " << result_all << std::endl;
}

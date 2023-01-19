/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:34:11
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-19 16:24:21
 * @Note: Note
 */

#include "dmoc_quadrotor.hpp"

DMOCUAV::DMOCUAV(const double dt, const double N, const double mass_quadrotor,
                 const std::vector<double> inertia_moment,
                 const double frame_size, const double motor_torque_const,
                 const double g)
    : DerivationUAV(g, mass_quadrotor, inertia_moment, frame_size,
                    motor_torque_const) {
  prediction_horizon_ = N;
  dt_ = dt;
}

void DMOCUAV::initialization_formulation() {
  ca::Function lagrangian_function = get_lagrangian_function();
  ca::Function d_lagrangian_function = get_derivative_lagrangian_function();
  ca::Function external_force_function = get_external_force_function();
  // MPC
  ca::MX U = ca::MX::sym("U", num_controls_, prediction_horizon_);
  ca::MX X = ca::MX::sym("X", num_dofs_, prediction_horizon_);
  // ca::MX current_state_ref = ca::MX::sym("current_state_ref", num_states_);
  ca::MX X_ref = ca::MX::sym("X_ref", num_states_, prediction_horizon_);

  ca::DM P_m = ca::DM::zeros(6, 6);
  P_m(0, 0) = 86.21;
  P_m(1, 1) = 86.21;
  P_m(2, 2) = 120.95;
  P_m(3, 3) = 6.94;
  P_m(4, 4) = 6.94;
  P_m(5, 5) = 11.04;
  P_m(0, 3) = 6.45;
  P_m(3, 0) = 6.45;
  P_m(1, 4) = 6.45;
  P_m(4, 1) = 6.45;
  P_m(2, 5) = 10.95;
  P_m(5, 2) = 10.95;

  ca::DM R_m = ca::DM::zeros(4, 4);
  R_m(0, 0) = 10.0;
  R_m(1, 1) = 10.0;
  R_m(2, 2) = 10.0;
  R_m(3, 3) = 10.0;

  ca::DM ref_u = ca::DM({1.0, 1.0, 1.0, 1.0}) *
                 mass_quadrotor_  * g_acceleration_ / 4.0;

  ca::MX obj = 0.0;

  // control cost
  for (int i = 0; i < prediction_horizon_; i++) {
    ca::MX temp_ = U(slice_all, i) - ref_u;
    obj += ca::MX::mtimes({temp_.T(), R_m, temp_});
  }
  // state cost
  for (int i = 0; i < prediction_horizon_; i++) {
    ca::MX temp_ = X(slice_state_, i) - X_ref(slice_state_, i);
    obj +=
        ca::MX::mtimes({temp_.T(), P_m, temp_});
  }

  // discrete lagrangian
  ca::MX q_nm1 = ca::MX::sym("q_nm1", num_dofs_);
  ca::MX q_n = ca::MX::sym("q_n", num_dofs_);
  ca::MX q_np1 = ca::MX::sym("q_np1", num_dofs_);
  ca::MX D2L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_nm1, q_n, lagrangian_function), q_n);
  ca::MX D1L_d = ca::MX::gradient(
      discrete_lagrange_verlet(dt_, q_n, q_np1, lagrangian_function), q_n);
  ca::Function d_EulerLagrange_function =
      ca::Function("dEL", {q_nm1, q_n, q_np1}, {D2L_d + D1L_d});

  ca::MX q_b = ca::MX::sym("q_b", num_dofs_);
  ca::MX q_b_dot = ca::MX::sym("q_b_dot", num_dofs_);
  std::vector<ca::MX> q_b_vec(2);
  q_b_vec[0] = q_b;
  q_b_vec[1] = q_b_dot;
  ca::MX D2L = d_lagrangian_function(q_b_vec).at(0);
  ca::Function d_EulerLagrange_init_function =
      ca::Function("dEl_init", {q_b, q_b_dot, q_n, q_np1}, {D2L + D1L_d});

  // constraints
  std::vector<ca::MX> constraint_vector;

  // constraint_vector.push_back(X(slice_state_, 0) -
  //                                 ca::MX::reshape(current_state_ref(slice_state_), -1, 1));
  constraint_vector.push_back(X(slice_state_, 0) - X_ref(slice_state_, 0));

  for (int i = 1; i < prediction_horizon_ - 1; i++) {
    ca::MX f_d_nm1 = discrete_forces(
        dt_, external_force_function, (ca::MX)X(slice_all, i - 1),
        (ca::MX)U(slice_all, i - 1), (ca::MX)U(slice_all, i));
    ca::MX f_d_n =
        discrete_forces(dt_, external_force_function, (ca::MX)X(slice_all, i),
                        (ca::MX)U(slice_all, i), (ca::MX)U(slice_all, i + 1));
    std::vector<ca::MX> input_vector(3);
    input_vector[0] = (ca::MX)X(slice_all, i - 1);
    input_vector[1] = (ca::MX)X(slice_all, i);
    input_vector[2] = (ca::MX)X(slice_all, i + 1);
    constraint_vector.push_back(f_d_nm1 + f_d_n +
                                d_EulerLagrange_function(input_vector).at(0));
  }

  // initialize condition
  ca::MX f_0 =
      discrete_forces(dt_, external_force_function, (ca::MX)X(slice_all, 0),
                      (ca::MX)U(slice_all, 0), (ca::MX)U(slice_all, 1));
  std::vector<ca::MX> input_vector(4);
  input_vector[0] = (ca::MX)X_ref(slice_state_, 0);
  input_vector[1] = (ca::MX)X_ref(ca::Slice(6, 12), 0);
  input_vector[2] = (ca::MX)X(slice_all, 0);
  input_vector[3] = (ca::MX)X(slice_all, 1);
  constraint_vector.push_back(
      d_EulerLagrange_init_function(input_vector).at(0) + f_0);

  number_equality_constraint_ = ca::MX::vertcat(constraint_vector).size().first;

  std::cout << "number of equal constraints: " << number_equality_constraint_
            << std::endl;

  // velocity constraints
  for (int i = 0; i < prediction_horizon_ - 1; ++i) {
    constraint_vector.push_back(average_velocity(
        dt_, (ca::MX)X(slice_all, i), (ca::MX)X(slice_all, i + 1), 6));
  }

  casadi::MXDict nlp = {
      {"x",
       ca::MX::vertcat({ca::MX::reshape(X, -1, 1), ca::MX::reshape(U, -1, 1)})},
      {"f", obj},
      {"p", ca::MX::reshape(X_ref, -1, 1)},
      {"g",
       ca::MX::vertcat(
           constraint_vector)}}; // {"p",
                                 // ca::MX::vertcat({ca::MX::reshape(current_state_ref,
                                 // -1, 1), ca::MX::reshape(X_ref, -1, 1)})},
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

void DMOCUAV::set_boundary(
    const std::vector<double> u_min, const std::vector<double> u_max,
    const std::vector<double> x_min,
    const std::vector<double> x_max, const std::vector<double> v,
    const std::vector<double> d_rpy) {
  for (int i = 0; i < prediction_horizon_; i++) {
    u_min_.insert(u_min_.end(), u_min.begin(), u_min.end());
    u_max_.insert(u_max_.end(), u_max.begin(), u_max.end());
    x_min_.insert(x_min_.end(), x_min.begin(), x_min.end());
    x_max_.insert(x_max_.end(), x_max.begin(), x_max.end());
  }

  for (int i = 0; i < number_equality_constraint_; i++) {
    lbg_.push_back(0.0);
    ubg_.push_back(0.0);
  }

  for (int i = 0; i < prediction_horizon_ - 1; i++) {
    for (int j = 0; j < v.size(); j++) {
      lbg_.push_back(-v[j]);
      ubg_.push_back(v[j]);
    }
    for (int j = 0; j < d_rpy.size(); j++) {
      lbg_.push_back(-d_rpy[j]);
      ubg_.push_back(d_rpy[j]);
    }
  }
}

void DMOCUAV::get_results(std::vector<double> init_value,
                          std::vector<double> desired_trajectory,
                          Eigen::MatrixXd &result_x_matrix,
                          Eigen::MatrixXd &result_u_matrix) {
  std::vector<double> lbx;
  std::vector<double> ubx;

  lbx.insert(lbx.end(), x_min_.begin(), x_min_.end());
  lbx.insert(lbx.end(), u_min_.begin(), u_min_.end());
  ubx.insert(ubx.end(), x_max_.begin(), x_max_.end());
  ubx.insert(ubx.end(), u_max_.begin(), u_max_.end());

  std::cout<< lbx.size() << ":" << ubx.size() << std::endl;

  ca::DMDict arg = {{"lbx", lbx},       {"ubx", ubx},
                    {"lbg", lbg_},       {"ubg", ubg_},
                    {"x0", init_value}, {"p", desired_trajectory}};

  opt_results_ = solver_(arg);

  std::vector<double> result_all(opt_results_.at("x"));
  std::vector<double> result_x, result_u;
  result_x.assign(result_all.begin(),
                  result_all.begin() + prediction_horizon_ * num_dofs_);
  result_x_matrix =
      Eigen::MatrixXd::Map(result_x.data(), num_dofs_, prediction_horizon_);

  result_u.assign(result_all.begin() + prediction_horizon_ * num_dofs_,
                  result_all.begin() + prediction_horizon_ * num_dofs_ +
                      +prediction_horizon_  * num_controls_);
  result_u_matrix = Eigen::MatrixXd::Map(result_u.data(), num_controls_,
                                         prediction_horizon_);
}

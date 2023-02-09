/*
 * @Author: Wei Luo
 * @Date: 2023-01-09 17:40:21
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-02-09 17:53:46
 * @Note: Note
 */

#include "mpc_quadrotor.hpp"

MPCQuadrotor::MPCQuadrotor(const double dt, const double N,
                           const double gravity_acceleration,
                           const double roll_tau, const double roll_gain,
                           const double pitch_tau, const double pitch_gain) {
  prediction_horizon_ = N;
  dt_ = dt;
  g_acceleration_ = gravity_acceleration;
  roll_tau_ = roll_tau;
  roll_gain_ = roll_gain;
  pitch_tau_ = pitch_tau;
  pitch_gain_ = pitch_gain;
}

MPCQuadrotor::~MPCQuadrotor() {}

void MPCQuadrotor::initialization_formulation() {
  // control parameter
  ca::MX roll_ref = ca::MX::sym("roll_ref");
  ca::MX pitch_ref = ca::MX::sym("pitch_ref");
  ca::MX thrust_ref = ca::MX::sym("thrust_ref");
  ca::MX controls = ca::MX::vertcat({roll_ref, pitch_ref, thrust_ref});
  num_controls = controls.size().first;

  ca::MX x = ca::MX::sym("x");
  ca::MX y = ca::MX::sym("y");
  ca::MX z = ca::MX::sym("z");
  ca::MX vx = ca::MX::sym("vx");
  ca::MX vy = ca::MX::sym("vy");
  ca::MX vz = ca::MX::sym("vz");
  ca::MX roll = ca::MX::sym("roll");
  ca::MX pitch = ca::MX::sym("pitch");
  ca::MX yaw = ca::MX::sym("yaw");
  ca::MX q = ca::MX::vertcat({x, y, z, vx, vy, vz, roll, pitch, yaw});
  num_states = q.size().first;

  std::vector<ca::MX> rhs;
  rhs.push_back(vx);
  rhs.push_back(vy);
  rhs.push_back(vz);
  rhs.push_back((ca::MX::cos(roll) * ca::MX::cos(yaw) * ca::MX::sin(pitch) +
                 ca::MX::sin(roll) * ca::MX::sin(yaw)) *
                thrust_ref);
  rhs.push_back((ca::MX::cos(roll) * ca::MX::sin(pitch) * ca::MX::sin(yaw) -
                 ca::MX::cos(yaw) * ca::MX::sin(roll)) *
                thrust_ref);
  rhs.push_back(
      (-g_acceleration_ + ca::MX::cos(pitch) * ca::MX::cos(roll) * thrust_ref));
  rhs.push_back((roll_gain_ * roll_ref - roll) / roll_tau_);
  rhs.push_back((pitch_gain_ * pitch_ref - pitch) / pitch_tau_);
  rhs.push_back(0.0);

  system_dynamics_ =
      ca::Function("dynamics", {q, controls}, {ca::MX::vertcat(rhs)});

  // MPC
  ca::MX U = ca::MX::sym("U", num_controls, prediction_horizon_ - 1);
  ca::MX X = ca::MX::sym("X", num_states, prediction_horizon_);
  ca::MX X_ref = ca::MX::sym("X_ref", num_states, prediction_horizon_);
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
  ca::DM R_m = ca::DM::zeros(3, 3);
  R_m(0, 0) = 50.0;
  R_m(1, 1) = 60.0;
  R_m(2, 2) = 1.0;
  ca::DM Q_m = ca::DM::zeros(8, 8);
  Q_m(0, 0) = 80.0;
  Q_m(1, 1) = 80.0;
  Q_m(2, 2) = 120.0;
  Q_m(3, 3) = 80.0;
  Q_m(4, 4) = 80.0;
  Q_m(5, 5) = 100.0;
  Q_m(6, 6) = 10.0;
  Q_m(7, 7) = 10.0;

  ca::MX obj;
  // end term
  obj = ca::MX::mtimes({(X(slice_state, -1) - X_ref(slice_state, -1)).T(), P_m,
                        X(slice_state, -1) - X_ref(slice_state, -1)});
  // control cost
  for (int i = 0; i < prediction_horizon_ - 1; i++) {
    ca::MX temp_ =
        ca::MX::vertcat({U(ca::Slice(0, 2), i),
                         ca::MX::cos(X(6, i)) * ca::MX::cos(X(7, i)) * U(2, i) -
                             g_acceleration_});
    obj += ca::MX::mtimes({temp_.T(), R_m, temp_});
  }
  // state cost
  for (int i = 0; i < prediction_horizon_ - 1; i++) {
    ca::MX temp_ = X(ca::Slice(0, -1), i) - X_ref(ca::Slice(0, -1), i + 1);
    obj += ca::MX::mtimes({temp_.T(), Q_m, temp_});
  }

  // constraints
  std::vector<ca::MX> constraint_vec;
  constraint_vec.push_back(X(slice_all, 0) - X_ref(slice_all, 0));

  for (int i = 0; i < prediction_horizon_ - 1; i++) {
    ca::MX temp_ = RK4_function((ca::MX)X(slice_all, i),
                                (ca::MX)U(slice_all, i), system_dynamics_);
    constraint_vec.push_back(temp_ - X(slice_all, i + 1));
  }

  casadi::MXDict nlp = {{"x", ca::MX::vertcat({ca::MX::reshape(X, -1, 1),
                                               ca::MX::reshape(U, -1, 1)})},
                        {"f", obj},
                        {"p", ca::MX::reshape(X_ref, -1, 1)},
                        {"g", ca::MX::vertcat(constraint_vec)}};
  std::string solver_name = "ipopt";
  casadi::Dict solver_opts;
  solver_opts["expand"] = true;
  solver_opts["ipopt.max_iter"] = 100;
  solver_opts["ipopt.print_level"] = 0;
  solver_opts["print_time"] = 0;
  solver_opts["ipopt.acceptable_tol"] = 1e-8;
  solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
  // std::string solver_name = "snopt";
  // casadi::Dict solver_opts;

  solver_ = casadi::nlpsol("nlpsol", solver_name, nlp, solver_opts);
}

void MPCQuadrotor::set_boundary(const std::vector<double> u_min,
                                const std::vector<double> u_max,
                                const std::vector<double> x_min,
                                const std::vector<double> x_max) {
  for (int i = 0; i < prediction_horizon_ - 1; i++) {
    u_min_.insert(u_min_.end(), u_min.begin(), u_min.end());
    u_max_.insert(u_max_.end(), u_max.begin(), u_max.end());
  }

  for (int i = 0; i < prediction_horizon_; i++) {
    x_min_.insert(x_min_.end(), x_min.begin(), x_min.end());
    x_max_.insert(x_max_.end(), x_max.begin(), x_max.end());
  }
}

void MPCQuadrotor::get_results(std::vector<double> init_value,
                               std::vector<double> desired_trajectory,
                               Eigen::MatrixXd &result_x_matrix,
                               Eigen::MatrixXd &result_u_matrix) {
  std::vector<double> lbx;
  std::vector<double> ubx;

  lbx.insert(lbx.end(), x_min_.begin(), x_min_.end());
  lbx.insert(lbx.end(), u_min_.begin(), u_min_.end());
  ubx.insert(ubx.end(), x_max_.begin(), x_max_.end());
  ubx.insert(ubx.end(), u_max_.begin(), u_max_.end());

  ca::DMDict arg = {{"lbx", lbx},       {"ubx", ubx},
                    {"lbg", 0.0},       {"ubg", 0.0},
                    {"x0", init_value}, {"p", desired_trajectory}};
  // std::cout << desired_trajectory << std::endl;
  opt_results_ = solver_(arg);

  std::vector<double> result_all(opt_results_.at("x"));
  std::vector<double> result_x, result_u;
  result_x.assign(result_all.begin(),
                  result_all.begin() + prediction_horizon_ * num_states);
  result_x_matrix =
      Eigen::MatrixXd::Map(result_x.data(), num_states, prediction_horizon_);

  result_u.assign(result_all.begin() + prediction_horizon_ * num_states,
                  result_all.begin() + prediction_horizon_ * num_states +
                      +(prediction_horizon_ - 1) * num_controls);
  result_u_matrix = Eigen::MatrixXd::Map(result_u.data(), num_controls,
                                         prediction_horizon_ - 1);
}

void MPCQuadrotor::model_based_movement(Eigen::VectorXd &state,
                                        Eigen::VectorXd control,
                                        Eigen::MatrixXd &guessed_state,
                                        Eigen::MatrixXd &guessed_control) {
  Eigen::VectorXd k1 =
      dyn_function_eigen(state, control);
  Eigen::VectorXd k2 = dyn_function_eigen(state + dt_ / 2.0 * k1, control);
  Eigen::VectorXd k3 = dyn_function_eigen(state + dt_ / 2.0 * k2, control);
  Eigen::VectorXd k4 = dyn_function_eigen(state + dt_ * k3, control);
  state +=
      dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  // std::cout << state.transpose() << std::endl;
  guessed_state(Eigen::all, Eigen::seqN(0, guessed_state.cols() - 1)) =
      guessed_state(Eigen::all, Eigen::seqN(1, guessed_state.cols() ));
  guessed_control(Eigen::all, Eigen::seqN(0, guessed_control.cols() - 1)) =
      guessed_control(Eigen::all, Eigen::seqN(1, guessed_control.cols()));
}
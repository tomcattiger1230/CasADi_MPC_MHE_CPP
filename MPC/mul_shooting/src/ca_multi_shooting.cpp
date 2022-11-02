/*
 * @Author: Wei Luo
 * @Date: 2022-10-31 14:49:55
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-11-02 17:26:24
 * @Note: Note
 */
#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <chrono>
#include <iostream>

#define PI 3.1415926

casadi::Function one_step_integration() {
  auto x = casadi::SX::sym("x");
  auto y = casadi::SX::sym("y");
  auto theta = casadi::SX::sym("theta");
  auto states = casadi::SX::vertcat({x, y, theta});
  auto v = casadi::SX::sym("v");
  auto omega = casadi::SX::sym("omega");
  auto controls = casadi::SX::vertcat({v, omega});

  auto rhs = casadi::SX::vertcat(
      {v * casadi::SX::cos(theta), v * casadi::SX::sin(theta)});
  rhs = casadi::SX::vertcat({rhs, omega});

  return casadi::Function("integrator", {states, controls}, {rhs});
}

std::vector<double> one_step_integration_std(double v, double omega,
                                             double theta) {
  std::vector<double> result_;
  result_.push_back(v * std::cos(theta));
  result_.push_back(v * std::sin(theta));
  result_.push_back(omega);
  return result_;
}

Eigen::Vector3d shift_movement(double dt, double &t, std::vector<double> x0,
                               std::vector<double> &u) {
  std::vector<double> dx;
  dx = one_step_integration_std(u[0], u[1], x0[2]);
  Eigen::Vector3d x_e(x0.data());
  Eigen::Vector3d dx_e(dx.data());
  x_e = x_e + dt * dx_e;
  u.erase(u.begin(), u.begin()+2);
  u.push_back(*(u.end() - 1));
  u.push_back(*(u.end() - 1));
  t = t + dt;
  return x_e;
}

int main() {
  double dT = 0.2;             // sampling time [s]
  int N = 100;                 // prediction steps
  double rob_diam = 0.3;       // [m]
  double velocity_max = 0.6;   // max translation velocity m/s
  double omega_max = PI / 4.0; // max rotation
                               //
  auto x = casadi::SX::sym("x");
  auto y = casadi::SX::sym("y");
  auto theta = casadi::SX::sym("theta");
  auto states = casadi::SX::vertcat({x, y, theta});
  int num_states = states.size().first;

  auto v = casadi::SX::sym("v");
  auto omega = casadi::SX::sym("omega");
  auto controls = casadi::SX::vertcat({v, omega});
  int num_controls = controls.size().first;

  // for MPC
  auto U = casadi::SX::sym("U", num_controls, N);
  auto X = casadi::SX::sym("X", num_states, N + 1);
  auto P = casadi::SX::sym("P", num_states * 2);

  casadi::Function integrator = one_step_integration();

  casadi::DM Q(3, 3);
  Q(0, 0) = 1.0;
  Q(1, 1) = 5.0;
  Q(2, 2) = 0.1;

  casadi::DM R(2, 2);
  R(0, 0) = 0.5;
  R(1, 1) = 0.05;

  // cost function
  casadi::SX obj = 0;
  // constraints
  casadi::SX g;
  g = X(casadi::Slice(), 0) - P(casadi::Slice(0, num_states));
  g = casadi::SX::reshape(g, -1, 1);
  //   std::cout << casadi::SX::reshape(X, -1, 1) << std::endl;
  for (int i = 0; i < N; i++) {
    // std::cout << U.nz(casadi::Slice(i * num_controls, (i + 1) * num_controls))
    //           << std::endl;
    obj = obj +
          // casadi::SX::mtimes(
          //     casadi::SX::mtimes(
          //         (X.nz(casadi::Slice(i * num_states, (i + 1) * num_states)) -
          //          P.nz(casadi::Slice(num_states, 2 * num_states)))
          //             .T(),
          //         Q),
          //     (X.nz(casadi::Slice(i * num_states, (i + 1) * num_states)) -
          //      P.nz(casadi::Slice(num_states, 2 * num_states)))) +
          casadi::SX::mtimes(
              casadi::SX::mtimes(
                  (X(casadi::Slice(), i) -
                   P(casadi::Slice(num_states, 2 * num_states)))
                      .T(),
                  Q),
              (X(casadi::Slice(), i) -
               P(casadi::Slice(num_states, 2 * num_states)))) +
          casadi::SX::mtimes(
              casadi::SX::mtimes(
                  U(casadi::Slice(), i)
                      .T(),
                  R),
              U(casadi::Slice(), i));
    std::vector<casadi::SX> input(2);
    input[0] = X(casadi::Slice(), i);
    input[1] = U(casadi::Slice(), i);
    auto x_next_ = integrator(input).at(0) * dT +
                   X(casadi::Slice(), i);
    g = casadi::SX::vertcat(
        {g, x_next_ - X(casadi::Slice(), i+1)});
  }
  casadi::SXDict nlp = {
      {"x", casadi::SX::vertcat({casadi::SX::reshape(X, -1, 1),
                                 casadi::SX::reshape(U, -1, 1)})},
      {"f", obj},
      {"p", P},
      {"g", g}};
  std::string solver_name = "ipopt";
  casadi::Dict solver_opts;
  solver_opts["expand"] = true;
  solver_opts["ipopt.max_iter"] = 100;
  solver_opts["ipopt.print_level"] = 0;
  // solver_opts["linear_solver"] = "ma57";
  solver_opts["print_time"] = 0;
  solver_opts["ipopt.acceptable_tol"] = 1e-8;
  solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;

  casadi::Function solver =
      casadi::nlpsol("nlpsol", solver_name, nlp, solver_opts);

  // constraints definition
  std::vector<double> lbg;
  std::vector<double> ubg;
  std::vector<double> lbx;
  std::vector<double> ubx;

  for (int i = 0; i < N + 1; i++) {
    lbg.push_back(0.0);
    lbg.push_back(0.0);
    lbg.push_back(0.0);
    ubg.push_back(0.0);
    ubg.push_back(0.0);
    ubg.push_back(0.0);
  }

  for (int i = 0; i < N + 1; i++) {
    lbx.push_back(-2.0);
    lbx.push_back(-2.0);
    lbx.push_back(-casadi::inf);
    ubx.push_back(2.0);
    ubx.push_back(2.0);
    ubx.push_back(casadi::inf);
  }

  for (int i = 0; i < N; i++) {
    lbx.push_back(-velocity_max);
    lbx.push_back(-omega_max);
    ubx.push_back(velocity_max);
    ubx.push_back(omega_max);
  }

  std::vector<double> x0{0.0, 0.0, 0.0};
  std::vector<double> xs{1.5, 1.5, 0.0};
  std::vector<double> control_params;
  std::vector<double> init_values;
  std::vector<double> u_init(N * num_controls, 0.0);
  std::vector<double> x_init((N + 1) * num_states, 0.0);

  std::map<std::string, casadi::DM> res;

  double distance_error = 1e4;
  int mpc_iter = 0;
  double sim_time = 20.0;
  int sim_iter = int(sim_time / dT);
  auto start_time = std::chrono::high_resolution_clock::now();
  double time_t = 0.0;
  Eigen::Vector3d target_x;
  target_x << xs[0], xs[1], xs[2];

  while (distance_error > 1e-2 && mpc_iter - sim_iter < 0) {
    control_params.clear();
    init_values.clear();
    // set parameters
    std::merge(x0.begin(), x0.end(), xs.begin(), xs.end(),
               std::back_inserter(control_params));
    // guess solution
    std::merge(x_init.begin(), x_init.end(), u_init.begin(), u_init.end(),
               std::back_inserter(init_values));
    casadi::DMDict arg = {{"lbx", lbx},          {"ubx", ubx},
                          {"lbg", lbg},          {"ubg", ubg},
                          {"p", control_params}, {"x0", init_values}};
    // auto start_time_i = std::chrono::high_resolution_clock::now();
    res = solver(arg);
    // auto stop_time_i = std::chrono::high_resolution_clock::now();
    // auto duration_i =
    //     duration_cast<std::chrono::microseconds>(stop_time_i - start_time_i);
    // time_vector_test.push_back(duration_i.count()/1e6);
    std::vector<double> result_all(res.at("x"));
    std::vector<double> result_x, result_u;
    result_x.assign(result_all.begin(),
                    result_all.begin() + (N + 1) * num_states);
    result_u.assign(result_all.begin() + (N + 1) * num_states,
                    result_all.end());

    // std::cout<< "result x: "<< result_x << std::endl;
    // std::cout<<
    // "================================================================" <<
    // std::endl;
    // std::cout<< "result u: "<< result_u << std::endl;
    auto new_x = shift_movement(dT, time_t, x0, result_u);

    distance_error = (new_x - target_x).norm();
    result_x.erase(result_x.begin(), result_x.begin() + num_states);
    for (int i = 0; i < num_states; i++) {
      result_x.push_back(*(result_x.end() - num_states + 1));
    }
    x0[0] = new_x[0];
    x0[1] = new_x[1];
    x0[2] = new_x[2];
    x_init = result_x;
    u_init = result_u;
    mpc_iter++;
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
  std::cout << "average calculation time for each iteration [s]: "
            << duration.count() / mpc_iter / 1e6 << std::endl;
  return 0;
}

/*
 * @Author: Wei Luo
 * @Date: 2022-10-31 14:49:55
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-10-31 18:45:14
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

void shift_movement(double dt, double t0, std::vector<double> x0, std::vector<double> u, casadi::Function f){
  std::vector<casadi::SX> input(2);
  input(0) = x0;
  input(1) = u[];
  auto x_next = f(input).at(0)
}

int main() {
  double T = 0.2;              // sampling time [s]
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
  g = X.nz(casadi::Slice(0, num_states)) - P.nz(casadi::Slice(0, num_states));
  //   std::cout << casadi::SX::reshape(X, -1, 1) << std::endl;
  for (int i = 0; i < N; i++) {
    obj = obj +
          casadi::SX::mtimes(
              casadi::SX::mtimes(
                  (X.nz(casadi::Slice(i * num_states, (i + 1) * num_states)) -
                   P.nz(casadi::Slice(num_states, 2 * num_states)))
                      .T(),
                  Q),
              (X.nz(casadi::Slice(i * num_states, (i + 1) * num_states)) -
               P.nz(casadi::Slice(num_states, 2 * num_states)))) +
          casadi::SX::mtimes(
              casadi::SX::mtimes(
                  U.nz(casadi::Slice(i * num_controls, (i + 1) * num_controls))
                      .T(),
                  R),
              U.nz(casadi::Slice(i * num_controls, (i + 1) * num_controls)));
    std::vector<casadi::SX> input(2);
    input[0] = X.nz(casadi::Slice(i * num_states, (i + 1) * num_states));
    input[1] = U.nz(casadi::Slice(i * num_controls, (i + 1) * num_controls));
    auto x_next_ = integrator(input).at(0) * T +
                   X.nz(casadi::Slice(i * num_states, (i + 1) * num_states));
    g = casadi::SX::vertcat(
        {g, x_next_ - X.nz(casadi::Slice((i + 1) * num_states,
                                         (i + 2) * num_states))});
  }
  casadi::SXDict nlp = {
      {"x", casadi::SX::vertcat({casadi::SX::reshape(X, -1, 1),
                                 casadi::SX::reshape(U, -1, 1)})},
      {"f", obj},
      {"p", P},
      {"g", g}};
  std::string solver_name = "ipopt";
  casadi::Dict solver_opts;
  solver_opts["ipopt.max_iter"] = 100;
  solver_opts["ipopt.print_level"] = 0;
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
    lbx.push_back(-PI);
    ubx.push_back(2.0);
    ubx.push_back(2.0);
    ubx.push_back(PI);
  }

  for (int i = 0; i < N; i++) {
    lbx.push_back(-velocity_max);
    lbx.push_back(-omega_max);
    ubx.push_back(velocity_max);
    ubx.push_back(omega_max);
  }

  std::vector<double> x0{0.0, 0.0, 0.0};
  std::vector<double> xs{1.5, 1.5, 0.0};
//   std::vector<double> control_params;
//   std::vector<double> init_values;
  std::vector<double> u_init(N*num_controls, 0.0);
  std::vector<double> x_init((N + 1) * num_states, 0.0);
//   std::merge(x0.begin(), x0.end(), xs.begin(), xs.end(),
//              std::back_inserter(control_params));
//   std::merge(x_init.begin(), x_init.end(), u_init.begin(), u_init.end(),
//              std::back_inserter(init_values));

  std::map<std::string, casadi::DM> res;
//   casadi::DMDict arg = {{"lbx", lbx},          {"ubx", ubx},
//                         {"lbg", lbg},          {"ubg", ubg},
//                         {"p", control_params}, {"x0", init_values}};


    // std::cout << "x solution: " << result_x << std::endl;
    // std::cout << "u solution: " << result_u << std::endl;

  double distance_error = 1e4;
  int mpc_iter = 0;
  double sim_time = 20.0;
  int sim_iter = int(sim_time / T);
  auto start_time = std::chrono::high_resolution_clock::now();
  while (distance_error > 1e-2 && mpc_iter - sim_iter < 0) {
  std::vector<double> control_params;
  std::vector<double> init_values;
    // set parameters
    std::merge(x0.begin(), x0.end(), xs.begin(), xs.end(),
               std::back_inserter(control_params));
    // guess solution
    std::merge(x_init.begin(), x_init.end(), u_init.begin(), u_init.end(),
               std::back_inserter(init_values));
    casadi::DMDict arg = {{"lbx", lbx},          {"ubx", ubx},
                          {"lbg", lbg},          {"ubg", ubg},
                          {"p", control_params}, {"x0", init_values}};

    res = solver(arg);
    std::vector<double> result_all(res.at("x"));
    std::vector<double> result_x, result_u;
    result_x.assign(result_all.begin(),
                    result_all.begin() + (N + 1) * num_states);
    result_u.assign(result_all.begin() + (N + 1) * num_states,
                    result_all.end());
    mpc_iter += 1;
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(stop_time - start_time);
  std::cout << "average calculation time for each iteration [s]: " << duration.count() / mpc_iter /1e6 << std::endl;
  return 0;
}

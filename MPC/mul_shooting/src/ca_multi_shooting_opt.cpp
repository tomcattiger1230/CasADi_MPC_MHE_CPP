/*
 * @Author: Wei Luo
 * @Date: 2022-11-01 21:34:03
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-12-12 08:11:15
 * @Note: Note
 */

#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <ctime>
#include <iostream>

#define PI 3.1415926

casadi::Function one_step_integration() {
  auto x = casadi::MX::sym("x");
  auto y = casadi::MX::sym("y");
  auto theta = casadi::MX::sym("theta");
  auto states = casadi::MX::vertcat({x, y, theta});
  auto v = casadi::MX::sym("v");
  auto omega = casadi::MX::sym("omega");
  auto controls = casadi::MX::vertcat({v, omega});

  auto rhs = casadi::MX::vertcat(
      {v * casadi::MX::cos(theta), v * casadi::MX::sin(theta)});
  rhs = casadi::MX::vertcat({rhs, omega});

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

void shift_movement(double dt, double &t, std::vector<double> &x0,
                    std::vector<double> & u) {
  std::vector<double> dx;
  dx = one_step_integration_std(u[0], u[1], x0[2]);
  Eigen::Vector3d x_e(x0.data());
  Eigen::Vector3d dx_e(dx.data());
  x_e = x_e + dt * dx_e;
  x0[0] = x_e[0];
  x0[1] = x_e[1];
  x0[2] = x_e[2];
  u.erase(u.begin(), u.begin() + 2);
  u.push_back(*(u.end() - 1));
  u.push_back(*(u.end() - 1));
  t = t + dt;
}

int main() {
  double dT = 0.2; // sampling time [s]
  const int N = 100;
  double velocity_max = 0.6;   // max translation velocity m/s
  double omega_max = PI / 4.0; // max rotation

  auto opti = casadi::Opti();
  casadi::Slice all;

  auto opt_states = opti.variable(N + 1, 3);
  auto opt_controls = opti.variable(N, 2);

  auto v = opt_controls(all, 0);
  auto omega = opt_controls(all, 1);
  auto x = opt_states(all, 0);
  auto y = opt_states(all, 1);
  auto theta = opt_states(all, 2);

  casadi::Function integrator_casadi = one_step_integration();

  casadi::MX opt_x0 = opti.parameter(3);
  casadi::MX opt_xs = opti.parameter(3);

  casadi::MX obj = 0;
  casadi::DM Q(3, 3);
  Q(0, 0) = 1.0;
  Q(1, 1) = 5.0;
  Q(2, 2) = 0.1;

  casadi::DM R(2, 2);
  R(0, 0) = 0.5;
  R(1, 1) = 0.05;

  for (int i = 0; i < N; i++) {
    obj = obj +
          casadi::MX::mtimes(
              casadi::MX::mtimes((opt_states(i, all) - opt_xs.T()), Q),
              (opt_states(i, all) - opt_xs.T()).T()) +
          casadi::MX::mtimes(casadi::MX::mtimes((opt_controls(i, all)), R),
                             opt_controls(i, all).T());
  }
  opti.minimize(obj);

  opti.subject_to(opt_states(0, all) == opt_x0.T());
  for (int i = 0; i < N; i++) {
    std::vector<casadi::MX> input(2);
    input[0] = opt_states(i, all);
    input[1] = opt_controls(i, all);
    auto x_next_ = integrator_casadi(input).at(0) * dT + opt_states(i, all).T();
    opti.subject_to(x_next_ == opt_states(i+1, all).T());
  }

  // boundary conditions
  opti.subject_to(opti.bounded (- 2.0 , x , 2.0));
  opti.subject_to(opti.bounded(-2.0, y, 2.0));
  opti.subject_to(opti.bounded(-velocity_max,v ,velocity_max));
  opti.subject_to(opti.bounded(-omega_max , omega, omega_max));

  std::string solver_name = "ipopt";
  casadi::Dict solver_opts;
  solver_opts["expand"] = true;
  solver_opts["ipopt.max_iter"] = 100;
  solver_opts["ipopt.print_level"] = 0;
  solver_opts["print_time"] = 0;
  solver_opts["ipopt.acceptable_tol"] = 1e-8;
  solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;

  opti.solver(solver_name, solver_opts);

  std::vector<double> x0{0.0, 0.0, 0.0};
  std::vector<double> xs{1.5, 1.5, 0.0};
  std::vector<double> control_params;
  std::vector<double> init_values;
  Eigen::Matrix<double, N, 2> u_init = Eigen::Matrix<double, N, 2>::Zero();
  Eigen::Matrix<double, N + 1, 3> x_init =
      Eigen::Matrix<double, N + 1, 3>::Zero();


  double distance_error = 1e4;
  int mpc_iter = 0;
  double sim_time = 20.0;
  int sim_iter = int(sim_time / dT);
  auto start_time = std::chrono::high_resolution_clock::now();
  double time_t = 0.0;

  while (distance_error > 1e-2 && mpc_iter - sim_iter < 0) {
    casadi::DM xs_{xs};
    casadi::DM x0_{x0};
    opti.set_value(opt_xs, xs_);
    opti.set_value(opt_x0, x0_);
    casadi::DM u_init_{
        std::vector<double>(u_init.data(), u_init.size() + u_init.data())};
    u_init_ = casadi::DM::reshape(u_init_, u_init.rows(), u_init.cols());
    casadi::DM x_init_{
        std::vector<double>(x_init.data(), x_init.size() + x_init.data())};
    x_init_ = casadi::DM::reshape(x_init_, x_init.rows(), x_init.cols());
    opti.set_initial(opt_controls, u_init_);
    opti.set_initial(opt_states, x_init_);
    auto res = opti.solve();

    auto result_x = res.value(opt_states);
    auto result_u = res.value(opt_controls);
    u_init = Eigen::Matrix<double, N, 2>(
        static_cast<std::vector<double>>(result_u).data());
    x_init = Eigen::Matrix<double, N + 1, 3>(
        static_cast<std::vector<double>>(result_x).data());
    std::vector<double> u_init_std(
        &u_init(0), u_init.data() + u_init.cols() * u_init.rows());

    shift_movement(dT, time_t, x0, u_init_std);
    mpc_iter++;
  }
  auto stop_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
      stop_time - start_time);
  std::cout << "average calculation time for each iteration [s]: "
            << duration.count() / mpc_iter / 1e6 << std::endl;

  return 0;
}
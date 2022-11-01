/*
 * @Author: Wei Luo
 * @Date: 2022-11-01 21:34:03
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-11-01 21:37:55
 * @Note: Note
 */

#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <ctime>
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
  u.erase(u.begin(), u.begin() + 2);
  u.push_back(*(u.end() - 1));
  u.push_back(*(u.end() - 1));
  t = t + dt;
  return x_e;
}

int main() {
  double dT = 0.2; // sampling time [s]
  int N = 100;

  auto opti = casadi :: Opti();
  casadi::Slice all;

  auto opt_states = opti.variable(3, N + 1);
//   auto Y =

  return 0;
}
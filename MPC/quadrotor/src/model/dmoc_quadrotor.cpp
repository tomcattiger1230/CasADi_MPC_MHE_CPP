/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:34:11
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-17 21:00:19
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
  // MPC
  ca::MX U = ca::MX::sym("U", num_controls_, prediction_horizon_);
  ca::MX X = ca::MX::sym("X", num_dofs_, prediction_horizon_);
  ca::MX X_ref = ca::MX::sym("X_ref", num_dofs_ + 3, prediction_horizon_);

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
  ca::DM ref_u = ca::DM::zeros(4);
  ref_u(0) = 0.25 * mass_quadrotor_ * g_acceleration_;
  ref_u(1) = 0.25 * mass_quadrotor_ * g_acceleration_;
  ref_u(2) = 0.25 * mass_quadrotor_ * g_acceleration_;
  ref_u(3) = 0.25 * mass_quadrotor_ * g_acceleration_;

  ca::MX obj;

  // control cost
  for (int i = 0; i < prediction_horizon_; i++) {
    ca::MX temp_ = U(ca::Slice(0, 4), i) - ref_u;
    obj += ca::MX::mtimes({temp_.T(), R_m, temp_});
  }
  // state cost
  for (int i = 0; i < prediction_horizon_; i++) {
    obj = ca::MX::mtimes({(X(slice_state_, i) - X_ref(slice_state_, i)).T(),
                          P_m, X(slice_state_, i) - X_ref(slice_state_, i)});
  }
}
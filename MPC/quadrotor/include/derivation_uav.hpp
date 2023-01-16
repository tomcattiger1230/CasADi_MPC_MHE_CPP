/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:29:49
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-16 23:43:36
 * @Note: Note
 */

#ifndef __DERIVATION_UAV__
#define __DERIVATION_UAV__

#include "dmcc_base.hpp"
#include <iostream>

namespace ca = casadi;

class DerivationUAV : public DMCCBase {
public:
  DerivationUAV(const double g, const double mass_quadrotor,
                const std::vector<double> I_quadrotor, const double frame_size,
                const double motor_torque_const) {
    g_acceleration_ = g;
    // I_quadrotor_ = Eigen::Vector3d::Map(I_quadrotor.data(), 3);
    I_quadrotor_ = I_quadrotor;
    frame_size_ = frame_size;
    motor_torque_const_ = motor_torque_const;

    get_Lagrangian_casadi();
  }

private:
  double g_acceleration_;
  double mass_quadrotor_;
//   Eigen::Vector3d I_quadrotor_;
  std::vector<double> I_quadrotor_;
  double frame_size_;
  double motor_torque_const_;

  int num_dofs_;
  int num_states_;

  void get_Lagrangian_casadi(){
    // state
    ca::MX x = ca::MX::sym("x");
    ca::MX y = ca::MX::sym("y");
    ca::MX z = ca::MX::sym("z");
    ca::MX phi = ca::MX::sym("phi");
    ca::MX theta = ca::MX::sym("theta");
    ca::MX psi = ca::MX::sym("psi");

    ca::MX q = ca::MX::vertcat({x, y, z, phi, theta, psi});
    num_dofs_ = q.size().first;
    ca::MX position_vector = ca::MX::vertcat({x, y, z});

    // state dot
    ca::MX d_x = ca::MX::sym("d_x");
    ca::MX d_y = ca::MX::sym("d_y");
    ca::MX d_z = ca::MX::sym("d_z");
    ca::MX d_phi = ca::MX::sym("d_phi");
    ca::MX d_theta = ca::MX::sym("d_theta");
    ca::MX d_psi = ca::MX::sym("d_psi");
    ca::MX d_q = ca::MX::vertcat({d_x, d_y, d_z, d_phi, d_theta, d_psi});

    ca::MX full_state = ca::MX::vertcat({q, d_q});
    num_states_ = full_state.size().first;
    ca::MX d_position_vector = ca::MX::vertcat({d_x, d_y, d_z});

    ca::MX Rz = ca::MX::zeros(3, 3);
    Rz(0, 0) = ca::MX::cos(psi);
    Rz(0, 1) = -ca::MX::sin(psi);
    Rz(1, 0) = ca::MX::sin(psi);
    Rz(1, 1) = ca::MX::cos(psi);
    Rz(2, 2) = 1.0;

    ca::MX Ry = ca::MX::zeros(3, 3);
    Ry(0, 0) = ca::MX::cos(theta);
    Ry(0, 2) = ca::MX::sin(theta);
    Ry(1, 1) = 1.0;
    Ry(2, 0) = -ca::MX::sin(theta);
    Ry(2, 2) = ca::MX::cos(theta);

    ca::MX Rx = ca::MX::zeros(3, 3);
    Rx(0, 0) = 1.0;
    Rx(1, 1) = ca::MX::cos(phi);
    Rx(1, 2) = -ca::MX::sin(phi);
    Rx(2, 1) = ca::MX::sin(phi);
    Rx(2, 2) = ca::MX::cos(phi);

    ca::MX eRb = ca::MX::mtimes({Rz, Ry, Rx});

    ca::DM e3 = ca::DM::zeros(3, 1);
    e3(2, 0) = 1.0;

    ca::MX T_matrix = ca::MX::zeros(3, 3);
    T_matrix(0, 0) = 1.0;
    T_matrix(0, 2) = -ca::MX::sin(theta);
    T_matrix(1, 1) = ca::MX::cos(phi);
    T_matrix(1, 2) = ca::MX::sin(phi) * ca::MX::cos(theta);
    T_matrix(2, 1) = -ca::MX::sin(phi);
    T_matrix(2, 2) = ca::MX::cos(phi) * ca::MX::cos(theta);

    ca::MX bW = ca::MX::mtimes({T_matrix, d_q(ca::Slice(3, 6))});
    ca::DM Ib = ca::DM::zeros(3, 3);
    Ib(0, 0) = I_quadrotor_[0];
    Ib(1, 1) = I_quadrotor_[1];
    Ib(2, 2) = I_quadrotor_[2];

    ca::MX K = 0.5 * mass_quadrotor_ * ca::MX::mtimes({q(ca::Slice(0, 3)).T(), q(ca::Slice(0, 3))}) + 0.5 * ca::MX::mtimes({bW.T(), Ib, bW});
  }
};

#endif
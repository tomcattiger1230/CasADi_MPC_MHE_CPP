/*
 * @Author: Wei Luo
 * @Date: 2022-11-30 15:49:20
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-02-06 07:55:55
 * @Note: Note
 */

#ifndef __DERIVATION_UAV_MANIPULATOR_POS_FORCE__
#define __DERIVATION_UAV_MANIPULATOR_POS_FORCE__

#include <dmcc_base.hpp>
#include <iostream>

namespace ca = casadi;

class DerivationUAVm : public DMCCBase {
public:
  DerivationUAVm(const double g, const double arm_length, const double mass_arm,
                 const double mass_quad, const std::vector<double> I_quadrotor,
                 const std::vector<double> I_arm, const double frame_size,
                 const double motor_torque_const,
                 const std::vector<double> montage_offset_b) {
    // set up parameters
    g_acceleration_ = g;
    arm_length_ = arm_length;
    mass_manipulator_ = mass_arm;
    mass_quadrotor_ = mass_quad;
    I_quadrotor_ = I_quadrotor;
    I_arm_ = I_arm;
    frame_size_ = frame_size;
    motor_torque_const_ = motor_torque_const;
    for (int i = 0; i < 3; ++i) {
      montage_offset_b_(i, 0) = montage_offset_b[i];
    }

    // calculate Lagrangian
    get_Lagrangian_casadi();
  };
  ~DerivationUAVm(){};

  ca::Function get_end_effector_position_function() {
    return end_position_function;
  }

  ca::Function get_end_effector_velocity_function() {
    return end_velocity_function;
  }



private:
  double motor_torque_const_;

  ca::Function end_position_function;
  ca::Function end_velocity_function;

  void get_Lagrangian_casadi() {
    ca::MX x = ca::MX::sym("x");
    ca::MX y = ca::MX::sym("y");
    ca::MX z = ca::MX::sym("z");
    ca::MX phi = ca::MX::sym("phi");
    ca::MX theta = ca::MX::sym("theta");
    ca::MX psi = ca::MX::sym("psi");
    ca::MX alpha = ca::MX::sym("alpha");

    ca::MX q = ca::MX::vertcat({x, y, z, phi, theta, psi, alpha});
    num_dofs_ = q.size().first;
    ca::MX pos = ca::MX::vertcat({x, y, z});

    // state dot
    ca::MX d_x = ca::MX::sym("d_x");
    ca::MX d_y = ca::MX::sym("d_y");
    ca::MX d_z = ca::MX::sym("d_z");
    ca::MX d_phi = ca::MX::sym("d_phi");
    ca::MX d_theta = ca::MX::sym("d_theta");
    ca::MX d_psi = ca::MX::sym("d_psi");
    ca::MX d_alpha = ca::MX::sym("d_alpha");
    ca::MX d_q =
        ca::MX::vertcat({d_x, d_y, d_z, d_phi, d_theta, d_psi, d_alpha});

    ca::MX full_state = ca::MX::vertcat({q, d_q});
    num_states_ = full_state.size().first;
    ca::MX d_pos = ca::MX::vertcat({d_x, d_y, d_z});

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

    ca::MX angle_rate = ca::MX::vertcat({d_phi, d_theta, d_psi});
    ca::MX bW = ca::MX::mtimes({T_matrix, angle_rate});

    ca::MX eW = ca::MX::mtimes({eRb, bW});
    ca::MX skew_eW_ca = ca::MX::skew(eW);

    ca::MX b0 = ca::MX::zeros(3, 1);
    b0(0, 0) = ca::MX::cos(alpha);
    b0(2, 0) = -ca::MX::sin(alpha);

    b0 = b0 * arm_length_ / 2.0 + montage_offset_b_;

    ca::MX e0 = pos + ca::MX::mtimes({eRb, b0});

    ca::MX Jl1 = ca::MX::zeros(3, 1);
    Jl1(0, 0) = -ca::MX::sin(alpha);
    Jl1(2, 0) = -ca::MX::cos(alpha);
    Jl1 = Jl1 * arm_length_ / 2.0;

    ca::MX b0_dot = Jl1 * d_alpha;

    ca::MX e0_dot = d_pos + ca::MX::mtimes({skew_eW_ca, eRb, b0}) +
                    ca::MX::mtimes({eRb, b0_dot});

    ca::MX bRm = ca::MX::zeros(3, 3);
    bRm(0, 0) = ca::MX::cos(alpha);
    bRm(0, 2) = ca::MX::sin(alpha);
    bRm(1, 1) = 1.0;
    bRm(2, 0) = -ca::MX::sin(alpha);
    bRm(2, 2) = ca::MX::cos(alpha);

    ca::DM Jol1 = ca::DM::zeros(3, 1);
    Jol1(1, 0) = 1.0;

    ca::MX bWl1 = Jol1 * d_alpha;
    ca::MX eWl1 = eW + ca::MX::mtimes({eRb, bWl1});

    ca::MX mani_angular_velocity_b = ca::MX::mtimes({bRm.T(), eRb.T(), eWl1});
    ca::MX quadrotor_inertia_moment = ca::DM::diag(I_quadrotor_);
    ca::MX manipulator_inertia_moment = ca::DM::diag(I_arm_);

    // kinetic energy
    ca::MX K_quad =
        0.5 * mass_quadrotor_ * ca::MX::mtimes({d_pos.T(), d_pos}) +
        0.5 * ca::MX::mtimes({bW.T(), quadrotor_inertia_moment, bW});
    ca::MX K_arm =
        0.5 * mass_manipulator_ * ca::MX::mtimes({e0_dot.T(), e0_dot}) +
        0.5 * ca::MX::mtimes({mani_angular_velocity_b.T(),
                              manipulator_inertia_moment,
                              mani_angular_velocity_b});
    ca::MX K_total = K_quad + K_arm;

    // potential energy
    ca::MX U_quad =
        mass_quadrotor_ * g_acceleration_ * ca::MX::mtimes({e3.T(), pos});
    ca::MX U_arm =
        mass_manipulator_ * g_acceleration_ * ca::MX::mtimes({e3.T(), e0});
    ca::MX U_total = U_quad + U_arm;

    ca::MX L = K_total - U_total;

    lagrangian_function_ = ca::Function("function_L", {q, d_q}, {L});

    ca::MX L_d_dot_q = ca::MX::gradient(L, d_q);

    derivative_lagrangian_function_ =
        ca::Function("function_lagrangian_derivative", {q, d_q}, {L_d_dot_q});

    // motor
    ca::MX U1 = ca::MX::sym("U1");
    ca::MX U2 = ca::MX::sym("U2");
    ca::MX U3 = ca::MX::sym("U3");
    ca::MX U4 = ca::MX::sym("U4");
    ca::MX tau_m = ca::MX::sym("tau_m");
    ca::MX input_vec = ca::MX::vertcat({U1, U2, U3, U4, tau_m});

    num_controls_ = input_vec.size().first;

    ca::MX total_force_local = U1 + U2 + U3 + U4;
    ca::MX total_force_global = ca::MX::mtimes({eRb, e3, total_force_local});

    ca::MX moment_x = (U2 + U3 - U1 - U4) * frame_size_ / 2.0 / std::sqrt(2);
    ca::MX moment_y = (U2 + U4 - U1 - U3) * frame_size_ / 2.0 / std::sqrt(2);
    ca::MX moment_z = (U3 + U4 - U1 - U2) * motor_torque_const_;

    ca::MX external_forces = ca::MX::vertcat(
        {total_force_global, moment_x, moment_y - tau_m, moment_z, tau_m});

    external_force_function_ = ca::Function("external_force_function",
                                           {q, input_vec}, {external_forces});

    ca::MX b_Mani_End = ca::MX::zeros(3, 1);
    b_Mani_End(0, 0) = ca::MX::cos(alpha) * arm_length_;
    b_Mani_End(2, 0) = -ca::MX::sin(alpha) * arm_length_;
    b_Mani_End += montage_offset_b_;
    ca::MX e_Mani_End = pos + ca::MX::mtimes({eRb, b_Mani_End});

    end_position_function =
        ca::Function("end_position_function", {q}, {e_Mani_End});

    ca::MX Jl1_Mani_End = ca::MX::zeros(3, 1);
    Jl1_Mani_End(0, 0) = -ca::MX::sin(alpha) * arm_length_;
    Jl1_Mani_End(2, 0) = -ca::MX::cos(alpha) * arm_length_;
    ca::MX b_Mani_End_dot = Jl1_Mani_End * d_alpha;
    ca::MX e_Mani_End_dot = d_pos +
                            ca::MX::mtimes({skew_eW_ca, eRb, b_Mani_End}) +
                            ca::MX::mtimes({eRb, b_Mani_End_dot});

    end_velocity_function =
        ca::Function("end_velocity_function", {q, d_q}, {e_Mani_End_dot});
  }
};
#endif
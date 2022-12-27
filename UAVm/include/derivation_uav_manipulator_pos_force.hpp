/*
 * @Author: Wei Luo
 * @Date: 2022-11-30 15:49:20
 * @LastEditors: Wei Luo
 * @LastEditTime: 2022-12-27 16:16:58
 * @Note: Note
 */

#ifndef __DERIVATION_UAV_MANIPULATOR_POS_FORCE__
#define __DERIVATION_UAV_MANIPULATOR_POS_FORCE__

#include <dmcc_base.hpp>

namespace ca = casadi;

class DerivationUAV : public DMCCBase {
public:
  DerivationUAV(const double g, const double arm_length, const double mass_arm,
                const double mass_quad, const std::vector<double> I_quadrotor,
                const std::vector<double> I_arm, const double frame_size,
                const double motor_torque_const,
                const std::vector<double> montage_offset_b) {
    // set up parameters
    g_ = g;
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
  ~DerivationUAV(){};

  ca::Function get_end_effector_position_function() {
    return end_position_function;
  }

  ca::Function get_end_effector_velocity_function() {
    return end_velocity_function;
  }

  ca::Function get_system_dynamics_function() {
    return system_dynamics_function;
  }

  ca::Function get_lagrangian_function() { return lagrangian_function; }

  ca::Function get_derivative_lagrangian_function() {
    return derivative_lagrangian_function;
  }



private:
  double motor_torque_const_;

  ca::Function end_position_function;
  ca::Function end_velocity_function;

  ca::Function system_dynamics_function;
  ca::Function derivative_lagrangian_function;
  ca::Function lagrangian_function;

  void get_Lagrangian_casadi() {
    auto x = ca::MX::sym("x");
    auto y = ca::MX::sym("y");
    auto z = ca::MX::sym("z");
    auto phi = ca::MX::sym("phi");
    auto theta = ca::MX::sym("theta");
    auto psi = ca::MX::sym("psi");
    auto alpha = ca::MX::sym("alpha");

    ca::MX q = ca::MX::vertcat({x, y, z, phi, theta, psi, alpha});
    num_dofs_ = q.size1();
    auto pos = ca::MX::vertcat({x, y, z});

    // state dot
    auto d_x = ca::MX::sym("d_x");
    auto d_y = ca::MX::sym("d_y");
    auto d_z = ca::MX::sym("d_z");
    auto d_phi = ca::MX::sym("d_phi");
    auto d_theta = ca::MX::sym("d_theta");
    auto d_psi = ca::MX::sym("d_psi");
    auto d_alpha = ca::MX::sym("d_alpha");
    ca::MX d_q =
        ca::MX::vertcat({d_x, d_y, d_z, d_phi, d_theta, d_psi, d_alpha});

    auto full_state = ca::MX::vertcat({q, d_q});
    num_states_ = full_state.size1();
    auto d_pos = ca::MX::vertcat({d_x, d_y, d_z});

    auto Rz = ca::MX(3, 3);
    Rz(0, 0) = ca::MX::cos(psi);
    Rz(0, 1) = -ca::MX::sin(psi);
    Rz(1, 0) = ca::MX::sin(psi);
    Rz(1, 1) = ca::MX::cos(psi);
    Rz(2, 2) = 1.0;

    auto Ry = ca::MX(3, 3);
    Ry(0, 0) = ca::MX::cos(theta);
    Ry(0, 2) = ca::MX::sin(theta);
    Ry(1, 1) = 1.0;
    Ry(2, 0) = -ca::MX::sin(theta);
    Ry(2, 2) = ca::MX::cos(theta);

    auto Rx = ca::MX(3, 3);
    Rx(0, 0) = 1.0;
    Rx(1, 1) = ca::MX::cos(phi);
    Rx(1, 2) = -ca::MX::sin(phi);
    Rx(2, 1) = ca::MX::sin(phi);
    Rx(2, 2) = ca::MX::cos(phi);

    auto eRb = ca::MX::mtimes({Rz, Ry, Rx});

    ca::DM e3(3, 1);
    e3(2, 0) = 1.0;

    auto T_matrix = ca::MX(3, 3);
    T_matrix(0, 0) = 1.0;
    T_matrix(0, 2) = -ca::MX::sin(theta);
    T_matrix(1, 1) = ca::MX::cos(phi);
    T_matrix(1, 2) = ca::MX::sin(phi) * ca::MX::cos(theta);
    T_matrix(2, 1) = -ca::MX::sin(phi);
    T_matrix(2, 2) = ca::MX::cos(phi) * ca::MX::cos(theta);

    auto angle_rate = ca::MX::vertcat({d_phi, d_theta, d_psi});
    auto bW = ca::MX::mtimes({T_matrix, angle_rate});

    auto eW = ca::MX::mtimes({eRb, bW});
    auto skew_eW_ca = ca::MX::skew(eW);

    auto b0 = ca::MX(3, 1);
    b0(0, 0) = ca::MX::cos(alpha);
    b0(1, 0) = 0.0;
    b0(2, 0) = -ca::MX::sin(alpha);

    b0 = b0 * arm_length_ / 2.0 + montage_offset_b_;

    auto e0 = pos + ca::MX::mtimes({eRb, b0});

    auto Jl1 = ca::MX(3, 1);
    Jl1(0, 0) = -ca::MX::sin(alpha);
    Jl1(1, 0) = 0.0;
    Jl1(2, 0) = -ca::MX::cos(alpha);
    Jl1 = Jl1 * arm_length_ / 2.0;

    auto b0_dot = Jl1 * d_alpha;

    auto e0_dot = d_pos + ca::MX::mtimes({skew_eW_ca, eRb, b0}) +
                  ca::MX::mtimes({eRb, b0_dot});

    auto bRm = ca::MX(3, 3);
    bRm(0, 0) = ca::MX::cos(alpha);
    bRm(0, 1) = 0.0;
    bRm(0, 2) = ca::MX::sin(alpha);
    bRm(1, 0) = 0.0;
    bRm(1, 1) = 1.0;
    bRm(1, 2) = 0.0;
    bRm(2, 0) = -ca::MX::sin(alpha);
    bRm(2, 1) = 0.0;
    bRm(2, 2) = ca::MX::cos(alpha);

    auto Jol1 = ca::DM(3, 1);
    Jol1(1, 0) = 1.0;

    auto bWl1 = Jol1 * d_alpha;
    auto eWl1 = eW + ca::MX::mtimes({eRb, bWl1});

    auto mani_angular_velocity_b = ca::MX::mtimes({bRm.T(), eRb.T(), eWl1});
    auto quadrotor_inertial_moment = ca::DM::diag(I_quadrotor_);
    auto manipulator_inertia_moment = ca::DM::diag(I_arm_);

    // kinetic energy
    auto K_quad = 0.5 * ca::MX::mtimes({mass_quadrotor_, d_pos.T(), d_pos}) +
                  0.5 * ca::MX::mtimes({bW.T(), quadrotor_inertial_moment, bW});
    auto K_arm = 0.5 * ca::MX::mtimes({mass_manipulator_, e0_dot.T(), e0_dot}) +
                 0.5 * ca::MX::mtimes({mani_angular_velocity_b.T(),
                                       manipulator_inertia_moment,
                                       mani_angular_velocity_b});
    auto K_total = K_quad + K_arm;

    // potential energy
    auto U_quad = ca::MX::mtimes({mass_quadrotor_, g_, e3.T(), pos});
    auto U_arm = ca::MX::mtimes({mass_manipulator_, g_, e3.T(), e0});
    auto U_total = U_quad + U_arm;

    ca::MX L = K_total - U_total;

    lagrangian_function = ca::Function("function_L", {q, d_q}, {L});

    auto L_d_dot_q = ca::MX::gradient(L, d_q);

    derivative_lagrangian_function =
        ca::Function("function_lagrangian_derivative", {q, d_q}, {L_d_dot_q});

    // motor
    auto U1 = ca::MX::sym("U1");
    auto U2 = ca::MX::sym("U2");
    auto U3 = ca::MX::sym("U3");
    auto U4 = ca::MX::sym("U4");

    auto tau_m = ca::MX::sym("tau_m");
    auto input_vec = ca::MX::vertcat({U1, U2, U3, U4, tau_m});

    num_controls_ = input_vec.size1();

    auto total_force_local = U1 + U2 + U3 + U4;
    auto total_force_global = ca::MX::mtimes({eRb, e3, total_force_local});

    auto moment_x = (U2 + U3 - U1 - U4) / frame_size_ / 2.0 / std::sqrt(2);
    auto moment_y = (U2 + U4 - U1 - U3) / frame_size_ / 2.0 / std::sqrt(2);
    auto moment_z = (U3 + U4 - U1 - U2) / motor_torque_const_;

    auto dynamics_rhs = ca::MX::vertcat(
        {total_force_global, moment_x, moment_y - tau_m, moment_z, tau_m});

    system_dynamics_function = ca::Function("system_dynamics_function",
                                            {q, input_vec}, {dynamics_rhs});

    auto b_Mani_End = ca::MX(3, 1);
    b_Mani_End(0, 0) = ca::MX::cos(alpha) * arm_length_;
    b_Mani_End(1, 0) = 0.0;
    b_Mani_End(2, 0) = -ca::MX::sin(alpha) * arm_length_;
    b_Mani_End += montage_offset_b_;

    end_position_function =
        ca::Function("end_position_function", {q}, {b_Mani_End});

    auto Jl1_Mani_End = ca::MX(3, 1);
    Jl1_Mani_End(0, 0) = -ca::MX::sin(alpha) * arm_length_;
    Jl1_Mani_End(1, 0) = 0.0;
    Jl1_Mani_End(2, 0) = -ca::MX::cos(alpha) * arm_length_;
    auto b_Mani_End_dot = Jl1_Mani_End * d_alpha;
    auto e_Mani_End_dot = d_pos +
                          ca::MX::mtimes({skew_eW_ca, eRb, b_Mani_End}) +
                          ca::MX::mtimes({eRb, b_Mani_End_dot});

    end_velocity_function =
        ca::Function("end_velocity_function", {q, d_q}, {e_Mani_End_dot});
  }
};
#endif
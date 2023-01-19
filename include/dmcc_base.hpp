/*
 * @Author: Wei Luo
 * @Date: 2022-12-12 17:55:51
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-19 14:43:46
 * @Note: Note
 */

#ifndef __DMCC_BASE_HPP__
#define __DMCC_BASE_HPP__

#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <math.h>


namespace ca = casadi;

class DMCCBase {
public:
  double g_acceleration_; // acceleration
  double arm_length_ = 0.0;
  double nominal_arm_distance_;
  double mass_quadrotor_;
  double mass_manipulator_;
  std::vector<double> I_quadrotor_;
  std::vector<double> I_arm_;
  ca::DM montage_offset_b_ = ca::DM(3, 1);
  double frame_size_;

  int num_controls_;
  int num_dofs_; // same as size of q
  int num_states_; // number of states,

  ca::Function derivative_lagrangian_function_;
  ca::Function lagrangian_function_;
  ca::Function external_force_function_;

  ca::Function get_external_force_function() { return external_force_function_; }

  ca::Function get_lagrangian_function() { return lagrangian_function_; }

  ca::Function get_derivative_lagrangian_function() {
    return derivative_lagrangian_function_;
  }

  template <typename T> T rotation_matrix(const T angle) {
    T phi = angle(0);
    T theta = angle(1);
    T psi = angle(2);

    T Rz = T(3, 3);
    Rz(0, 0) = T::cos(psi);
    Rz(0, 1) = -T::sin(psi);
    Rz(0, 2) = 0.0;
    Rz(1, 0) = T::sin(psi);
    Rz(1, 1) = T::cos(psi);
    Rz(1, 2) = 0.0;
    Rz(2, 0) = 0.0;
    Rz(2, 1) = 0.0;
    Rz(2, 2) = 1.0;

    T Ry = T(3, 3);
    Ry(0, 0) = T::cos(theta);
    Ry(0, 1) = 0.0;
    Ry(0, 2) = T::sin(theta);
    Ry(1, 0) = 0.0;
    Ry(1, 1) = 1.0;
    Ry(1, 2) = 0.0;
    Ry(2, 0) = -T::sin(theta);
    Ry(2, 1) = 0.0;
    Ry(2, 2) = T::cos(theta);

    T Rx = T(3, 3);
    Rx(0, 0) = 1.0;
    Rx(0, 1) = 0.0;
    Rx(0, 2) = 0.0;
    Rx(1, 0) = 0.0;
    Rx(1, 1) = T::cos(phi);
    Rx(1, 2) = -T::sin(phi);
    Rx(2, 0) = 0.0;
    Rx(2, 1) = T::sin(phi);
    Rx(2, 2) = T::cos(phi);

    return T::mtimes({Rz, Ry, Rx});
  }

  template <typename T> T average_rpy(T rpy_1, T rpy_2) {
    assert((rpy_1.size(0) == 3 and rpy_2.size(0) == 3));

    T rpy_sum = rpy_1 + rpy_2;
    for (int i = 0; i < 3; ++i) {
      T s_ = T::sin(rpy_1(i)) + T::sin(rpy_2(i));
      T c_ = T::cos(rpy_1(i)) + T::cos(rpy_2(i));
      rpy_sum(i) = T::atan2(s_, c_);
    }
    return rpy_sum;
  }

  template <typename T> T difference_rpy(T rpy_1, T rpy_2, T dt) {
    // assert((rpy_1.size(0) == 3 and rpy_2.size(0) == 3));

    T r_1 = rotation_matrix(rpy_1);
    T r_2 = rotation_matrix(rpy_2);

    auto unit_matrix = T(3, 3);
    unit_matrix(0, 0) = 1.0;
    unit_matrix(1, 1) = 1.0;
    unit_matrix(2, 2) = 1.0;

    return T::inv_skew(T::mtimes({r_1.T(), r_2}) - unit_matrix) / dt;
  }

  template <typename T> T difference_rpy(T rpy_1, T rpy_2, double dt) {
    T r_1 = rotation_matrix(rpy_1);
    T r_2 = rotation_matrix(rpy_2);

    ca::DM unit_matrix = ca::DM::zeros(3, 3);
    unit_matrix(0, 0) = 1.0;
    unit_matrix(1, 1) = 1.0;
    unit_matrix(2, 2) = 1.0;

    return T::inv_skew(T::mtimes({r_1.T(), r_2}) - unit_matrix) / dt;
  }

  template <typename T>
  T discrete_lagrange(T dt, T q_n, T q_np1, ca::Function fct_L) {
    T q = (q_n + q_np1) / 2.0;
    q(slice_rpy) = average_rpy(q_n(slice_rpy), q_np1(slice_rpy));
    T q_dot = (q_np1 - q_n) / dt;
    q_dot(slice_rpy) =
        difference_rpy((T)q_n(slice_rpy), (T)q_np1(slice_rpy), dt);

    if (num_dofs_ > 6) {
      for (int i = 6; i < num_dofs_; ++i) {
        T s_ = T::sin(q_n(i)) + T::sin(q_np1(i));
        T c_ = T::cos(q_n(i)) + T::cos(q_np1(i));
        q(i) = T::atan2(s_, c_);
        auto diff_ = q_np1(i) - q_n(i);
        q_dot(i) = T::atan2(T::sin(diff_), T::cos(diff_)) / dt;
      }
    }

    std::vector<T> input(2);
    input[0] = q;
    input[1] = q_dot;
    return dt + fct_L(input).at(0);
  }

  template <typename T>
  T discrete_lagrange_verlet(T dt, T q_n, T q_np1, ca::Function fct_L) {
    T q_dot = (q_np1 - q_n) / dt;
    q_dot(slice_rpy) =
        difference_rpy((T)q_n(slice_rpy), (T)q_np1(slice_rpy), dt);

    if (num_dofs_ > 6) {
      for (int i = 6; i < num_dofs_; ++i) {
        auto diff_ = q_np1(i) - q_n(i);
        q_dot(i) = T::atan2(T::sin(diff_), T::cos(diff_)) / dt;
      }
    }
    std::vector<T> input_1(2);
    input_1[0] = q_n;
    input_1[1] = q_dot;
    std::vector<T> input_2(2);
    input_2[0] = q_np1;
    input_2[1] = q_dot;

    return 0.5 * dt * fct_L(input_1).at(0) + 0.5 * dt * fct_L(input_2).at(0);
  }

  template <typename T>
  T discrete_lagrange_verlet(double dt, T q_n, T q_np1, ca::Function fct_L) {
    T q_dot = (q_np1 - q_n) / dt;
    q_dot(slice_rpy) =
        difference_rpy((T)q_n(slice_rpy), (T)q_np1(slice_rpy), dt);

    if (num_dofs_ > 6) {
      for (int i = 6; i < num_dofs_; ++i) {
        auto diff_ = q_np1(i) - q_n(i);
        q_dot(i) = T::atan2(T::sin(diff_), T::cos(diff_)) / dt;
      }
    }
    std::vector<T> input_1(2);
    input_1[0] = q_n;
    input_1[1] = q_dot;
    std::vector<T> input_2(2);
    input_2[0] = q_np1;
    input_2[1] = q_dot;

    return 0.5 * dt * fct_L(input_1).at(0) + 0.5 * dt * fct_L(input_2).at(0);
  }

  template <typename T>
  T discrete_forces(T dt, ca::Function f, T q, T u_n, T u_np1) {
    std::vector<T> input_1(2);
    input_1[0] = q;
    input_1[1] = u_n;
    std::vector<T> input_2(2);
    input_2[0] = q;
    input_2[1] = u_np1;
    return 0.25 * dt * (f(input_1).at(0) + f(input_2).at(0));
  }

  template <typename T>
  T discrete_forces(double dt, ca::Function f, T q, T u_n, T u_np1) {
    std::vector<T> input_1(2);
    input_1[0] = q;
    input_1[1] = u_n;
    std::vector<T> input_2(2);
    input_2[0] = q;
    input_2[1] = u_np1;
    return 0.25 * dt * (f(input_1).at(0) + f(input_2).at(0));
  }

  template <typename T>
  T average_velocity(T dt, T state_1, T state_2, int num_dof = 6) {
    auto d_q = (state_2 - state_1);

    for (int i = 0; i < 3; i++) {
      d_q(i) = d_q(i) / dt;
    }

    if (num_dof > 3) {
      d_q(slice_rpy) = difference_rpy((T) state_1(ca::Slice(3, 6)),
                                      (T) state_2(ca::Slice(3, 6)), dt);
    //   if (num_dof == 7)
    //     d_q(6) = T::atan2(T::sin(d_q(6)), T::cos(d_q(6))) / dt;
    }

    return d_q;
  }

  template <typename T>
  T average_velocity(double dt, T state_1, T state_2, int num_dof = 6) {
    auto d_q = (state_2 - state_1);

    for (int i = 0; i < 3; i++) {
      d_q(i) = d_q(i) / dt;
    }

    if (num_dof > 3) {
      d_q(slice_rpy) = difference_rpy((T)state_1(ca::Slice(3, 6)),
                                      (T)state_2(ca::Slice(3, 6)), dt);
      //   if (num_dof == 7)
      //     d_q(6) = T::atan2(T::sin(d_q(6)), T::cos(d_q(6))) / dt;
    }

    return d_q;
  }

private:
  ca::Slice all;
  ca::Slice slice_rpy = ca::Slice(3, 6);
};

#endif
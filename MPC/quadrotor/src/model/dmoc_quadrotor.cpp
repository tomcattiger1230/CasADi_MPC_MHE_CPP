/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:34:11
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-16 23:51:02
 * @Note: Note
 */

#include "dmoc_quadrotor.hpp"

DMOCUAV::DMOCUAV(const double mass_quadrotor,
                 const std::vector<double> inertia_moment,
                 const double frame_size, const double motor_torque_const,
                 const double g)
    : DerivationUAV(g, mass_quadrotor, inertia_moment, frame_size, motor_torque_const) {}
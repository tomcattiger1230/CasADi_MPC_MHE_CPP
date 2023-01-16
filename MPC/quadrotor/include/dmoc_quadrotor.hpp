/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:33:44
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-16 23:47:38
 * @Note: Note
 */

#ifndef __DMOC_QUADROTOR__
#define __DMOC_QUADROTOR__
#include <chrono>
#include "derivation_uav.hpp"
#include <iostream>

class DMOCUAV : public DerivationUAV{
    public:
      DMOCUAV(const double mass_quadrotor,
             const std::vector<double> inertia_moment, const double frame_size,
             const double motor_torque_const, const double g);
};

#endif
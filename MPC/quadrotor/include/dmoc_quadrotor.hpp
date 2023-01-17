/*
 * @Author: Wei Luo
 * @Date: 2023-01-05 22:33:44
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-17 20:47:00
 * @Note: Note
 */

#ifndef __DMOC_QUADROTOR__
#define __DMOC_QUADROTOR__
#include "derivation_uav.hpp"
#include <chrono>
#include <iostream>

class DMOCUAV : public DerivationUAV {
public:
  DMOCUAV(const double dt, const double N, const double mass_quadrotor,
          const std::vector<double> inertia_moment, const double frame_size,
          const double motor_torque_const, const double g);

  void initialization_formulation();

private:
  int prediction_horizon_;
  double dt_;
  ca::Slice slice_state_ = ca::Slice(0, 6);
};

#endif
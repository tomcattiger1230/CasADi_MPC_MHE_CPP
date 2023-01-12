/*
 * @Author: Wei Luo
 * @Date: 2023-01-06 21:14:04
 * @LastEditors: Wei Luo
 * @LastEditTime: 2023-01-12 18:09:40
 * @Note: Note
 */


#include "mpc_quadrotor.hpp"

int main() {
  constexpr double quadrotor_mass = 1.659;
  constexpr double g_acceleration = 9.8066;
  double roll_tau = 0.257;
  double roll_gain = 0.75;
  double pitch_tau = 0.259;
  double pitch_gain = 0.78;
  auto mpc_quadrotor_handle =
      std::make_shared<MPCQuadrotor>(MPCQuadrotor(0.1, 20, g_acceleration, roll_tau, roll_gain, pitch_tau, pitch_gain));
  mpc_quadrotor_handle->initialization_formulation();
  return 0;
}
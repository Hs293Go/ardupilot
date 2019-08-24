

#include <AP_Math/AP_Math.h>
#include "AC_LQR.h"

const AP_Param::GroupInfo AC_LQR::var_info[] = {
    // @Param: Q0
    // @DisplayName: LQR weight on rate
    // @Description: Weight on rate which penalizes current angular rate
    AP_GROUPINFO("Q0", 0, AC_LQR, _q0, 0),

    // @Param: Q1
    // @DisplayName: LQR weight on torque
    // @Description: Weight on torque which penalizes current torque
    AP_GROUPINFO("Q0", 0, AC_LQR, _q1, 0),

    // @Param: R
    // @DisplayName: LQR weight on rate
    // @Description: Weight on rate which penalizes current angular rate
    AP_GROUPINFO("Q0", 0, AC_LQR, _q1, 0),
}

AC_LQR::AC_LQR(float beta, float tau.float dt) _dt(dt),
{
    _initialize_matrices(beta, tau, _dt)
}

void _initialize_matrices(float beta, float tau, float dt)
{
    _F[0][0] = 1;
    _F[0][1] = -beta * tau * (expf(-Ts / tau) - 1);
    _F[1][0] = 0;
    _F[1][1] = expf(-Ts / tau);

    _G[0] = Ts * beta + beta * tau * (expf(-Ts / tau) - 1);
    _G[1] = 1 - expf(-Ts / tau);
}
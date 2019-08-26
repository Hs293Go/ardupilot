

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

void AC_LQR::_getStateInput(float beta, float tau, float dt)
{
    _F[0][0] = 1;
    _F[0][1] = -beta * tau * (expf(-Ts / tau) - 1);
    _F[1][0] = 0;
    _F[1][1] = expf(-Ts / tau);

    _G[0] = Ts * beta + beta * tau * (expf(-Ts / tau) - 1);
    _G[1] = 1 - expf(-Ts / tau);
}

/*
	LQR covariance cycle.

	P = A'PA - A'PB (R+B'PB)^-1 B'PA + Q

	K = (R + B'PB)^-1 B'PA

	TODO: Should redo and optimize the math of this, after folding gains calculation into this.
	      Because of obvious substitutions.
*/
bool AC_LQR::covarianceStep()
{
    float nP[2][2];

    float B0B0 = powf(_G[0], 2);
    float B1B1 = powf(_G[1], 2);
    float B0B1 = B0 * B1;
    float A00A00 = powf(_F[0][0], 2);
    float A01A01 = powf(_F[0][1], 2);
    float A11A11 = powf(_F[1][1], 2);
    float P01P10 = P01 * P10;
    float P00P11 = P00 * P11;

    float div = (R + B0B0 * P00 + B1B1 * P11 + B0B1 * (P01 + P10));

    nP[0][0] = (Q00 * R + P00 * (A00A00 * R + B0B0 * Q00) +
                B1B1 * (P11 * Q00 + A00A00 * (P00P11 - P01P10)) +
                B0B1 * Q00 * (P01 + P10)) /
               div;

    float common = A01 * (P00 * R + B1B1 * (P00P11 - P01P10)) - A11 * B0B1 * (P00P11 + P01P10);

    nP[1][0] = (A00 * (A11 * P10 * R + common)) /
               div;

    nP[0][1] = (A00 * (A11 * P01 * R + common)) /
               div;

    nP[1][1] = (Q11 * R + A01A01 * P00 * R + B0B0 * P00 * Q11 +
                A11A11 * (P11 * R + B0B0 * (P00P11 - P01P10)) +
                B1B1 * (P11 * Q11 + A01A01 * P00P11 - A01A01 * P01P10) +
                A01 * A11 * P01 * R + A01 * A11 * P10 * R + B0B1 * P01 * Q11 +
                B0B1 * (P10 * Q11 - 2 * A01 * A11 * P00P11 + 2 * A01 * A11 * P01P10)) /
               div;

    P00 = nP[0][0];
    P01 = nP[0][1];
    P10 = nP[1][0];
    P11 = nP[1][1];

    bool solved = false;

    div = (R + B1 * B1 * P11 + B0 * (B0 * P00 + B1 * (P01 + P10)));
    float nK[2];

    nK[0] = (A00 * (B0 * P00 + B1 * P10)) / div;
    nK[1] = (A01 * (B0 * P00 + B1 * P10) + A11 * (B0 * P01 + B1 * P11)) / div;

    if (fabsf(K0 - nK[0]) <= SOLVER_LQR_RATE_EPSILON && fabsf(K1 - nK[1]) <= SOLVER_LQR_TORQUE_EPSILON)
        solved = true;

    K0 = nK[0];
    K1 = nK[1];

    return solved;
}
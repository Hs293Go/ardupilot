#pragma once

/// @file	AC_LQG.h
/// @brief	Generic PID algorithm, with EEPROM-backed storage of constants.

#include <stdlib.h>
#include <cmath>
#include <AP_Param/AP_Param.h>
#include <AP_Math/vectorN.h>
#include <AP_Math/AP_Math.h>
#include <AP_Logger/AP_Logger.h>
#include <AP_Common/AP_Common.h>

class AC_LQR
{

    // Constructor for LQG
    AC_LQR(float beta, float tau, float dT);

    // get accessors
    AP_Float &Q0() { return _q0; }
    AP_Float &Q1() { return _q1; }
    AP_Float &R() { return _r; }

    // set accessors
    void Q0(const float v) { _q0.set(v); }
    void Q1(const float v) { _q1.set(v); }
    void R(const float v) { _r.set(v); }

private:
#if MATH_CHECK_INDEXES
    typedef VectorN<ftype, 2> Vector2;
    typedef VectorN<ftype, 3> Vector3;
    typedef VectorN<VectorN<ftype, 2>, 2> Matrix2;
    typedef VectorN<VectorN<ftype, 3>, 3> Matrix3;
#else
    typedef ftype Vector2[2];
    typedef ftype Vector3[3];
    typedef ftype Matrix2[2][2];
    typedef ftype Matrix3[3][3];
#endif
protected:
    void _initialize_matrices(float beta, float tau, float dT);

    // parameters
    AP_Float _q0;
    AP_Float _q1;
    AP_Float _r;

    // matrices and internal variables

    Vector2 _G; // Input matrix/vector, also denoted B
    Matrix2 _Q; // Weight on error matrix
    Matrix2 _P; // Coefficient matrix for performance index
    Matrix2 _F; // State matrix, also denoted A
    float _R;   // Weight on control effort
}
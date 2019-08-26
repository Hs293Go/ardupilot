#pragma once

#include <AP_Common/AP_Common.h>
#include <AP_Param/AP_Param.h>
#include <AP_Math/AP_Math.h>
#include "AC_LQR.h"

class RTKF
{
public:
    // Constructor for RTKF
    RTKF(float beta, float tau, float dt, float R, float q1, float q2, float q3, float biaslim);

    virtual ~RTKF(){};

    void rtkf_prediction_step(float signal, float input);
    void rtkf_predict_axis(float signal, float input);

    // get accessors
    AP_Float &q1() { return _q1; }
    AP_Float &q2() { return _q2; }
    AP_Float &q3() { return _q3; }

    // parameter var table
    //static const struct AP_Param::GroupInfo var_info[];

private:
#if MATH_CHECK_INDEXES
    typedef VectorN<ftype, 3> Vector3;
    typedef VectorN<VectorN<ftype, 3>, 3> Matrix3;
#else
    typedef ftype Vector3[3];
    typedef ftype Matrix3[3][3];
#endif
protected:
    void _getStateInput(float beta, float tau, float dt);

    Vector3 _X;
    Vector3 _K;
    Vector3 _G;
    Matrix3 _Q;
    Matrix3 _P;
    Matrix3 _F;

    //Parameters
    AP_Float _q1;
    AP_Float _q2;
    AP_Float _q3;
    AP_Float _r;

    float _dt;
    Matrix3 _Q;
    float _R;
    float _biaslim;
};
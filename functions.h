#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "Einsum.h"
#include <math.h>

// diagonal matrix elements
double fd(double *params)
{
    return params[0] + params[1];
}

//non-diagonal matrix elements 1
double fnd1(double *params)
{
    return params[1] - params[0];
}

//non-diagonal matrix elements 2
double fnd2(double *params)
{
    return pow(params[1] - params[0], 2);
}

double One(double *params)
{
    return 1;
}

double Zero(double *params)
{
    return 0;
}

#endif
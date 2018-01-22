#include "stdafx.h"
#include "MathFunction.h"
#include <float.h>

InterpolativeMathFunction::InterpolativeMathFunction(const MathFunction function)
        :z_(function)
{
    double h = 0.001;
    auto& z = z_;

    dx_ = [=](double x, double y)
    {
        return (z(x + h, y) - z(x, y)) / h;
    };
    dy_ = [=](double x, double y)
    {
        return (z(x, y + h) - z(x, y)) / h;
    };

    auto& dx = dx_;
    dxy_ = [=](double x, double y)
    {
        return (dx(x, y + h) - dx(x, y)) / h;
    };
}

const MathFunction& InterpolativeMathFunction::Z()
const
{
    return z_;
}

const MathFunction& InterpolativeMathFunction::Dx()
const
{
    return dx_;
}

const MathFunction& InterpolativeMathFunction::Dy()
const
{
    return dy_;
}

const MathFunction& InterpolativeMathFunction::Dxy() const
{
    return dxy_;
}

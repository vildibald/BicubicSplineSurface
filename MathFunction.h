#pragma once
#include <functional>

typedef std::function<double(double, double)> MathFunction;

class InterpolativeMathFunction
{
    MathFunction z_;
    MathFunction dx_;
    MathFunction dy_;
    MathFunction dxy_;
public:
    InterpolativeMathFunction(const MathFunction function);

    const MathFunction& Z() const;

    const MathFunction& Dx() const;

    const MathFunction& Dy() const;

    const MathFunction& Dxy() const;
};

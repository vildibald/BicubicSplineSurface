#pragma once

#include "KnotVector.h"
#include "KnotMatrix.h"
#include "MathFunction.h"

class InterpolativeMathFunction;

class Spline
{
    KnotMatrix z_;
    KnotMatrix zt_;
    KnotMatrix dx_;
    KnotMatrix dy_;
    KnotMatrix dxy_;
    KnotVector x_;
    KnotVector y_;

    Spline();
public:
    static Spline EmptySpline();

    bool IsEmpty() const;

    Spline(KnotVector rowVector, KnotVector columnVector);

    size_t RowsCount() const
    {
        return x_.size();
    }

    size_t ColumnsCount() const
    {
        return y_.size();
    }


    KnotMatrix Z() const
    {
        return z_;
    }

    KnotMatrix Zt() const
    {
        return z_;
    }


    KnotMatrix Dx() const
    {
        return dx_;
    }

    KnotMatrix Dy() const
    {
        return dy_;
    }

    KnotMatrix Dxy() const
    {
        return dxy_;
    }

    const KnotVector& X() const
    {
        return x_;
    }

    const KnotVector& Y() const
    {
        return y_;
    }

    const double X(size_t i) const
    {
        return x_[i];
    }

    const double Y(size_t j) const
    {
        return y_[j];
    }

    double Z(const size_t i, const size_t j) const
    {
        return z_[i][j];
    }

    double Zt(const size_t i, const size_t j) const
    {
        return z_[j][i];
    }

    double Dx(const size_t i, const size_t j) const
    {
        return dx_[j][i];
    }

    double Dy(const size_t i, const size_t j) const
    {
        return dy_[i][j];
    }

    double Dxy(const size_t i, const size_t j) const
    {
        return dxy_[i][j];
    }

    void SetZ(const size_t i, const size_t j, const double value)
    {
        z_[i][j] = value;
        zt_[j][i] = value;
    }

    void SetDx(const size_t i, const size_t j, const double value)
    {
        dx_[j][i] = value;
    }

    void SetDy(const size_t i, const size_t j, const double value)
    {
        dy_[i][j] = value;
    }

    void SetDxy(const size_t i, const size_t j, const double value)
    {
        dxy_[i][j] = value;
    }

    void Initialize(InterpolativeMathFunction mathFunction);

    void Print();
};

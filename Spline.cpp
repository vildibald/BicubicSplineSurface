#include "stdafx.h"
#include "Spline.h"
#include "utils.h"
#include <iostream>
#include "MathFunction.h"

void Spline::Print()
{
    using namespace std;
    cout << "---------- Knot matrix ----------" << endl;
    for (size_t i = 0; i < RowsCount(); i++)
    {
        cout << "Row " << i << " :\n";
        for (size_t j = 0; j < ColumnsCount(); j++)
        {
            cout << j << ":\n"
                 << "z: " << z_[i][j] << '\n'
                 << "dx: " << dx_[j][i] << '\n'
                 << "dy: " << dy_[i][j] << '\n'
                 << "dxy: " << dxy_[i][j] << '\n';
        }
        cout << endl;
    }
    cout << "-------------------------------" << endl;
}


Spline Spline::EmptySpline()
{
    Spline nullval;
    return nullval;
}

bool Spline::IsEmpty() const
{
    if (RowsCount() < 1 || ColumnsCount() < 1)
        return true;
    return false;
}

Spline::Spline()
        : z_(), dx_(), dy_(), dxy_(), x_(), y_()
{
}

Spline::Spline(KnotVector rowVector, KnotVector columnVector)
        : z_(), dx_(), dy_(), dxy_(), x_(std::move(rowVector)), y_(std::move(columnVector))
{
    z_.resize(x_.size());
    dx_.resize(x_.size());
    dy_.resize(x_.size());
    dxy_.resize(x_.size());

    for (int i = 0; i < x_.size(); ++i) {
        z_[i].resize(y_.size());
        dx_[i].resize(y_.size());
        dy_[i].resize(y_.size());
        dxy_[i].resize(y_.size());
    }
}

void Spline::Initialize(InterpolativeMathFunction mathFunction) {

    // Init Z
    for (size_t i = 0; i < RowsCount(); i++) {
        for (size_t j = 0; j < ColumnsCount(); j++) {
            const auto z = mathFunction.Z()(X(i), Y(j));
            //mathFunction.Z(u,v); //Z(u, v);
            SetZ(i, j, z);
        }

    }
    // Init Dx
    const auto uKnotCountMin1 = RowsCount() - 1;
    for (size_t j = 0; j < ColumnsCount(); j++) {
        SetDx(0, j, mathFunction.Dx()(X(0), Y(j))); //mathFunction.Dx(values(0,j).X, values(0,j).Y);
        SetDx(uKnotCountMin1, j, mathFunction.Dx()(X(uKnotCountMin1),
                                                   Y(j)));
    }
    // Init Dy
    const auto vKnotCountMin1 = ColumnsCount() - 1;
    for (size_t i = 0; i < RowsCount(); i++) {
        SetDy(i, 0, mathFunction.Dy()(X(i), Y(0)));
        SetDy(i, vKnotCountMin1,
              mathFunction.Dy()(X(i), Y(vKnotCountMin1))
        );
    }
    // Init Dxy
    SetDxy(0, 0, mathFunction.Dxy()(X(0), Y(0)));
    SetDxy(uKnotCountMin1, 0, mathFunction.Dxy()(X(uKnotCountMin1),
                                                 Y(0)));
    SetDxy(0, vKnotCountMin1, mathFunction.Dxy()(X(0),
                                                 Y(vKnotCountMin1)));
    SetDxy(uKnotCountMin1, vKnotCountMin1, mathFunction.Dxy()(X(uKnotCountMin1),
                                                              Y(vKnotCountMin1)));
}
#include "stdafx.h"
#include "FullAlgorithm.h"
#include "utils.h"
#include <array>

FullAlgorithm::FullAlgorithm(const InterpolativeMathFunction f)
        : function_(f), xTridiagonals_(),
          yTridiagonals_(), isParallel_(false)
{
}

Spline FullAlgorithm::Calculate(const KnotVector xVector,
                                const KnotVector yVector)
{
    Spline spline{std::move(xVector), std::move(yVector)};
    Initialize(spline);
    timer_.Reset();
    timer_.Start();
    FillDx(spline);
    FillDy(spline);
    FillDxy(spline);
    FillDyx(spline);
    timer_.Stop();

    return spline;
}

void FullAlgorithm::FillDx(Spline& spline)
{
    auto differenceGetter = [&](size_t i)
    {
        return spline.X(i);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Z(i, j);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dx(i, j);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDx(i, j, value);
    };
    FillD(spline.ColumnsCount(), spline.RowsCount(), xTridiagonals_, differenceGetter, parameterGetter,
          derivationGetter, derivationSetter);
}

void FullAlgorithm::FillDy(Spline& spline)
{
    auto differenceGetter = [&](size_t i)
    {
        return spline.Y(i);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Z(j, i);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dy(j, i);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDy(j, i, value);
    };
    FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, differenceGetter, parameterGetter,
          derivationGetter, derivationSetter);
}

void FullAlgorithm::FillDxy(Spline& spline)
{
    auto differenceGetter = [&](size_t i)
    {
        return spline.X(i);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Dy(i, j);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dxy(i, j);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDxy(i, j, value);
    };

    size_t loop[] = {0, spline.ColumnsCount() - 1};
    for (auto j : loop)
    {
        Solve(spline.ColumnsCount(), spline.RowsCount(), xTridiagonals_, differenceGetter,
              parameterGetter, derivationGetter, derivationSetter, j);
    }
}

void FullAlgorithm::FillDyx(Spline& spline)
{
    auto differenceGetter = [&](size_t i)
    {
        return spline.Y(i);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Dx(j, i);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dxy(j, i);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDxy(j, i, value);
    };
    FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, differenceGetter, parameterGetter,
          derivationGetter, derivationSetter);
}

void FullAlgorithm::Initialize(Spline& spline)
{
    InitializeTridiagonals(spline);
    spline.Initialize(function_);
}

bool FullAlgorithm::IsParallel() const
{
    return isParallel_;
}

CalculationMode FullAlgorithm::GetCalculationMode() const
{
    return calculationMode_;
}

void FullAlgorithm::SetCalculationMode(const CalculationMode calculationMode)
{
    calculationMode_ = calculationMode;
}

void FullAlgorithm::InParallel(const bool value)
{
    isParallel_ = value;
}

void FullAlgorithm::Parallelize(const bool inParallel)
{
    xTridiagonals_.Parallelize(inParallel);
    yTridiagonals_.Parallelize(inParallel);
}

void FullAlgorithm::InitializeTridiagonals(Spline& spline)
{
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().clear();

    switch (calculationMode_)
    {
        case NON_OPTIMIZED:
            xTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateFullWithDivisionsTridiagonal(spline.X(),
                                                                             spline.RowsCount()));
            yTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateFullWithDivisionsTridiagonal(spline.Y(),
                                                                             spline.ColumnsCount()));
            break;
        case OPTIMIZED_DIVISIONS:
        case OPTIMIZED_DIVISIONS_BUFFERED:
        default:
            xTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateFullTridiagonal(spline.X(), spline.RowsCount()));
            yTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateFullTridiagonal(spline.Y(), spline.ColumnsCount()));
            break;
    }


    Parallelize(isParallel_);
}

double FullAlgorithm::ExecutionTime()
{
    return timer_.ExecutionTime();
}

double FullAlgorithm::AllTime()
{
    return timer_.AllTime() + xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].
            AllTime();
}

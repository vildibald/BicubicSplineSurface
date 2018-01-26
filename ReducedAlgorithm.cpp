#include "stdafx.h"
#include "utils.h"
#include "ReducedAlgorithm.h"

ReducedAlgorithm::ReducedAlgorithm(const InterpolativeMathFunction f)
        : function_(f), xTridiagonals_(),
          yTridiagonals_(), isParallel_(false)
{
}

Spline ReducedAlgorithm::Calculate(KnotVector xVector,
                                   KnotVector yVector)
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

void ReducedAlgorithm::FillDx(Spline& spline)
{
    auto differenceGetter = [&](size_t i)
    {
        return spline.X(i);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Zt(i, j);
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

void ReducedAlgorithm::FillDy(Spline& spline)
{
    auto differenceGetter = [&](size_t j)
    {
        return spline.Y(j);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Z(i, j);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dy(i, j);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDy(i, j, value);
    };
    FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, differenceGetter, parameterGetter,
          derivationGetter, derivationSetter);
}

void ReducedAlgorithm::FillDxy(Spline& spline)
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

void ReducedAlgorithm::FillDyx(Spline& spline)
{
    auto differenceGetter = [&](size_t j)
    {
        return spline.Y(j);
    };
    auto parameterGetter = [&](size_t i, size_t j)
    {
        return spline.Dx(i, j);
    };
    auto derivationGetter = [&](size_t i, size_t j)
    {
        return spline.Dxy(i, j);
    };
    auto derivationSetter = [&](size_t i, size_t j, double value)
    {
        spline.SetDxy(i, j, value);
    };
    FillD(spline.RowsCount(), spline.ColumnsCount(), yTridiagonals_, differenceGetter, parameterGetter,
          derivationGetter, derivationSetter);
}

void ReducedAlgorithm::Initialize(Spline& spline)
{
    InitializeTridiagonals(spline);
    spline.Initialize(function_);
}

bool ReducedAlgorithm::IsParallel() const
{
    return isParallel_;
}

CalculationMode ReducedAlgorithm::GetCalculationMode() const
{
    return calculationMode_;
}

void ReducedAlgorithm::SetCalculationMode(const CalculationMode calculationMode)
{
    calculationMode_ = calculationMode;
}

void ReducedAlgorithm::InParallel(const bool value)
{
    isParallel_ = value;
}

void ReducedAlgorithm::Parallelize(const bool inParallel)
{
    xTridiagonals_.Parallelize(inParallel);
    yTridiagonals_.Parallelize(inParallel);
}

void ReducedAlgorithm::InitializeTridiagonals(Spline& spline)
{
    xTridiagonals_.GetAll().clear();
    xTridiagonals_.GetAll().clear();
    switch (calculationMode_)
    {
        case NON_OPTIMIZED:
            xTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateReducedWithDivisionsTridiagonal(spline.X(),
                                                                                spline.RowsCount()));
            yTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateReducedWithDivisionsTridiagonal(spline.Y(),
                                                                                spline.ColumnsCount()));
            break;
        case OPTIMIZED_DIVISIONS:
        case OPTIMIZED_DIVISIONS_BUFFERED:
        default:
            xTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateReducedTridiagonal(spline.X(), spline.RowsCount()));
            yTridiagonals_.GetAll().emplace_back(
                    Tridiagonal::Factory::CreateReducedTridiagonal(spline.Y(), spline.ColumnsCount()));
            break;
    }

    Parallelize(isParallel_);
}

double ReducedAlgorithm::ExecutionTime()
{
    return timer_.ExecutionTime();
}

double ReducedAlgorithm::AllTime()
{
    return timer_.AllTime() + xTridiagonals_.GetAll()[0].AllTime() + yTridiagonals_.GetAll()[0].
            AllTime();
}

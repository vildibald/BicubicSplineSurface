//
// Created by Viliam on 15.2.2017.
//
#include "MathFunction.h"
#include "Spline.h"
#include <functional>
#include <vector>
#include "Tridiagonal.h"
#include <omp.h>
#include "Timer.h"
#include "utils.h"
#include "CalculationMode.h"


class FullAlgorithm
{
    InterpolativeMathFunction function_;
    Tridiagonals xTridiagonals_;
    Tridiagonals yTridiagonals_;
    Timer timer_;
    bool isParallel_;
    CalculationMode calculationMode_;

public:
    FullAlgorithm(const InterpolativeMathFunction f);

    Spline Calculate(const KnotVector xVector,
                     const KnotVector yVector);

    double ExecutionTime();

    double AllTime();

    void InParallel(bool value);

    bool IsParallel() const;

    CalculationMode GetCalculationMode() const;

    void SetCalculationMode(const CalculationMode calculationMode);

private:
    void Initialize(Spline& values);

    void FillDx(Spline& values);

    void FillDy(Spline& values);

    void FillDxy(Spline& values);

    void FillDyx(Spline& values);

    void Parallelize(bool inParallel);

    void InitializeTridiagonals(Spline& spline);

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void FillD(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
               const DifferenceGetter h, const ParameterGetter& p, const DerivationGetter& dget,
               DerivationSetter& dset)
    {
        utils::For(0, static_cast<int>(systemCount), 1, isParallel_, [&](int j)
        {
            Solve(systemCount, derivationCount, tridiagonals, h, p, dget, dset, j);
        });
    }

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void Solve(const size_t systemCount, const size_t derivationCount, Tridiagonals& tridiagonals,
               const DifferenceGetter h, const ParameterGetter& p, const DerivationGetter& dget,
               DerivationSetter& dset, size_t systemIdx)
    {
        switch (calculationMode_)
        {
            case NON_OPTIMIZED:
                SolveNonOptimized(systemCount, derivationCount, tridiagonals, h, p, dget, dset, systemIdx);
                break;
            case OPTIMIZED_DIVISIONS:
                SolveOptimized(systemCount, derivationCount, tridiagonals, h, p, dget, dset, systemIdx);
                break;
            case OPTIMIZED_DIVISIONS_BUFFERED:
            default:
                SolveBuffered(systemCount, derivationCount, tridiagonals, h, p, dget, dset, systemIdx);
                break;
        }
    }

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void SolveBuffered(const size_t systemCount, const size_t derivationCount,
                       Tridiagonals& tridiagonals,
                       const DifferenceGetter h, const ParameterGetter& p,
                       const DerivationGetter& dget,
                       DerivationSetter& dset, size_t systemIdx)
    {
        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto& r2 = tridiagonal.RhsCoeficients()[2];
        const auto& r1 = tridiagonal.RhsCoeficients()[1];
        const auto& r0 = tridiagonal.RhsCoeficients()[0];
        const auto l2 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < rightSide.size(); ++i)
        {
            rightSide[i] = r2[i] * p(i + 2, systemIdx) + r1[i] * p(i + 1, systemIdx) + r0[i] * p(
                    i, systemIdx);
        }

        rightSide.front() -= l0 * dget(0, systemIdx);;
        rightSide.back() -= l2 * dget(derivationCount - 1, systemIdx);
        tridiagonal.Solve();

        for (auto i = 0; i < rightSide.size(); ++i)
        {
            dset(i + 1, systemIdx, rightSide[i]);
        }
    }

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void SolveOptimized(const size_t systemCount, const size_t derivationCount,
                        Tridiagonals& tridiagonals,
                        const DifferenceGetter h, const ParameterGetter& p,
                        const DerivationGetter& dget,
                        DerivationSetter& dset, size_t systemIdx)
    {
        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto l2 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < rightSide.size(); ++i)
        {
            auto h0 = h(i + 1) - h(i);
            auto h1 = h(i + 2) - h(i + 1);
            auto h00 = h0 * h0;
            auto h11 = h1 * h1;
            rightSide[i] = 3 * (h00 * p(i + 2, systemIdx) + (h11 - h00) * p(i + 1, systemIdx) - h11 *
                                                                                                p(i, systemIdx));
        }

        rightSide.front() -= l0 * dget(0, systemIdx);;
        rightSide.back() -= l2 * dget(derivationCount - 1, systemIdx);
        tridiagonal.Solve();

        for (auto i = 0; i < rightSide.size(); ++i)
        {
            dset(i + 1, systemIdx, rightSide[i]);
        }
    }

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void SolveNonOptimized(const size_t systemCount, const size_t derivationCount,
                           Tridiagonals& tridiagonals,
                           const DifferenceGetter h, const ParameterGetter& p,
                           const DerivationGetter& dget,
                           DerivationSetter& dset, size_t systemIdx)
    {
        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto l2 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < rightSide.size(); ++i)
        {
            auto h0 = h(i + 1) - h(i);
            auto h1 = h(i + 2) - h(i + 1);
            auto h00 = h0 * h0;
            auto h11 = h1 * h1;
            auto h10 = h1 * h0;
            rightSide[i] = 3 * (h0 / h1 * (p(i + 2, systemIdx) - p(i + 1, systemIdx))
                                + h1 / h0 * (p(i + 1, systemIdx)) - p(i, systemIdx));
        }

        rightSide.front() -= l0 * dget(0, systemIdx);;
        rightSide.back() -= l2 * dget(derivationCount - 1, systemIdx);
        tridiagonal.Solve();

        for (auto i = 0; i < rightSide.size(); ++i)
        {
            dset(i + 1, systemIdx, rightSide[i]);
        }
    }
};

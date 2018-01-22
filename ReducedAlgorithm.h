#pragma once
#include "MathFunction.h"
#include "Spline.h"
#include <vector>
#include "Tridiagonal.h"
#include "CalculationMode.h"

class ReducedAlgorithm
{
    InterpolativeMathFunction function_;
    Tridiagonals xTridiagonals_;
    Tridiagonals yTridiagonals_;
    Timer timer_;
    bool isParallel_;
    CalculationMode calculationMode_ = OPTIMIZED_DIVISIONS_BUFFERED;

public:
    ReducedAlgorithm(const InterpolativeMathFunction f);

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
        const auto unknownsCount = derivationCount;
        const auto even = unknownsCount % 2 == 0;
        const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;

        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto& r4 = tridiagonal.RhsCoeficients()[4];
        const auto& r3 = tridiagonal.RhsCoeficients()[3];
        const auto& r2 = tridiagonal.RhsCoeficients()[2];
        const auto& r1 = tridiagonal.RhsCoeficients()[1];
        const auto& r0 = tridiagonal.RhsCoeficients()[0];
        const auto l4 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < equationsCount - 1; ++i)
        {
            auto i21 = 2 * (i + 1);
            rightSide[i] = r4[i] * p(i21 + 2, systemIdx) + r3[i] * p(i21 + 1, systemIdx) + r2[i] *
                                                                                           p(i21, systemIdx) + r1[i] *
                                                                                                               p(i21 - 1, systemIdx) + r0[i] * p(i21 - 2, systemIdx);
        }
        rightSide.front() -= l4 * dget(0, systemIdx);
        {
            const size_t i = equationsCount - 1;
            const auto i21 = 2 * (i + 1);
            if (even)
            {
                // TODO: Needs support for even count.
                rightSide.back() = 0;
            }
            else
            {
                rightSide.back() = r4[i] * p(i21 + 2, systemIdx) + r3[i] * p(i21 + 1, systemIdx) + r2[i] *
                                                                                                   p(i21, systemIdx) + r1[i] *
                                                                                                                       p(i21 - 1, systemIdx) + r0[i] * p(i21 - 2, systemIdx);
            }
            rightSide.back() -= l0 * dget(derivationCount - 1, systemIdx);
        }
        tridiagonal.Solve();

        const auto& rst2 = tridiagonal.RhsCoeficients()[7];
        const auto& rst1 = tridiagonal.RhsCoeficients()[6];
        const auto& rst0 = tridiagonal.RhsCoeficients()[5];
        const auto& rstd2 = tridiagonal.RhsCoeficients()[9];
        const auto& rstd0 = tridiagonal.RhsCoeficients()[8];

        for (size_t i = 0; i < equationsCount; ++i)
        {
            const auto evenI = 2 * (i + 1);
            dset(evenI, systemIdx, rightSide[i]);

            const auto oddI = 2 * i + 1;
            dset(oddI, systemIdx,
                 rst2[i] * p(oddI + 1, systemIdx) + rst1[i] * p(oddI, systemIdx) + rst0[i] * p(
                         oddI - 1, systemIdx) +
                 rstd2[i] * dget(oddI + 1, systemIdx) + rstd0[i] * dget(oddI - 1, systemIdx)
            );
        }
        {
            const size_t i = equationsCount;
            const auto oddI = 2 * i + 1;
            dset(oddI, systemIdx,
                 rst2[i] * p(oddI + 1, systemIdx) + rst1[i] * p(oddI, systemIdx) + rst0[i] * p(
                         oddI - 1, systemIdx) +
                 rstd2[i] * dget(oddI + 1, systemIdx) + rstd0[i] * dget(oddI - 1, systemIdx)
            );
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
        const auto unknownsCount = derivationCount;
        const auto even = unknownsCount % 2 == 0;
        const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto l4 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < equationsCount - 1; ++i)
        {
            const auto i2 = 2 * (i + 1);
            const auto i1 = i2 - 1;
            const auto i0 = i2 - 2;
            const auto i3 = i2 + 1;
            const auto i4 = i2 + 2;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto h2 = h(i3) - h(i2);
            const auto h3 = h(i4) - h(i3);
            const auto h11 = h1 * h1;
            const auto h22 = h2 * h2;
            const auto h3Plush2 = h3 + h2;
            const auto h1Plush0 = h1 + h0;


            rightSide[i] = 3 * (
                    h3 * (h3Plush2 * h22 * (h1Plush0 * h1Plush0 * p(i1, systemIdx) - h11 * p(i0, systemIdx)) + h0 *
                                                                                                               (h3Plush2 * h0 * h22 -
                                                                                                                h1Plush0 * (2 * h3Plush2 * (h22 - h11) + h3 * h11)) * p(i2, systemIdx))
                    - h1Plush0 * h11 * h0 * (h3Plush2 * p(i3, systemIdx) - h22 * p(i4, systemIdx))
            );
        }
        rightSide.front() -= l4 * dget(0, systemIdx);
        {
            const int i = equationsCount - 1;
            const auto i2 = 2 * (i + 1);
            const auto i1 = i2 - 1;
            const auto i0 = i2 - 2;
            const auto i3 = i2 + 1;
            const auto i4 = i2 + 2;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto h2 = h(i3) - h(i2);
            const auto h3 = h(i4) - h(i3);
            const auto h11 = h1 * h1;
            const auto h22 = h2 * h2;
            const auto h3Plush2 = h3 + h2;
            const auto h1Plush0 = h1 + h0;
            if (even)
            {
                // TODO: Needs support for even count.
                rightSide.back() = 0;
            }
            else
            {
                rightSide.back() = 3 * (
                        h3 * (h3Plush2 * h22 * (h1Plush0 * h1Plush0 * p(i1, systemIdx) - h11 * p(i0, systemIdx)) + h0 *
                                                                                                                   (h3Plush2 * h0 * h22 -
                                                                                                                    h1Plush0 * (2 * h3Plush2 * (h22 - h11) + h3 * h11)) * p(i2, systemIdx))
                        - h1Plush0 * h11 * h0 * (h3Plush2 * p(i3, systemIdx) - h22 * p(i4, systemIdx))
                );
            }
            rightSide.back() -= l0 * dget(derivationCount - 1, systemIdx);
        }
        tridiagonal.Solve();

        for (size_t i = 0; i < equationsCount; ++i)
        {
            const auto evenI = 2 * (i + 1);
            dset(evenI, systemIdx, rightSide[i]);


            const auto i1 = 2 * i + 1;
            const auto i2 = i1 + 1;
            const auto i0 = i1 - 1;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto threeDivh0 = 3 / h0;
            const auto threeDivh1 = 3 / h1;
            const auto h1Plush0 = h1 + h0;
            const auto rth1Plush0h10 = 1 / (2 * h1Plush0 * h1 * h0);

            dset(i1, systemIdx,
                 rth1Plush0h10 * (
                         h0 * (threeDivh1 * (p(i2, systemIdx) - p(i1, systemIdx)) - dget(i2, systemIdx)) +
                         h1 * (threeDivh0 * (p(i1, systemIdx) - p(i0, systemIdx)) - dget(i0, systemIdx))
                 )
            );
        }
        const auto i = equationsCount;
        const auto i1 = 2 * i + 1;
        const auto i2 = i1 + 1;
        const auto i0 = i1 - 1;

        const auto h0 = h(i1) - h(i0);
        const auto h1 = h(i2) - h(i1);
        const auto threeDivh0 = 3 / h0;
        const auto threeDivh1 = 3 / h1;
        const auto h1Plush0 = h1 + h0;
        const auto rth1Plush0h10 = 1 / (2 * h1Plush0 * h1 * h0);

        dset(i1, systemIdx,
             rth1Plush0h10 * (
                     h0 * (threeDivh1 * (p(i2, systemIdx) - p(i1, systemIdx)) - dget(i2, systemIdx)) +
                     h1 * (threeDivh0 * (p(i1, systemIdx) - p(i0, systemIdx)) - dget(i0, systemIdx))
             )
        );
    }

    template <typename DifferenceGetter, typename ParameterGetter, typename DerivationGetter, typename
    DerivationSetter>
    void SolveNonOptimized(const size_t systemCount, const size_t derivationCount,
                           Tridiagonals& tridiagonals,
                           const DifferenceGetter h, const ParameterGetter& p,
                           const DerivationGetter& dget,
                           DerivationSetter& dset, size_t systemIdx)
    {
        const auto unknownsCount = derivationCount;
        const auto even = unknownsCount % 2 == 0;
        const auto equationsCount = even ? unknownsCount / 2 - 2 : unknownsCount / 2 - 1;
        auto& tridiagonal = tridiagonals.Get();
        auto& rightSide = tridiagonal.RightSideBuffer();
        const auto l4 = tridiagonal.Lhs2Coeficients().front();
        const auto l0 = tridiagonal.Lhs0Coeficients().back();

        for (size_t i = 0; i < equationsCount - 1; ++i)
        {
            const auto i2 = 2 * (i + 1);
            const auto i1 = i2 - 1;
            const auto i0 = i1 - 1;
            const auto i3 = i2 + 1;
            const auto i4 = i3 + 1;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto h2 = h(i3) - h(i2);
            const auto h3 = h(i4) - h(i3);
            const auto h11 = h1 * h1;
            const auto h22 = h2 * h2;
            const auto h311 = h3 * h11;
            const auto h3Plush2 = h3 + h2;
            const auto h1Plush0 = h1 + h0;
            const auto h3Plush2Pow2 = h3Plush2 * h3Plush2;
            const auto h1Plush0Pow2 = h1Plush0 * h1Plush0;

            rightSide[i] = 3 * (
                    h1Plush0 / h2 * (
                            h1 * (h22 * p(i4, systemIdx) - h3Plush2Pow2 * p(i3, systemIdx)) / h3 -
                            (2 * h3Plush2 * (h22 - h11) + h311) / h1 * p(i2, systemIdx)
                    ) + h3Plush2 * h2 / h1 * (
                            (h1Plush0Pow2 * p(i1, systemIdx) - h11 * p(i0, systemIdx)) / h0 + h0 * p(i2, systemIdx)
                    )
            );
        }
        rightSide.front() -= l4 * dget(0, systemIdx);
        {
            const size_t i = equationsCount - 1;
            const auto i2 = 2 * (i + 1);
            const auto i1 = i2 - 1;
            const auto i0 = i1 - 1;
            const auto i3 = i2 + 1;
            const auto i4 = i3 + 1;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto h2 = h(i3) - h(i2);
            const auto h3 = h(i4) - h(i3);
            const auto h11 = h1 * h1;
            const auto h22 = h2 * h2;
            const auto h311 = h3 * h11;
            const auto h3Plush2 = h3 + h2;
            const auto h1Plush0 = h1 + h0;
            const auto h3Plush2Pow2 = h3Plush2 * h3Plush2;
            const auto h1Plush0Pow2 = h1Plush0 * h1Plush0;
            if (even)
            {
                // TODO: Needs support for even count.
                rightSide.back() = 0;
            }
            else
            {
                rightSide.back() = 3 * (
                        h1Plush0 / h2 * (
                                h1 * (h22 * p(i4, systemIdx) - h3Plush2Pow2 * p(i3, systemIdx)) / h3 -
                                (2 * h3Plush2 * (h22 - h11) + h311) / h1 * p(i2, systemIdx)
                        ) + h3Plush2 * h2 / h1 * (
                                (h1Plush0Pow2 * p(i1, systemIdx) - h11 * p(i0, systemIdx)) / h0 + h0 * p(i2, systemIdx)
                        )
                );
            }
            rightSide.back() -= l0 * dget(derivationCount - 1, systemIdx);
        }
        tridiagonal.Solve();

        for (size_t i=0; i < equationsCount; ++i)
        {
            const auto evenI = 2 * (i + 1);
            dset(evenI, systemIdx, rightSide[i]);


            const auto i1 = 2 * i + 1;
            const auto i2 = i1 + 1;
            const auto i0 = i1 - 1;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto threeDivh0 = 3 / h0;
            const auto threeDivh1 = 3 / h1;
            const auto h1Plush0 = h1 + h0;
            const auto rth1Plush0h10 = 1 / (2 * h1Plush0);

            dset(i1, systemIdx,
                 rth1Plush0h10 * (h0 * (threeDivh1 * (p(i2, systemIdx) - p(i1, systemIdx) - dget(
                         i2, systemIdx))) +
                                  h1 * (threeDivh0 * (p(i1, systemIdx) - p(i0, systemIdx) - dget(
                                          i0, systemIdx)))
                 )
            );
        }
        {
            const size_t i = equationsCount;
            const auto i1 = 2 * i + 1;
            const auto i2 = i1 + 1;
            const auto i0 = i1 - 1;

            const auto h0 = h(i1) - h(i0);
            const auto h1 = h(i2) - h(i1);
            const auto threeDivh0 = 3 / h0;
            const auto threeDivh1 = 3 / h1;
            const auto h1Plush0 = h1 + h0;
            const auto rth1Plush0h10 = 1 / (2 * h1Plush0);

            dset(i1, systemIdx,
                 rth1Plush0h10 * (h0 * (threeDivh1 * (p(i2, systemIdx) - p(i1, systemIdx) - dget(
                         i2, systemIdx))) +
                                  h1 * (threeDivh0 * (p(i1, systemIdx) - p(i0, systemIdx) - dget(
                                          i0, systemIdx)))
                 )
            );
        }
    }
};

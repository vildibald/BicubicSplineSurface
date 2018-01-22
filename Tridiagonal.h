#pragma once

#include "stdafx.h"
#include <vector>
#include "Spline.h"
#include "KnotVector.h"
#include "MultithreadPreparator.h"
#include <numeric>
#include "Timer.h"


class Tridiagonal;

class Tridiagonals;

class Tridiagonals final
{
    std::vector<Tridiagonal> tridiagonals_;

public:
    std::vector<Tridiagonal>& GetAll();

    Tridiagonal& Get();

    void Parallelize(const bool inParallel)
    {
        MultithreadPreparator multithreadPreparator;
        multithreadPreparator.PrepareVector(inParallel, GetAll());
    }

    void Initialize(Tridiagonal tridiagonal);
};

class Tridiagonal final
{
    std::vector<double> luBuffer_;
    std::vector<double> rightSideBuffer_;
    std::vector<double> lhs0Coeficients_;
    std::vector<double> lhs1Coeficients_;
    std::vector<double> lhs2Coeficients_;
    std::vector<std::vector<double>> rhsCoeficients_;
    size_t numUnknowns_;
    size_t problemSize_;
    Timer timer_;

    Tridiagonal(std::vector<double> lhs0Coeficients,
                std::vector<double> lhs1Coeficients,
                std::vector<double> lhs2Coeficients,
                std::vector<std::vector<double>> rhsCoeficients,
                size_t numUnknowns,
                size_t problemSize
    );

public:

    std::vector<double>& Solve();

    std::vector<double>& ResetBufferAndGet();

    std::vector<double>& Buffer();

    const std::vector<double>& Lhs0Coeficients() const;

    const std::vector<double>& Lhs1Coeficients() const;

    const std::vector<double>& Lhs2Coeficients() const;

    const std::vector<std::vector<double>>& RhsCoeficients() const;

    std::vector<double>& RightSideBuffer();

    size_t NumUnknowns() const;

    size_t ProblemSize() const;

    double ExecutionTime() {
        return timer_.ExecutionTime();
    }

    double AllTime() {
        return timer_.AllTime();
    }

    class Factory final
    {
    public:

        static Tridiagonal
        CreateEmptyTridiagonal();

        static Tridiagonal
        CreateFullTridiagonal(const KnotVector& knotVector, size_t numUnknowns);

        static Tridiagonal
        CreateReducedTridiagonal(const KnotVector& knotVector, size_t numUnknowns);

        static Tridiagonal
        CreateFullWithDivisionsTridiagonal(const KnotVector& knotVector, size_t numUnknowns);

        static Tridiagonal
        CreateReducedWithDivisionsTridiagonal(const KnotVector& knotVector, size_t numUnknowns);
    };

    static double AccumulateAllTimes(Tridiagonals& tridiagonals) {
        return std::accumulate(tridiagonals.GetAll().begin(), tridiagonals.GetAll().end(),
                               tridiagonals.GetAll()[0].AllTime(),
                               [](double time, Tridiagonal& tridiagonal) {
                                   return time + tridiagonal.AllTime();
                               });
    }

    static double AccumulateExecutionTimes(Tridiagonals& tridiagonals);
};

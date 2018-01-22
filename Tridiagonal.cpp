#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>

std::vector<Tridiagonal>& Tridiagonals::GetAll()
{
    return tridiagonals_;
}

Tridiagonal& Tridiagonals::Get()
{
    if (tridiagonals_.size() == 1)
    {
        return tridiagonals_[0];
    }
    return tridiagonals_[omp_get_thread_num()];
}

void Tridiagonals::Initialize(Tridiagonal tridiagonal)
{
    GetAll().clear();
    GetAll().emplace_back(std::move(tridiagonal));
}

std::vector<double>&
Tridiagonal::ResetBufferAndGet()
{
    auto& buffer = luBuffer_;
    std::fill(buffer.begin(), buffer.end(), 1);
    return buffer;
}

std::vector<double>&
Tridiagonal::Buffer()
{
    return luBuffer_;
}

std::vector<double>&
Tridiagonal::Solve()
{
    auto& buffer = Buffer();
    utils::SolveTridiagonalSystemBuffered(&lhs2Coeficients_.front(),
                                          &lhs1Coeficients_.front(),
                                          &lhs0Coeficients_.front(),
                                          &rightSideBuffer_.front(), numUnknowns_,
                                          &buffer.front());
    return rightSideBuffer_;
}

Tridiagonal::Tridiagonal(std::vector<double> lhs0Coeficients,
                         std::vector<double> lhs1Coeficients,
                         std::vector<double> lhs2Coeficients,
                         std::vector<std::vector<double>> rhsCoeficients,
                         const size_t numUnknowns,
                         const size_t problemSize)
        : luBuffer_(numUnknowns),
          rightSideBuffer_(numUnknowns),
          lhs0Coeficients_(std::move(lhs0Coeficients)),
          lhs1Coeficients_(std::move(lhs1Coeficients)),
          lhs2Coeficients_(std::move(lhs2Coeficients)),
          rhsCoeficients_(std::move(rhsCoeficients)),
          numUnknowns_(numUnknowns),
          problemSize_(problemSize)
{
    luBuffer_.assign(numUnknowns, 0);
    rightSideBuffer_.assign(numUnknowns, 0);
}


std::vector<double>& Tridiagonal::RightSideBuffer()
{
    return rightSideBuffer_;
}


const std::vector<double>& Tridiagonal::Lhs0Coeficients() const
{
    return lhs0Coeficients_;
}

const std::vector<double>& Tridiagonal::Lhs1Coeficients() const
{
    return lhs1Coeficients_;
}

const std::vector<double>& Tridiagonal::Lhs2Coeficients() const
{
    return lhs2Coeficients_;
}

const std::vector<std::vector<double>>& Tridiagonal::RhsCoeficients() const
{
    return rhsCoeficients_;
}

size_t Tridiagonal::NumUnknowns() const
{
    return numUnknowns_;
}

size_t Tridiagonal::ProblemSize() const
{
    return problemSize_;
}


Tridiagonal Tridiagonal::Factory::CreateEmptyTridiagonal()
{
    Tridiagonal tridiagonal(std::vector<double>(),
                            std::vector<double>(),
                            std::vector<double>(),
                            std::vector<std::vector<double>>(),
    0,
            0
    );
    return tridiagonal;
}

Tridiagonal
Tridiagonal::Factory::CreateFullTridiagonal(const KnotVector& knotVector,
                                            const size_t numKnots)
{
    std::vector<double> lhs0Coeficients;
    std::vector<double> lhs1Coeficients;
    std::vector<double> lhs2Coeficients;
    std::vector<double> rhs0Coeficients;
    std::vector<double> rhs1Coeficients;
    std::vector<double> rhs2Coeficients;

    const auto numUnknowns = numKnots - 2;
    lhs0Coeficients.reserve(numUnknowns);
    lhs1Coeficients.reserve(numUnknowns);
    lhs2Coeficients.reserve(numUnknowns);
    rhs0Coeficients.reserve(numUnknowns);
    rhs1Coeficients.reserve(numUnknowns);
    rhs2Coeficients.reserve(numUnknowns);

    for (int i = 1; i < numUnknowns + 1; ++i)
    {
        const auto h1 = knotVector[i + 1] - knotVector[i];
        const auto h0 = knotVector[i] - knotVector[i - 1];
        const auto h11 = h1 * h1;
        const auto h00 = h0 * h0;
        const auto h10 = h1 * h0;
        lhs2Coeficients.emplace_back(
                h1 * h00
        );

        lhs1Coeficients.emplace_back(
                2 * h10 * (h1 + h0)
        );

        lhs0Coeficients.emplace_back(
                h11 * h0
        );

        rhs2Coeficients.emplace_back(
                3 * h00
        );

        rhs1Coeficients.emplace_back(
                3 * (h11 - h00)
        );

        rhs0Coeficients.emplace_back(
                -3 * h11
        );
    }

    std::vector<std::vector<double>> rhsCoeficients =
            {
                    std::move(rhs0Coeficients),
                    std::move(rhs1Coeficients),
                    std::move(rhs2Coeficients)
            };
    Tridiagonal tridiagonal(
            std::move(lhs0Coeficients),
            std::move(lhs1Coeficients),
            std::move(lhs2Coeficients),
            std::move(rhsCoeficients),
            numUnknowns,
            knotVector.size()
    );
    return tridiagonal;
}

Tridiagonal
Tridiagonal::Factory::CreateReducedTridiagonal(const KnotVector& knotVector,
                                               const size_t numKnots)
{
    std::vector<double> lhs0Coeficients;
    std::vector<double> lhs1Coeficients;
    std::vector<double> lhs2Coeficients;
    std::vector<double> rhs0Coeficients;
    std::vector<double> rhs1Coeficients;
    std::vector<double> rhs2Coeficients;
    std::vector<double> rhs3Coeficients;
    std::vector<double> rhs4Coeficients;
    std::vector<double> restZ0Coeficients;
    std::vector<double> restZ1Coeficients;
    std::vector<double> restZ2Coeficients;
    std::vector<double> restD0Coeficients;
    std::vector<double> restD2Coeficients;

    const auto even = numKnots % 2 == 0;
    const auto numUnknowns = even ? numKnots / 2 - 2 : numKnots / 2 - 1;

    lhs0Coeficients.reserve(numUnknowns);
    lhs1Coeficients.reserve(numUnknowns);
    lhs2Coeficients.reserve(numUnknowns);
    rhs0Coeficients.reserve(numUnknowns);
    rhs1Coeficients.reserve(numUnknowns);
    rhs2Coeficients.reserve(numUnknowns);
    rhs3Coeficients.reserve(numUnknowns);
    rhs4Coeficients.reserve(numUnknowns);
    restZ0Coeficients.reserve(numUnknowns + 1);
    restZ1Coeficients.reserve(numUnknowns + 1);
    restZ1Coeficients.reserve(numUnknowns + 1);
    restD0Coeficients.reserve(numUnknowns + 1);
    restD2Coeficients.reserve(numUnknowns + 1);

    const auto until = even ? knotVector.size() - 2 : knotVector.size() - 1;
    auto i = 2;
    for (; i < until; i += 2)
    {
        const auto h3 = knotVector[i + 2] - knotVector[i + 1];
        const auto h2 = knotVector[i + 1] - knotVector[i];
        const auto h1 = knotVector[i] - knotVector[i - 1];
        const auto h0 = knotVector[i - 1] - knotVector[i - 2];
        const auto h00 = h0 * h0;
        const auto h11 = h1 * h1;
        const auto h22 = h2 * h2;
        const auto h20 = h0 * h2;
        const auto h21 = h1 * h2;
        const auto h31 = h1 * h3;
        const auto h100 = h1 * h00;
        const auto h110 = h11 * h0;
        const auto h322 = h3 * h22;
        const auto h3210 = h20 * h31;
        const auto h322110 = h3210 * h21;
        const auto h3Plush2 = h3 + h2;
        const auto h2Plush1 = h2 + h1;
        const auto h1Plush0 = h1 + h0;
        const auto rth1Plush0h10 = 1 / (2 * h1Plush0 * h1 * h0);

        lhs2Coeficients.emplace_back(
                h1Plush0 * h322110
        );

        lhs1Coeficients.emplace_back(
                h3210 * (h31 * h1Plush0 + h3Plush2 * (h20 - 4 * h1Plush0 * h2Plush1))

        );

        lhs0Coeficients.emplace_back(
                h3Plush2 * h322110
        );

        rhs4Coeficients.emplace_back(
                3 * h1Plush0 * h22 * h110
        );

        rhs3Coeficients.emplace_back(
                -3 * h1Plush0 * h3Plush2 * h3Plush2 * h110
        );

        rhs2Coeficients.emplace_back(
                3 * h3 * h0 * (h3Plush2 * h0 * h22 - h1Plush0 * (2 * h3Plush2 * (h2 - h1) * h2Plush1 + h3 * h11))
        );

        rhs1Coeficients.emplace_back(
                3 * h3Plush2 * h1Plush0 * h1Plush0 * h322
        );

        rhs0Coeficients.emplace_back(
                -3 * h3Plush2 * h322 * h11
        );

        restZ2Coeficients.emplace_back(
                3 * h00 * rth1Plush0h10
        );

        restZ1Coeficients.emplace_back(
                -3 * (h00 - h11) * rth1Plush0h10
        );

        restZ0Coeficients.emplace_back(
                -3 * h11 * rth1Plush0h10
        );

        restD2Coeficients.emplace_back(
                -h100 * rth1Plush0h10
        );

        restD0Coeficients.emplace_back(
                -h110 * rth1Plush0h10
        );
    }

    const auto h1 = knotVector[i] - knotVector[i - 1];
    const auto h0 = knotVector[i - 1] - knotVector[i - 2];
    const auto h00 = h0 * h0;
    const auto h11 = h1 * h1;
    const auto h100 = h1 * h00;
    const auto h110 = h11 * h0;
    const auto h1Plush0 = h1 + h0;
    const auto rth1Plush0h10 = 1 / (2 * h1Plush0 * h1 * h0);

    restZ2Coeficients.emplace_back(
            3 * h00 * rth1Plush0h10
    );

    restZ1Coeficients.emplace_back(
            -3 * (h00 - h11) * rth1Plush0h10
    );

    restZ0Coeficients.emplace_back(
            -3 * h11 * rth1Plush0h10
    );

    restD2Coeficients.emplace_back(
            -h100 * rth1Plush0h10
    );

    restD0Coeficients.emplace_back(
            -h110 * rth1Plush0h10
    );

    std::vector<std::vector<double>> rhsCoeficients =
            {
                    std::move(rhs0Coeficients),
                    std::move(rhs1Coeficients),
                    std::move(rhs2Coeficients),
                    std::move(rhs3Coeficients),
                    std::move(rhs4Coeficients),
                    std::move(restZ0Coeficients),
                    std::move(restZ1Coeficients),
                    std::move(restZ2Coeficients),
                    std::move(restD0Coeficients),
                    std::move(restD2Coeficients)
            };
    Tridiagonal tridiagonal(
            std::move(lhs0Coeficients),
            std::move(lhs1Coeficients),
            std::move(lhs2Coeficients),
            std::move(rhsCoeficients),
            numUnknowns,
            knotVector.size()
    );
    return tridiagonal;
}

Tridiagonal Tridiagonal::Factory::CreateFullWithDivisionsTridiagonal(const KnotVector& knotVector,
                                                                     const size_t numKnots)
{
    std::vector<double> lhs0Coeficients;
    std::vector<double> lhs1Coeficients;
    std::vector<double> lhs2Coeficients;

    const auto numUnknowns = numKnots - 2;
    lhs0Coeficients.reserve(numUnknowns);
    lhs1Coeficients.reserve(numUnknowns);
    lhs2Coeficients.reserve(numUnknowns);

    for (int i = 1; i < numUnknowns + 1; ++i)
    {
        const auto h1 = knotVector[i + 1] - knotVector[i];
        const auto h0 = knotVector[i] - knotVector[i - 1];

        lhs2Coeficients.emplace_back(
                h0
        );

        lhs1Coeficients.emplace_back(
                2 * (h1 + h0)
        );

        lhs0Coeficients.emplace_back(
                h1
        );
    }

    std::vector<std::vector<double>> rhsCoeficients;
    Tridiagonal tridiagonal(
            std::move(lhs0Coeficients),
            std::move(lhs1Coeficients),
            std::move(lhs2Coeficients),
            std::move(rhsCoeficients),
            numUnknowns,
            knotVector.size()
    );
    return tridiagonal;
}

Tridiagonal Tridiagonal::Factory::CreateReducedWithDivisionsTridiagonal(
        const KnotVector& knotVector, const size_t numKnots)
{
    std::vector<double> lhs0Coeficients;
    std::vector<double> lhs1Coeficients;
    std::vector<double> lhs2Coeficients;

    const auto even = numKnots % 2 == 0;
    const auto numUnknowns = even ? numKnots / 2 - 2 : numKnots / 2 - 1;

    lhs0Coeficients.reserve(numUnknowns);
    lhs1Coeficients.reserve(numUnknowns);
    lhs2Coeficients.reserve(numUnknowns);

    const auto until = even ? knotVector.size() - 2 : knotVector.size() - 1;
    auto i = 2;
    for (; i < until; i += 2)
    {
        const auto h3 = knotVector[i + 2] - knotVector[i + 1];
        const auto h2 = knotVector[i + 1] - knotVector[i];
        const auto h1 = knotVector[i] - knotVector[i - 1];
        const auto h0 = knotVector[i - 1] - knotVector[i - 2];
        const auto h20 = h0 * h2;
        const auto h21 = h1 * h2;
        const auto h31 = h1 * h3;
        const auto h3Plush2 = h3 + h2;
        const auto h2Plush1 = h2 + h1;
        const auto h1Plush0 = h1 + h0;

        lhs2Coeficients.emplace_back(
                h1Plush0
        );

        lhs1Coeficients.emplace_back(
                (h31*h1Plush0+h3Plush2*(h20-4*h1Plush0*h2Plush1))/h21
        );

        lhs0Coeficients.emplace_back(
                h3Plush2
        );

    }

    std::vector<std::vector<double>> rhsCoeficients;
    Tridiagonal tridiagonal(
            std::move(lhs0Coeficients),
            std::move(lhs1Coeficients),
            std::move(lhs2Coeficients),
            std::move(rhsCoeficients),
            numUnknowns,
            knotVector.size()
    );
    return tridiagonal;
}

double Tridiagonal::AccumulateExecutionTimes(Tridiagonals& tridiagonals)
{
    return 0;
}

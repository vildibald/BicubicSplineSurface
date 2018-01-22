#include "stdafx.h"
#include "utils.h"
#include <omp.h>
#include <vector>

namespace utils
{
    unsigned int logicalThreadCount = std::thread::hardware_concurrency();

    std::vector<double> SolveCsabaDeboorTridiagonalSystem(
            double main_diagonal_value, double right_side[], unsigned int
    num_equations, double last_main_diagonal_value)
    {
        std::vector<double> L(num_equations);
        std::vector<double> U(num_equations);
        std::vector<double> Y(num_equations);
        std::vector<double> D(num_equations);
        SolveCsabaDeboorTridiagonalSystemBuffered(main_diagonal_value,
                                                  right_side, num_equations, &L.front(),
                                                  &U.front(), &Y.front(), &D.front(),
                                                  last_main_diagonal_value);
        return D;
    }

    void SolveCsabaDeboorTridiagonalSystemBuffered(double main_diagonal_value,
                                                   double right_side[], unsigned int num_equations,
                                                   double L[],
                                                   double U[], double Y[], double D[],
                                                   double last_main_diagonal_value)
    {
        if (last_main_diagonal_value == DBL_MIN)
            last_main_diagonal_value = main_diagonal_value;
        U[0] = main_diagonal_value;
        for (size_t i = 1; i < num_equations - 1; i++)
        {
            L[i] = 1 / U[i - 1];
            U[i] = main_diagonal_value - L[i];
        }
        L[num_equations - 1] = 1 / U[num_equations - 1 - 1];
        U[num_equations - 1] = last_main_diagonal_value - L[num_equations - 1];
        Y[0] = right_side[0];
        for (int i = 1; i < num_equations; i++)
        {
            // Fw
            Y[i] = right_side[i] - L[i] * Y[i - 1];
        }
        D[num_equations - 1] = Y[num_equations - 1] / U[num_equations - 1];
        for (int i = num_equations - 2; 0 <= i; i--)
        {
            // Bw
            D[i] = (Y[i] - D[i + 1]) / U[i];
        }
    }

    void SolveTridiagonalSystemBuffered(double* lowerDiagonal, double*
    mainDiagonal, double* upperDiagonal,
                                        double* rightSide, const size_t numEquations,
                                        double* buffer)
    {
        std::copy(upperDiagonal, upperDiagonal + numEquations, buffer);
//		memcpy(buffer, upperDiagonal, numEquations);
        buffer[0] /= mainDiagonal[0];
        rightSide[0] /= mainDiagonal[0];
        for (size_t i = 1; i < numEquations; i++)
        {
            const auto m = 1 / (mainDiagonal[i] - lowerDiagonal[i] *
                                                  buffer[i - 1]);
            buffer[i] *= m;
            rightSide[i] = (rightSide[i] - lowerDiagonal[i] *
                                           rightSide[i - 1]) * m;
        }
        for (size_t i = numEquations - 1; i-- > 0;)
        {
            rightSide[i] -= buffer[i] * rightSide[i + 1];
        }
    }

    void SolveDeboorTridiagonalSystemBuffered(const double lowerDiagonalValue,
                                              const double mainDiagonalValue,
                                              const double upperDiagonalValue, double*
    rightSide, const size_t numEquations, double* buffer,
                                              double
                                              lastMainDiagonalValue)
    {
        if (lastMainDiagonalValue == DBL_MIN)
            lastMainDiagonalValue = mainDiagonalValue;
        auto m0 = 1 / mainDiagonalValue;
        buffer[0] = upperDiagonalValue * m0;
        rightSide[0] *= m0;
        const auto lastindex = numEquations - 1;
        for (size_t i = 1; i < lastindex; i++)
        {
            auto m = 1 / (mainDiagonalValue - lowerDiagonalValue *
                                              buffer[i - 1]);
            buffer[i] = upperDiagonalValue * m;
            rightSide[i] = (rightSide[i] - lowerDiagonalValue *
                                           rightSide[i - 1]) * m;
        }

        m0 = 1 / (lastMainDiagonalValue - lowerDiagonalValue *
                                          buffer[lastindex - 1]);
        buffer[lastindex] = upperDiagonalValue * m0;
        rightSide[lastindex] = (rightSide[lastindex] - lowerDiagonalValue *
                                                       rightSide[lastindex - 1]) * m0;

        for (size_t i = numEquations - 1; i-- > 0;)
        {
            rightSide[i] -= buffer[i] * rightSide[i + 1];
        }
    }

    void SolveDeboorTridiagonalSystemBuffered(double main_diagonal_value,
                                              double right_side[], size_t num_equations,
                                              double buffer[], double
                                              last_main_diagonal_value)
    {
        if (last_main_diagonal_value == DBL_MIN)
            last_main_diagonal_value = main_diagonal_value;
        auto m0 = 1 / main_diagonal_value;
        buffer[0] = m0;
        right_side[0] *= m0;
        auto lastindex = num_equations - 1;
        for (size_t i = 1; i < lastindex; i++)
        {
            auto m = 1 / (main_diagonal_value - buffer[i - 1]);
            buffer[i] = m;
            right_side[i] = (right_side[i] - right_side[i - 1]) * m;
        }

        m0 = 1 / (last_main_diagonal_value - buffer[lastindex - 1]);
        buffer[lastindex] = m0;
        right_side[lastindex] = (right_side[lastindex]
                                 - right_side[lastindex - 1]) * m0;

        for (size_t i = num_equations - 1; i-- > 0;)
        {
            right_side[i] -= buffer[i] * right_side[i + 1];
        }
    }

    void SolveDeboorTridiagonalSystem(double main_diagonal_value,
                                      double right_side[], size_t num_equations,
                                      double last_main_diagonal_value)
    {
        auto buffer = new double[num_equations];
        SolveDeboorTridiagonalSystemBuffered(main_diagonal_value, right_side,
                                             num_equations, buffer, last_main_diagonal_value);
        delete[] buffer;
    }
}

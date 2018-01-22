#pragma once

#include <omp.h>
#include <thread>
#include <vector>
#include <cfloat>

namespace utils
{
    extern unsigned int logicalThreadCount;

    template <typename Iterator, typename Function>
    void For(Iterator from, Iterator to, Iterator incrementBy, const bool inParallel,
             Function function)
    {
        if (inParallel)
        {
#pragma omp parallel for
            for (Iterator i = from; i < to; i += incrementBy)
            {
                function(i);
            }
        }
        else
        {
            // Loop for sequential computations.
            // '#pragma omp parallel for if(inParallel)' in case of inParallel==false will execute
            // such loop in one thread, but still with overhead of OpenMP thread creation.
            for (Iterator i = from; i < to; i += incrementBy)
            {
                function(i);
            }
        }
    }

    template <typename T>
    void DeleteJaggedArray(T** jaggedArray, const size_t rows, size_t columns)
    {
        for (size_t i = 0; i < rows; i++)
        {
            delete[] jaggedArray[i];
            jaggedArray[i] = nullptr;
        }
        delete[] jaggedArray;
        jaggedArray = nullptr;
    }

    template <typename T>
    T** CreateJaggedArray(const size_t rows, const size_t columns)
    {
        auto res = new T*[rows];
        for (size_t i = 0; i < columns; i++)
        {
            res[i] = new T[columns];
        }

        return res;
    }

    std::vector<double> SolveCsabaDeboorTridiagonalSystem(double b,
                                                          double rightSide[],
                                                          unsigned int numEquations,
                                                          double lastMainDiagonalValue = DBL_MIN);

    void SolveCsabaDeboorTridiagonalSystemBuffered(double mainDiagonalValue,
                                                   double rightSide[], unsigned int numEquations,
                                                   double l[],
                                                   double u[], double y[], double d[],
                                                   double lastMainDiagonalValue = DBL_MIN);

    void SolveTridiagonalSystem(double lower_diagonal[],
                                double main_diagonal[], double upper_diagonal[],
                                double right_side[],
                                size_t num_equations);

    void SolveTridiagonalSystemBuffered(double lower_diagonal[],
                                        double main_diagonal[],
                                        double upper_diagonal[], double right_side[],
                                        size_t num_equations,
                                        double buffer[]);

    void SolveDeboorTridiagonalSystem(double lower_diagonal_value,
                                      double main_diagonal_value, double upper_diagonal_value,
                                      double right_side[], size_t num_equations,
                                      double last_main_diagonal_value = DBL_MIN);

    void SolveDeboorTridiagonalSystemBuffered(double lower_diagonal_value,
                                              double main_diagonal_value,
                                              double upper_diagonal_value,
                                              double right_side[], size_t num_equations,
                                              double buffer[],
                                              double last_main_diagonal_value = DBL_MIN);

    void SolveDeboorTridiagonalSystem(double main_diagonal_value,
                                      double right_side[], size_t num_equations,
                                      double last_main_diagonal_value = DBL_MIN);

    void SolveDeboorTridiagonalSystemBuffered(double main_diagonal_value,
                                              double right_side[], size_t num_equations,
                                              double buffer[], double
                                              last_main_diagonal_value = DBL_MIN);

    template <typename T>
    double Average(T arr[], size_t arr_size);
};

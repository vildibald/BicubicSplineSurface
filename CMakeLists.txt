cmake_minimum_required(VERSION 3.9)
project(BicubicSplineSurface)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -fopenmp")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=core2")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")

add_executable(BicubicSplineSurface main.cpp

        stdafx.cpp stdafx.h SplineKnots.h Spline.cpp Spline.h KnotVector.h KnotMatrix.h utils.cpp utils.h ComparisonBenchmarkResult.cpp ComparisonBenchmarkResult.h StopWatch.h MathFunction.cpp MathFunction.h FullAlgorithm.cpp FullAlgorithm.h Tridiagonal.cpp Tridiagonal.h MultithreadPreparator.h Timer.h CalculationMode.h ReducedAlgorithm.cpp ReducedAlgorithm.h)
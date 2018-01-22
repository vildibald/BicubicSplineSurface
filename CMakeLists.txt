cmake_minimum_required(VERSION 3.9)
project(BicubicSplineSurface)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")

#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=core2")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=nehalem")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=sandybridge")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=haswell")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=skylake")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=barcelona")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=bdver1")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=bdver4")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=znver1")

add_executable(BicubicSplineSurface main.cpp

        stdafx.cpp stdafx.h SplineKnots.h Spline.cpp Spline.h KnotVector.h KnotMatrix.h utils.cpp utils.h ComparisonBenchmarkResult.cpp ComparisonBenchmarkResult.h StopWatch.h MathFunction.cpp MathFunction.h FullAlgorithm.cpp FullAlgorithm.h Tridiagonal.cpp Tridiagonal.h MultithreadPreparator.h Timer.h CalculationMode.h ReducedAlgorithm.cpp ReducedAlgorithm.h)
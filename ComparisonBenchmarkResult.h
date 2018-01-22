#pragma once

#include <vector>

class ComparisonBenchmarkResult {
    std::vector<unsigned long long> algorithmTimes;

public:
    ComparisonBenchmarkResult();

    ComparisonBenchmarkResult& Add(unsigned long long time);

    unsigned long long operator[](size_t i) const;

    double Ratio(size_t i, size_t j) const;
};

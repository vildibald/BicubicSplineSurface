#include "stdafx.h"
#include "ComparisonBenchmarkResult.h"


ComparisonBenchmarkResult::ComparisonBenchmarkResult()
        : algorithmTimes() {
}

ComparisonBenchmarkResult& ComparisonBenchmarkResult::Add(unsigned long long time) {
    algorithmTimes.emplace_back(time);
    return *this;
}

unsigned long long ComparisonBenchmarkResult::operator[](const size_t i) const {
    return algorithmTimes[i];
}

double ComparisonBenchmarkResult::Ratio(const size_t i, const size_t j) const {
    return static_cast<double>(algorithmTimes[i]) / static_cast<double>(algorithmTimes[j]);
}

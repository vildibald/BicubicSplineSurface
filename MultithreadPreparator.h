#pragma once
#include <vector>
#include <thread>

class MultithreadPreparator {

public:

    template<typename T>
    void PrepareVector(std::vector<T>& vector) {
        auto threads = std::thread::hardware_concurrency();
        vector.reserve(threads);
        for (auto i = vector.size(); i < threads; i++) {
            // create copy of tridiagonal solver
            vector.emplace_back(vector[0]);
        }
    }

    template<typename T>
    void UnprepareVector(std::vector<T>& vector) {
        for (int i = 0; i < vector.size() - 1; ++i) {
            vector.pop_back();
        }
    }

    template<typename T>
    void PrepareVector(const bool inParallel, std::vector<T>& vector) {
        if (inParallel) {
            PrepareVector(vector);
        }
        else {
            UnprepareVector(vector);
        }
    }
};

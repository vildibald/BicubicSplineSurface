#pragma once
enum CalculationMode
{
    OPTIMIZED_DIVISIONS_BUFFERED,
    OPTIMIZED_DIVISIONS,
    NON_OPTIMIZED
};

inline CalculationMode ToCalculationMode(int value)
{
    if(value < OPTIMIZED_DIVISIONS_BUFFERED || value > NON_OPTIMIZED)
    {
        return OPTIMIZED_DIVISIONS_BUFFERED;
    }
    return static_cast<CalculationMode>(value);
}
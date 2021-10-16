#pragma once

#include "resistor.h"
#include <array>
#include <fmt/format.h>
#include <fmt/compile.h>
#include <iostream>

#define FORMAT(fmtStr, ...) fmt::format(FMT_COMPILE(fmtStr), __VA_ARGS__)

struct identity
{
    template<typename T>
    T&& operator() (T&& val) const
    {
        return std::forward<T>(val);
    }
};

template<typename T, std::size_t N, typename Projection>
void formatArray(std::ostream& os, const std::array<T, N>& arr, Projection proj)
{
    os << '[';

    if (arr.size() > 0)
    {
        for (std::size_t i = 0; i < arr.size() - 1; ++i)
            os << FORMAT("{: 8.4f}, ", proj(arr[i]));
    }

    os << FORMAT("{: 8.4f}]", proj(arr.back()));
}

template<std::size_t N>
std::ostream& operator<< (std::ostream& os, const std::array<double, N>& arr)
{
    formatArray(os, arr, identity());
    return os;
}

template<std::size_t N>
std::ostream& operator<< (std::ostream& os, const std::array<Resistor, N>& arr)
{
    formatArray(os, arr, [](const Resistor& r) { return r.getValue(); });
    return os;
}


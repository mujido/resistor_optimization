#pragma once

#include <array>
#include <cmath>

unsigned log2(unsigned value);

template<typename Container, typename T, std::size_t N, typename Projection>
std::size_t closest(Container&& cont, const T& val, Projection&& proj)
{
    static_assert(N > 0);

    T minDiff = std::abs(proj(cont[0]) - val);
    std::size_t minIndex = 0;

    for (std::size_t i = 1; i < cont.size(); ++i)
    {
        double diff = std::abs(proj(cont[i]) - val);
        if (diff < minDiff)
        {
            minDiff = diff;
            minIndex = i;
        }
    }

    return minIndex;
}
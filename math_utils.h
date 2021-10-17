#pragma once

#include <array>
#include <cmath>

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

template<typename T, std::size_t N>
std::array<T, N> operator- (const std::array<T, N>& lhs, const std::array<T, N>& rhs)
{
    auto result = lhs;
    for (unsigned i = 0; i < lhs.size(); ++i)
        result[i] -= rhs[i];

    return result;
}


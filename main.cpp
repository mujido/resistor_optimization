#include <boost/format.hpp>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <random>
#include <valarray>

#include "math_utils.h"
#include "model.h"
#include "resistor_series.h"
#include "resistor.h"
#include "signals.h"

using boost::format;

constexpr std::array<double, 13> targetPotentials{
    2.0722,
    2.2126,
    2.2848,
    2.7874,
    2.9278,
    3,
    3.0722,
    3.2126,
    3.2848,
    3.7152,
    3.7874,
    3.9278,
    4,
};

template<typename T, std::size_t N>
std::ostream& operator<< (std::ostream& os, const std::array<T, N>& arr)
{
    auto fmt = boost::format("% 8.4f");

    os << '[';

    if (arr.size() > 0)
    {
        for (std::size_t i = 0; i < arr.size() - 1; ++i)
            os << fmt % arr[i] << ", ";
    }

    return os << fmt % arr.back() << ']';
}

template<unsigned thisCount>
void runArbitraryAt(unsigned rCount)
{
    if (thisCount <= 0)
        throw std::runtime_error("reached 0");

    if (thisCount == rCount)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r()};
        std::mt19937 rnd(seed);

        const auto& resistorSeries = e12Series;
        auto model = makeModel<thisCount, resistorSeries.decadeSize_>(targetPotentials, resistorSeries, rnd);
        // auto startingResistors = model.valuesToResistors({ 2700, 1200, 3900, 680, 68, 220, 270 });
        auto startingResistors = model.getRandomResistors();

        model.run(startingResistors);
    }
    else if constexpr (thisCount > 1)
        runArbitraryAt<thisCount - 1>(rCount);
}

constexpr unsigned MAX_RESISTOR_COUNT = 20;

void runArbitrary(unsigned rCount)
{
    runArbitraryAt<MAX_RESISTOR_COUNT>(rCount);
}

int main(int argc, const char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: resistors <resistor-count>" << std::endl;
        return 1;
    }

    int rCount = std::atoi(argv[1]);
    if (rCount <= 0 || rCount > MAX_RESISTOR_COUNT)
    {
        std::cerr << "resistor count out of range" << std::endl;
        return 1;
    }

    setupSignalHandler();

    runArbitrary(static_cast<unsigned>(rCount));
    return 0;
}
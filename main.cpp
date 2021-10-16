#include <fmt/format.h>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <random>
#include <valarray>

#include "math_utils.h"
#include "model.h"
#include "resistor_series.h"
#include "resistor.h"
#include "signals.h"

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

constexpr unsigned MAX_RESISTOR_COUNT = 20;

template<typename Series>
std::vector<Resistor> getRandomResistors(Series& series, unsigned count, RandomGen& rnd)
{
    std::vector<Resistor> resistors;
    resistors.reserve(count);

    std::uniform_int_distribution<std::size_t> resistorDomain(0, series.fullRange_.size() - 1);

    for (unsigned i = 0; i < count; ++i)
        resistors.push_back(series.fullRange_[resistorDomain(rnd)]);

	return resistors;
}

template<typename RSeries>
using ModelCreatorPtr = std::unique_ptr<IModel>(*)(const RSeries&, RandomGen&);

template<typename RSeries>
using ModelCreators = std::vector<ModelCreatorPtr<RSeries>>;

template<typename Series, unsigned... resistorCounts>
ModelCreators<Series> makeModelCreators(std::integer_sequence<unsigned, 0, resistorCounts...>)
{
    return ModelCreators<Series>{
        [](const Series& s, RandomGen& r) { return makeModel<resistorCounts, Series::decadeSize_>(targetPotentials, s, r); }...
    };
}

template<typename Series>
std::unique_ptr<IModel> getModelFor(unsigned resistorCount, const Series& resistorSeries, RandomGen& rnd)
{
    using CountSequence = std::make_integer_sequence<unsigned, MAX_RESISTOR_COUNT>;
    static const auto modelCreators = makeModelCreators<Series>(CountSequence());

    return modelCreators.at(resistorCount - 1)(resistorSeries, rnd);
}

void runArbitrary(unsigned rCount)
{
	std::random_device r;
	std::seed_seq seed{r(), r(), r(), r(), r(), r()};
	std::mt19937 rnd(seed);

	const auto& resistorSeries = e12Series;
	//auto model = makeModel<thisCount, resistorSeries.decadeSize_>(targetPotentials, resistorSeries, rnd);
	auto model = getModelFor(rCount, resistorSeries, rnd);
	// auto startingResistors = model.valuesToResistors({ 2700, 1200, 3900, 680, 68, 220, 270 });
	auto startingResistors = getRandomResistors(resistorSeries, rCount, rnd);

	model->run(startingResistors);
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
#pragma once

#include "math_utils.h"
#include "resistor.h"
#include <array>

template<std::size_t decadeSize, std::size_t decadeCount>
class ResistorSeries
{
public:
    using Base = std::array<double, decadeSize>;
    using FullRange = std::array<Resistor, decadeSize * decadeCount>;

    static constexpr std::size_t decadeSize_ = decadeSize;

    Base base_;
    FullRange fullRange_;
    FullRange reciprocals_;

    explicit ResistorSeries(const Base& base) :
        base_(base),
        fullRange_(createFullRange(base))
    {

    }

    const Resistor& getClosest(double value)
    {
        return closest(fullRange_, value, [](const Resistor& r) { return r.getValue(); });
    }

private:
    static FullRange createFullRange(const Base& bases)
    {
        double mult = 1.0;

        FullRange result{};

        for (std::size_t i = 0; i < decadeCount; ++i)
        {
            for (std::size_t j = 0; j < decadeSize; ++j)
            {
                auto index = decadeSize * i + j;
                result[index] = Resistor(index, bases[j] * mult);
            }

            mult *= 10.0;
        }

        return result;
    }

	static FullRange reciprocate(const FullRange& src)
	{
		FullRange result{};

		for (std::size_t i = 0; i < N; ++i)
			result[i] = 1.0 / src[i];

		return result;
	}
};

using E12SeriesType = ResistorSeries<12, 6>;
using E48SeriesType = ResistorSeries<48, 6>;

extern E12SeriesType e12Series;
extern E48SeriesType e48Series;

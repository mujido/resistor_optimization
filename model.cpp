#include "model.h"

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

auto model = makeModel<5, 12>(targetPotentials, e12Series, RandomGen());

#include "math_utils.h"
#include <array>

constexpr std::array<unsigned, 256> createLog2LookupTable()
{
    std::array<unsigned, 256> table{};

    for (unsigned i = 0; i < 8; ++i)
        for (unsigned j = 1U << i; j < (1U << (i + 1)); ++j)
            table[j] = i;

    return table;
}

unsigned log2(unsigned value)
{
    static constexpr auto lookupTable = createLog2LookupTable();

    unsigned tt;
    if ((tt = value >> 24) != 0)
        return 24 + lookupTable[tt];
    else if ((tt = value >> 16) != 0)
        return 16 + lookupTable[tt];
    else if ((tt = value >> 8) != 0)
        return 8 + lookupTable[tt];
    else
        return lookupTable[value];
}


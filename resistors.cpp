#include <boost/format.hpp>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <random>
#include <valarray>

#include "signals.h"

using boost::format;

static constexpr std::array<double, 13> targetPotentials{
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

static constexpr std::array<double, 12> e12Bases = { 1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2 };
static constexpr std::array<double, 48> e48Bases = {
    1.00, 1.05, 1.10, 1.15, 1.21, 1.27, 1.33, 1.40,
    1.47, 1.54, 1.62, 1.69, 1.78, 1.87, 1.96, 2.05, 
    2.15, 2.26, 2.37, 2.49, 2.61, 2.74, 2.87, 3.01, 
    3.16, 3.32, 3.48, 3.65, 3.83, 4.02, 4.22, 4.42, 
    4.64, 4.87, 5.11, 5.36, 5.62, 5.90, 6.19, 6.49,
    6.81, 7.15, 7.50, 7.87, 8.25, 8.66, 9.09, 9.53
};

using Ratios = std::array<double, targetPotentials.size()>;

template<std::size_t N>
constexpr std::array<double, N*6> createFullRange(const std::array<double, N>& bases)
{
    double mult = 1.0;

    std::array<double, N*6> result{};

    for (std::size_t i = 0; i < 6; ++i)
    {
        for (std::size_t j = 0; j < N; ++j)
            result[N*i + j] = bases[j] * mult;

        mult *= 10.0;
    }

    return result;
}

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

template<std::size_t N>
constexpr std::array<double, N> reciprocate(const std::array<double, N>& src)
{
    std::array<double, N> result{};

    for (std::size_t i = 0; i < N; ++i)
        result[i] = 1.0 / src[i];

    return result;
}

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

// const auto resistorDomain = createFullRange(e12Bases);
constexpr auto& resistorDomainBases = e12Bases;
constexpr auto resistorDomain = createFullRange(resistorDomainBases);
constexpr auto resistorReciprocals = reciprocate(resistorDomain);

class Resistor
{
public:
    constexpr Resistor() : index_(std::numeric_limits<std::size_t>::max()), value_(std::numeric_limits<double>::quiet_NaN())
    {}

    Resistor(std::size_t index) : index_(index), value_(resistorDomain[index_])
    {}

    std::size_t getIndex() const { return index_; }

    double getValue() const { return value_; }

    double getReciprocal() const
    {
        return resistorReciprocals[getIndex()];
    }

private:
    std::size_t index_;
    double value_;
};

std::ostream& operator<< (std::ostream& os, const Resistor& r)
{
    return os << r.getValue();
}

template<typename T, std::size_t N>
std::size_t closest(const std::array<T, N>& arr, const T& val)
{
    static_assert(N > 0);

    T minDiff = std::abs(arr[0] - val);
    std::size_t minIndex = 0;

    for (std::size_t i = 1; i < arr.size(); ++i)
    {
        double diff = std::abs(arr[i] - val);
        if (diff < minDiff)
        {
            minDiff = diff;
            minIndex = i;
        }
    }

    return minIndex;
}

template<unsigned resistorCount, unsigned bestMatchSearchDistance>
struct Model 
{
    using RandomGen = std::mt19937;

    RandomGen& rnd_;
    std::uniform_int_distribution<std::size_t> resistorCountDist_;
    std::uniform_int_distribution<std::size_t> resistorIndexDist_;
    std::uniform_int_distribution<std::size_t> resistorDomainIndexDist_;

//     constexpr auto resistorDomain = [](const auto& bases)
// {
//     double mult = 1.0;

//     std::array<double, 12*6> result;
//     auto resultIter = result.begin();

//     for (int i = 0; i < 6; ++i)
//     {
//         for (auto val : bases)
//             *resultIter++ = val * mult;

//         mult *= 10.0;
//     }

//     return result;
// }(e12Bases);

    struct AssignedRatio
    {
        double calculatedRatio_ = 0;
        std::bitset<resistorCount> resistorStates_{};
    };

    using ResistorArray = std::array<Resistor, resistorCount>;
    using AssignedRatios = std::array<AssignedRatio, std::tuple_size<Ratios>::value>;

    Model(RandomGen& rnd) :
        rnd_(rnd),
        resistorCountDist_(resistorCount - 2, resistorCount),
        resistorIndexDist_(0, resistorCount - 1),
        resistorDomainIndexDist_(0, resistorDomain.size() - 1)
    {}

    static AssignedRatios getBestRatios(const Ratios& desiredRatios, const ResistorArray& resistors)
    {
        const unsigned allOn = (1 << resistors.size()) - 1;

        AssignedRatios assigned{};

        for (unsigned i = 0; i <= allOn; ++i)
        {
            double ratio = 0.0;
            if (i == 0)
                ratio = 0.0;
            else if (i == allOn)
                ratio = 1.0;
            else
            {
                double top = 0.0;
                double bottom = 0.0;

                for (unsigned j = 0; j < resistors.size(); ++j)
                {
                    if ((i & (1 << j)) != 0)
                        top += resistors[j].getReciprocal();
                    else
                        bottom += resistors[j].getReciprocal();
                }

                ratio = top / bottom;
            }

            for (unsigned j = 0; j < desiredRatios.size(); ++j)
                if (std::abs(desiredRatios[j] - ratio) < std::abs(desiredRatios[j] - assigned[j].calculatedRatio_))
                {
                    assigned[j].calculatedRatio_ = ratio;
                    assigned[j].resistorStates_ = std::bitset<resistorCount>(i);
                }
        }

        return assigned;
    }

    static AssignedRatios getBestRatiosGrey(const Ratios& desiredRatios, const ResistorArray& resistors)
    {
        static_assert(resistorCount < std::numeric_limits<unsigned>::digits);

        const unsigned comboCount = 1 << resistors.size();

        AssignedRatios assigned{};
        unsigned prevBits = 0;
        double top = 0.0;
        double bottom = 0.0;

        // Starting with all resistors on bottom
        for (auto& r : resistors)
            bottom += r.getReciprocal();

        for (unsigned comboIndex = 0; comboIndex < comboCount; ++comboIndex)
        {
            // Use greycode to reduce number of calculations per iteration
            unsigned bits = comboIndex ^ (comboIndex >> 1);

            double ratio = 0.0;
            if (bits == 0)
                ratio = 0.0;
            else if (~bits == 0)
                ratio = 1.0;
            else
            {
                unsigned changedBit = prevBits ^ bits;
                unsigned resistorPos = log2(changedBit);
                double v = resistors[resistorPos].getReciprocal();

                // std::cout << boost::format("%8x %8x %d") % prevBits % bits % resistorPos << std::endl;

                if ((bits & (1 << resistorPos)) != 0)
                {
                    top += v;
                    bottom -= v;
                }
                else
                {
                    top -= v;
                    bottom += v;
                }

                ratio = top / bottom;
            }

            for (unsigned j = 0; j < desiredRatios.size(); ++j)
            {
                if (std::abs(desiredRatios[j] - ratio) < std::abs(desiredRatios[j] - assigned[j].calculatedRatio_))
                {
                    assigned[j].calculatedRatio_ = ratio;
                    assigned[j].resistorStates_ = std::bitset<resistorCount>(bits);
                }
            }

            prevBits = bits;
        }

        return assigned;
    }

    static double checkAccuracy(const Ratios& desiredRatios, const ResistorArray& resistors)
    {
        auto ratios = getBestRatiosGrey(desiredRatios, resistors);

        double error = 0.0;
        for (unsigned i = 0; i < ratios.size(); ++i)
            error += std::abs((desiredRatios[i] - ratios[i].calculatedRatio_) / desiredRatios[i]);

        return 100.0 * error;
    }

    static std::tuple<ResistorArray, double> findBiggestChange(const Ratios& desiredRatios, const ResistorArray& resistors, double currentScore)
    {
        int bestIndex = -1;
        Resistor bestResistor;
        double bestScore = currentScore;

        for (unsigned i = 0; i < resistors.size(); ++i)
        {
            const auto& r = resistors[i];
            ResistorArray rcopy = resistors;

            const auto evalChange = [&, i](int offsetLow, int offsetHigh)
            {
                for (int idxOffset = offsetLow; idxOffset <= offsetHigh; ++idxOffset)
                {
                    int newIdx = static_cast<int>(r.getIndex()) + idxOffset;
                    if (newIdx < 0 || static_cast<unsigned>(newIdx) >= resistorDomain.size())
                        continue;

                    Resistor newR(newIdx);
                    rcopy[i] = newR;

                    double newScore = checkAccuracy(desiredRatios, rcopy);
                    if (newScore < bestScore)
                    {
                        bestIndex = i;
                        bestResistor = newR;
                        bestScore = newScore;
                    }
                }
            };

            evalChange(-static_cast<int>(bestMatchSearchDistance), -1);
            evalChange(1, bestMatchSearchDistance);
        }

        if (bestIndex != -1)
        {
            auto rcopy = resistors;
            rcopy[bestIndex] = bestResistor;
            return {std::move(rcopy), bestScore};
        }
        else
            return {resistors, bestScore};
    }

    ResistorArray getRandomResistors()
    {
        ResistorArray resistors{};

        for (auto& r : resistors)
            r = Resistor(resistorDomainIndexDist_(rnd_));

        return resistors;
    }

    static ResistorArray valuesToResistors(const std::array<double, resistorCount>& resistorValues)
    {
        ResistorArray resistors{};

        for (std::size_t i = 0; i < resistorCount; ++i)
            resistors[i] = Resistor(closest(resistorDomain, resistorValues[i]));

        return resistors;
    }

    void printTable(const Ratios& targetRatios, AssignedRatios& bestRatios, const ResistorArray& bestResistors)
    {
        std::sort(bestRatios.begin(), bestRatios.end(), [](const AssignedRatio& lhs, const AssignedRatio& rhs)
        {
            return lhs.calculatedRatio_ > rhs.calculatedRatio_;
        });

        for (unsigned i = 0; i < targetRatios.size(); ++i)
        {
            std::cout << boost::format("% 12.7f  % 12.7f  % 12.7f  ") 
                % bestRatios[i].calculatedRatio_ 
                % (5.0 / (1.0 + bestRatios[i].calculatedRatio_))
                % (targetRatios[i] - bestRatios[i].calculatedRatio_);
            for (int r = 0; r < resistorCount; ++r)
                std::cout << ' ' << (bestRatios[i].resistorStates_.test(r) ? 1 : 0);

            std::cout << std::endl;
        }
    }

    void run(const Ratios& targetRatios, const ResistorArray& startingResistors)
    {
        double startingScore = checkAccuracy(targetRatios, startingResistors); 

        auto bestResistors = startingResistors;
        auto bestScore = startingScore;

        auto currentResistors = bestResistors;
        double currentScore = bestScore;
        unsigned randomizations = 0;

        auto lastPrint = std::chrono::steady_clock::now();

        auto printProgress = [&]()
        {
            std::cout << format("[%s] Score: %8.4f (%+8.4f) best: %8.4f") 
                % randomizations
                % currentScore
                % (currentScore - bestScore)
                % bestScore << "\n";
        };

        while (!exitedRequested())
        {
            auto [newResistors, newScore] = findBiggestChange(targetRatios, currentResistors, currentScore);

            if (newScore < currentScore)
            {
                currentResistors = newResistors;
                currentScore = newScore;
            }
            else
            {
                if (currentScore < bestScore)
                {
                    bestResistors = currentResistors;
                    bestScore = currentScore;
                }

                // printProgress();
                currentResistors = bestResistors;

                auto randCount = resistorIndexDist_(rnd_);
                std::bitset<resistorCount> switched{};
                while (randCount > 0)
                {
                    auto idx = resistorIndexDist_(rnd_);
                    if (!switched.test(idx))
                    {
                        currentResistors[idx] = Resistor(resistorDomainIndexDist_(rnd_));
                        switched.set(idx);
                        --randCount;
                    }
                }

                randomizations += 1;
                currentScore = 1000000;
            }

            auto now = std::chrono::steady_clock::now();
            if (now - lastPrint > std::chrono::seconds(2))
            {
                lastPrint = now;
                printProgress();
            }
        }

        std::sort(bestResistors.begin(), bestResistors.end(), [](const Resistor& lhs, const Resistor& rhs)
        {
            return lhs.getValue() > rhs.getValue();
        });

        auto bestRatios = getBestRatios(targetRatios, bestResistors);
        bestScore = checkAccuracy(targetRatios, bestResistors);

        auto percentDiffs = targetRatios;
        auto pdMean = 0.0;
        for (unsigned i = 0; i < percentDiffs.size(); ++i)
        {
            percentDiffs[i] -= bestRatios[i].calculatedRatio_;
            percentDiffs[i] = 100.0 * std::abs(percentDiffs[i]) / targetRatios[i];
            pdMean += percentDiffs[i];
        }

        pdMean /= percentDiffs.size();

        double pdStdDev = 0.0;
        for (auto pd : percentDiffs)
        {
            double x = pd - pdMean;
            pdStdDev += x * x;
        }

        pdStdDev = std::sqrt(pdStdDev / percentDiffs.size());

        Ratios bestRatioValues;
        for (unsigned i = 0; i < bestRatioValues.size(); ++i)
            bestRatioValues[i] = bestRatios[i].calculatedRatio_;

        std::cout << "Targets: " << targetRatios << std::endl;
        std::cout << "Best:    " << bestRatioValues << std::endl;
        std::cout << "% diff:  " << percentDiffs << std::endl;
        std::cout << (format("Score: %10.4f (%10.4f)  stddev: %.3g") % bestScore % (bestScore - startingScore) % pdStdDev) << std::endl;
        std::cout << bestResistors << std::endl;

        printTable(targetRatios, bestRatios, bestResistors);
    }

};

template<typename T, std::size_t N>
std::array<T, N> operator- (const std::array<T, N>& lhs, const std::array<T, N>& rhs)
{
    auto result = lhs;
    for (unsigned i = 0; i < lhs.size(); ++i)
        result[i] -= rhs[i];

    return result;
}

// template<typename T, std::size_t N>
// std::size_t argmin(const std::array<T, N>& arr)
// {
//     static_assert(N > 0);

//     T minValue = arr[0];
//     std::size_t minIndex = 0;

//     for (std::size_t i = 1; i < arr.size(); ++i)
//     {
//         if (arr[i] < minValue)
//         {
//             minValue = arr[i];
//             minIndex = i;
//         }
//     }

//     return minIndex;
// }

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

        Model<thisCount, resistorDomainBases.size()> model(rnd);
        // auto startingResistors = model.valuesToResistors({ 2700, 1200, 3900, 680, 68, 220, 270 });
        auto startingResistors = model.getRandomResistors();

        Ratios targetRatios{};
        for (unsigned i = 0; i < targetPotentials.size(); ++i)
            targetRatios[i] = 5.0 / targetPotentials[i] - 1.0;

        model.run(targetRatios, startingResistors);
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
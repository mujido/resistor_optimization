#pragma once

#include "formatters.h"
#include "resistor.h"
#include "resistor_series.h"
#include "signals.h"
#include <bitset>
#include <functional>
#include <random>

template<std::size_t N>
using Potentials = std::array<double, N>;

template<std::size_t N>
using Ratios = std::array<double, N>;

using RandomGen = std::mt19937;

struct IModel
{
    virtual ~IModel() = default;
    virtual void run(const std::vector<Resistor>& startingResistors) = 0;
};

template<
    unsigned networkResistorCount,
    unsigned bestMatchSearchDistance,
    std::size_t targetCount,
    typename ResistorSeriesType>
struct Model : IModel
{
    struct AssignedRatio
    {
        double calculatedRatio_ = 0;
        std::bitset<networkResistorCount> resistorStates_{};
    };

    using ResistorArray = std::array<Resistor, networkResistorCount>;
    using AssignedRatios = std::array<AssignedRatio, targetCount>;
    using TargetPotentials = Potentials<targetCount>;
    using TargetRatios = Ratios<targetCount>;


    TargetRatios targetRatios_;
    const ResistorSeriesType& resistorSeries_;
    RandomGen& rnd_;
    std::uniform_int_distribution<std::size_t> resistorCountDist_;
    std::uniform_int_distribution<std::size_t> resistorIndexDist_;
    std::uniform_int_distribution<std::size_t> resistorDomainIndexDist_;

    Model(const TargetPotentials& targetPotentials, const ResistorSeriesType& resistorSeries, RandomGen& rnd) : 
        targetRatios_(calcTargetRatios(targetPotentials)),
        resistorSeries_(resistorSeries),
        rnd_(rnd),
        resistorCountDist_(networkResistorCount - 2, networkResistorCount),
        resistorIndexDist_(0, networkResistorCount - 1),
        resistorDomainIndexDist_(0, resistorSeries.fullRange_.size() - 1)
    {
    }

    static TargetRatios calcTargetRatios(const TargetPotentials& targetPotentials)
    {
        TargetRatios targetRatios{};

        for (unsigned i = 0; i < targetPotentials.size(); ++i)
            targetRatios[i] = 5.0 / targetPotentials[i] - 1.0;

        return targetRatios;
    }

    AssignedRatios getBestRatios(const ResistorArray& resistors)
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

            for (unsigned j = 0; j < targetRatios_.size(); ++j)
                if (std::abs(targetRatios_[j] - ratio) < std::abs(targetRatios_[j] - assigned[j].calculatedRatio_))
                {
                    assigned[j].calculatedRatio_ = ratio;
                    assigned[j].resistorStates_ = std::bitset<networkResistorCount>(i);
                }
        }

        return assigned;
    }

    static AssignedRatios getBestRatiosGrey(const TargetRatios& desiredRatios, const ResistorArray& resistors)
    {
        static_assert(networkResistorCount < std::numeric_limits<unsigned>::digits);

        constexpr unsigned comboCount = 1 << networkResistorCount;
        constexpr unsigned resistorBitMask = (1 << networkResistorCount) - 1;

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
            else if ((bits & resistorBitMask) == resistorBitMask)
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

            for (unsigned j = 0; j < targetCount; ++j)
            {
                if (std::abs(desiredRatios[j] - ratio) < std::abs(desiredRatios[j] - assigned[j].calculatedRatio_))
                {
                    assigned[j].calculatedRatio_ = ratio;
                    assigned[j].resistorStates_ = std::bitset<networkResistorCount>(bits);
                }
            }

            prevBits = bits;
        }

        return assigned;
    }

    double checkAccuracy(const ResistorArray& resistors)
    {
        auto ratios = getBestRatiosGrey(targetRatios_, resistors);

        double error = 0.0;
        for (unsigned i = 0; i < ratios.size(); ++i)
            error += std::abs((targetRatios_[i] - ratios[i].calculatedRatio_) / targetRatios_[i]);

        return 100.0 * error;
    }

    std::tuple<ResistorArray, double> findBiggestChange(const ResistorArray& resistors, double currentScore)
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
                    if (newIdx < 0 || static_cast<unsigned>(newIdx) >= resistorSeries_.fullRange_.size())
                        continue;

                    auto newR = resistor(newIdx);
                    rcopy[i] = newR;

                    double newScore = checkAccuracy(rcopy);
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

    Resistor resistor(std::size_t index)
    {
        return resistorSeries_.fullRange_[index];
    }

    ResistorArray valuesToResistors(const std::array<double, networkResistorCount>& resistorValues)
    {
        ResistorArray resistors{};

        for (std::size_t i = 0; i < networkResistorCount; ++i)
            resistors[i] = resistorSeries_.getClosest(resistorValues[i]);

        return resistors;
    }

    void printTable(AssignedRatios& bestRatios, const ResistorArray& bestResistors)
    {
        std::sort(bestRatios.begin(), bestRatios.end(), [](const AssignedRatio& lhs, const AssignedRatio& rhs)
        {
            return lhs.calculatedRatio_ > rhs.calculatedRatio_;
        });

        for (unsigned i = 0; i < targetRatios_.size(); ++i)
        {
            std::cout << FORMAT("{: 12.7f}  {: 12.7f}  {: 12.7f}  ",
                bestRatios[i].calculatedRatio_,
                5.0 / (1.0 + bestRatios[i].calculatedRatio_),
                targetRatios_[i] - bestRatios[i].calculatedRatio_);
            for (int r = 0; r < networkResistorCount; ++r)
                std::cout << ' ' << (bestRatios[i].resistorStates_.test(r) ? 1 : 0);

            std::cout << std::endl;
        }
    }

    virtual void run(const std::vector<Resistor>& startingResistorsVect) final
    {
        if (startingResistorsVect.size() != networkResistorCount)
            throw std::logic_error("Starting resistors size doesn't match model parameters");

        ResistorArray startingResistors;
        std::copy(startingResistorsVect.begin(), startingResistorsVect.end(), startingResistors.begin());

        double startingScore = checkAccuracy(startingResistors); 

        auto bestResistors = startingResistors;
        auto bestScore = startingScore;

        auto currentResistors = bestResistors;
        double currentScore = bestScore;
        unsigned randomizations = 0;

        auto lastPrint = std::chrono::steady_clock::now();

        auto printProgress = [&]()
        {
            std::cout << FORMAT("[{}] Score: {:8.4f} ({:+8.4f}) best: {:8.4f}",
                randomizations,
                currentScore,
                currentScore - bestScore,
                bestScore);
            std::cout << "\n";
        };

        while (!exitedRequested())
        {
            auto [newResistors, newScore] = findBiggestChange(currentResistors, currentScore);

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
                std::bitset<networkResistorCount> switched{};
                while (randCount > 0)
                {
                    auto idx = resistorIndexDist_(rnd_);
                    if (!switched.test(idx))
                    {
                        currentResistors[idx] = resistor(resistorDomainIndexDist_(rnd_));
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

        auto bestRatios = getBestRatios(bestResistors);
        bestScore = checkAccuracy(bestResistors);

        auto percentDiffs = targetRatios_;
        auto pdMean = 0.0;
        for (unsigned i = 0; i < percentDiffs.size(); ++i)
        {
            percentDiffs[i] -= bestRatios[i].calculatedRatio_;
            percentDiffs[i] = 100.0 * std::abs(percentDiffs[i]) / targetRatios_[i];
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

        TargetRatios bestRatioValues;
        for (unsigned i = 0; i < bestRatioValues.size(); ++i)
            bestRatioValues[i] = bestRatios[i].calculatedRatio_;

        std::cout << "Targets: " << targetRatios_ << std::endl;
        std::cout << "Best:    " << bestRatioValues << std::endl;
        std::cout << "% diff:  " << percentDiffs << std::endl;
        std::cout << FORMAT("Score: {:10.4f} ({:10.4f})  stddev: {:.3g}\n",
            bestScore,
            bestScore - startingScore,
            pdStdDev);
        std::cout << bestResistors << std::endl;

        printTable(bestRatios, bestResistors);
    }
};

template<typename ResistorSeriesType>
struct ModelParameters
{
    unsigned networkResistorCount_;
    unsigned searchDistance_;
    unsigned targetCount_;
    
    using ResistorSeries = ResistorSeriesType;
};

template<
    unsigned networkResistorCount,
    unsigned searchDistance,
    std::size_t targetCount, 
    typename ResistorSeriesType>
std::unique_ptr<IModel> makeModel(
    const Potentials<targetCount>& targetPotentials,
    const ResistorSeriesType& resistorSeries,
    RandomGen& gen)
{
    return std::make_unique<Model<networkResistorCount, searchDistance, targetCount, ResistorSeriesType>>(targetPotentials, resistorSeries, gen);
}



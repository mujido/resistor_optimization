#pragma once

#include <limits>
#include <iostream>

class Resistor
{
public:
	constexpr Resistor() : 
		index_(std::numeric_limits<std::size_t>::max()), 
		value_(std::numeric_limits<double>::quiet_NaN()),
		reciprocal_(std::numeric_limits<double>::quiet_NaN())
	{}

	Resistor(std::size_t index, double value) :
		index_(index),
		value_(value),
		reciprocal_(1 / value_)
	{}

	std::size_t getIndex() const { return index_; }

	double getValue() const { return value_; }

	double getReciprocal() const
	{
		return reciprocal_;
	}

private:
	std::size_t index_;
	double value_;
	double reciprocal_;
};

inline std::ostream& operator<< (std::ostream& os, const Resistor& r)
{
	return os << r.getValue();
}
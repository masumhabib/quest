/* 
 * File:   random.h
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 6, 2014, 10:30 AM
 */

#ifndef RANDOM_H
#define RANDOM_H


#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <vector>
#include "maths/constants.h"
#include "utils/stringutils.h"
#include "utils/Printable.hpp"
#include "utils/serialize.hpp"
#include "utils/vout.h"
#include "utils/std.hpp"
#include "maths/arma.hpp"
#include "maths/linspace.hpp"



using utils::Printable;
using utils::strings::trim;
using namespace utils::stds;
using namespace maths::armadillo;
using maths::armadillo::vec;
using maths::linspace;
using maths::constants::pi;

namespace utils{
namespace random{

using std::vector;
using maths::constants::NaN;
namespace rnd = boost::random;

class Random {
public: 
    typedef std::unique_ptr<Random> Uptr;

public:
    virtual void reset() = 0;
    virtual double generate () = 0;
    virtual std::vector<double> generate (size_t count);
    virtual std::vector<double> generate (int count);

    virtual ~Random() = default;

protected:
    Random();

protected:
    rnd::random_device device;
    rnd::mt19937_64 generator;
};

class UniformRandom : public Random {
public:
    UniformRandom (double min, double max);

    void reset();

    using Random::generate;
    double generate ();

private:
    double min;
    double max;
    rnd::uniform_real_distribution<double> distribution;
};

class NormalRandom : public Random {
public:
    NormalRandom(double std_dev = 1.0, double mean = 0.0);
    NormalRandom(double std_dev, double mean, double min, double max);

    void reset();
    double generate();
    using Random::generate;

private:
    double min = NaN;
    double max = NaN;
    rnd::normal_distribution<double> distribution;
    static constexpr size_t MAX_RETRY = 1000;
};

typedef NormalRandom GaussianRandom;

class DiscreteRandom : public Random {
public:
    struct WeightFunction {
        virtual double operator () (const double& value) const = 0;
    };

public:
    DiscreteRandom(const std::vector<double>& values, 
            const std::vector<double>& weights);
    DiscreteRandom(const std::vector<double>& values, 
            const WeightFunction& weight_function);
    DiscreteRandom(double xmin, double xmax, std::size_t count, 
            const WeightFunction& weight_function);
    DiscreteRandom(double xmin, double xmax, double step, 
            const WeightFunction& weight_function);

    void reset();
    double generate ();
    using Random::generate;

private:
    static std::vector<double> calculate_weights (
            const std::vector<double>& domain, 
            const WeightFunction& weight_function);
    std::vector<double> create_domain (double xmin, double xmax, size_t count);
    std::vector<double> create_domain (double xmin, double xmax, double step);

private:
    std::vector<double> domain;
    rnd::discrete_distribution<unsigned int> distribution;
};

class DiscreteCosineRandom : public DiscreteRandom {
public:
    DiscreteCosineRandom(const std::vector<double>& angles);
    DiscreteCosineRandom(double angle_min, double angle_max, std::size_t count);
    DiscreteCosineRandom(double angle_min, double angle_max, double step);

public:
    struct CosineFunction : public WeightFunction {
        double operator () (const double& angle) const {
            return std::abs (std::cos (angle));
        }
    };
};


}}

#endif  /* RANDOM_H */


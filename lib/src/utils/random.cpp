/* 
 * File:   random.cpp
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 6, 2014, 10:30 AM
 */

#include "utils/random.h"

namespace utils{ namespace random{

Random::Random() : generator(device) { 
}

std::vector<double> Random::generate (int count) {
    if (count < 0) {
        throw std::invalid_argument("Random::generate(): count must be >= 0");
    }
    
    return generate((size_t) count);
}

std::vector<double> Random::generate (size_t count) {
    vector<double> numbers(count);
    for (size_t i = 0; i < count; i += 1){
        numbers[i] = generate();
    }
    return numbers;
}

//-----------------------------------------------------------------------------
// Uniform Random
//-----------------------------------------------------------------------------

UniformRandom::UniformRandom (double min, double max) 
: max(max), min(min), distribution(min, max) { 
}


double UniformRandom::generate () {
    return distribution(generator);
}

void UniformRandom::reset() {
    distribution.reset();
}

//-----------------------------------------------------------------------------
// Normal Random
//-----------------------------------------------------------------------------

NormalRandom::NormalRandom(double std_dev, double mean)
: distribution (mean, std_dev) {
}

NormalRandom::NormalRandom(double std_dev, double mean, double min, double max) 
: min (min), max (max), distribution (mean, std_dev) {
}

void NormalRandom::reset() {
    distribution.reset();
}

double NormalRandom::generate() {
    double random_number;

    for (size_t retry = 0; retry < MAX_RETRY; retry += 1) {
        random_number = distribution(generator);
        // check the range
        if ((std::isnan(min) || random_number >= min)
        &&  (std::isnan(max) || random_number <= max)) {
            break;
        }
    }

    return random_number;
}

//-----------------------------------------------------------------------------
// Discrete Random
//-----------------------------------------------------------------------------

DiscreteRandom::DiscreteRandom (const std::vector<double>& values, 
const std::vector<double>& weights) : domain (values), distribution (weights) {
    if (values.size() != weights.size()) {
        throw invalid_argument ("DiscreteRandom::DiscreteRandom(): values.size() != weights.size()");
    }
}

DiscreteRandom::DiscreteRandom(const std::vector<double>& values, 
const WeightFunction& weight_function): domain (values), distribution (
calculate_weights(values, weight_function)) {
}

DiscreteRandom::DiscreteRandom(double xmin, double xmax, std::size_t count, 
const WeightFunction& weight_function) : domain (create_domain (xmin, xmax, 
count)), distribution (calculate_weights (domain, weight_function)) {
}

DiscreteRandom::DiscreteRandom(double xmin, double xmax, double step, 
const WeightFunction& weight_function) : domain (create_domain (xmin, xmax, 
step)), distribution (calculate_weights (domain, weight_function)) {
}

void DiscreteRandom::reset() {
    distribution.reset();
}

double DiscreteRandom::generate () {
    int random_number = distribution(generator);
    assert (random_number < domain.size());
    return domain[random_number];
}

std::vector<double> DiscreteRandom::calculate_weights (
const std::vector<double>& domain, const WeightFunction& weight_function) {
    std::vector<double> weights;
    // for performance boost
    weights.reserve(domain.size());

    for (auto x : domain) {
        weights.push_back(weight_function(x));
    }

    return weights;
}

std::vector<double> DiscreteRandom::create_domain (
double xmin, double xmax, size_t count) {
    std::vector<double> new_domain = utils::stds::linspace(xmin, xmax, count);

    return new_domain;
}

std::vector<double> DiscreteRandom::create_domain (
double xmin, double xmax, double step) {
    std::vector<double> new_domain = utils::stds::linspace(xmin, xmax, step);

    return new_domain;
}

//-----------------------------------------------------------------------------
// Discrete Cosine Random
//-----------------------------------------------------------------------------
//

DiscreteCosineRandom::DiscreteCosineRandom (const std::vector<double>& angles) :
DiscreteRandom (angles, CosineFunction()) {
}

DiscreteCosineRandom::DiscreteCosineRandom(double angle_min, double angle_max, 
std::size_t count) : DiscreteRandom (angle_min, angle_max, count, 
CosineFunction ()) {
}

DiscreteCosineRandom::DiscreteCosineRandom(double angle_min, double angle_max, 
double step) : DiscreteRandom (angle_min, angle_max, step, CosineFunction()) {
}


}}

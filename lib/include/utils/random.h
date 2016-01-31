/* 
 * File:   random.h
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 6, 2014, 10:30 AM
 */

#ifndef RANDOM_H
#define	RANDOM_H

#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <vector>

namespace utils{namespace random{
using std::vector;
namespace br = boost::random;
void genNormalDist(vector<double> &result, double sigma = 1, double mean = 0);   
double getGaussianRand(double sigma, double mean, double min, double max, 
        bool reset = false);
double getGaussianRand(double sigma, double mean, bool reset = false);
double getUniformRand(double min, double max, bool reset = false);

}}

#endif	/* RANDOM_H */


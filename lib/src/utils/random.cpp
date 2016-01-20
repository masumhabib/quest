/* 
 * File:   random.cpp
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 6, 2014, 10:30 AM
 */

#include "utils/random.h"

namespace utils{namespace random{

void genNormalDist(vector<double> &result, double sigma, double mean){
    br::random_device device;
    br::normal_distribution<> distribution(mean, sigma);
    
    br::variate_generator<br::random_device&, 
                           br::normal_distribution<> > generator(device, distribution);
    
    for (auto it = result.begin(); it != result.end(); ++it){
        *it = generator();
    }
}

double getGaussianRand(double sigma, double mean, double min, double max, 
        bool reset){
    double number;
    do {
       number = getGaussianRand(sigma, mean, reset);
    } while (number < min || number > max); 

    return number;
}

double getGaussianRand(double sigma, double mean, bool reset){
    static br::random_device device;
    static br::normal_distribution<> distribution(mean, sigma);
    
    static br::variate_generator<br::random_device&, 
        br::normal_distribution<> > generator(device, distribution);
 
    if (reset) {
        distribution.reset();
    }
 
    return generator();
}

double getUniformRand(double min, double max, bool reset){
    static br::random_device device;
    static br::uniform_real_distribution<> distribution(min, max);

    if (reset) {
        distribution.reset();
    }
 
    return distribution(device);
}


}}

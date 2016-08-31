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

////!< Constructor - Gaussian Random (boost) or Uniform Random or Normal Dist
//Distribution::Distribution( DistributionType distributionType, double sigma_or_min,
//		double mean_or_max) : mNormal_dist_engn( mean_or_max, sigma_or_min ),
//			mUniform_dist_engn(  sigma_or_min, mean_or_max ),
//			mGenerator(mDevice, mNormal_dist_engn){
//	mDistributionType = distributionType;
//	if( distributionType==DistributionType::GAUSSIAN_RANDOM ||
//			distributionType==DistributionType::NORMAL){
//		mSigma = sigma_or_min;
//		mMean = mean_or_max;
//	}else if( distributionType==DistributionType::UNIFORM_RANDOM ){
//		mMin = sigma_or_min;
//		mMax = mean_or_max;
//	}else{
//		throw invalid_argument("No constructor for this type of Distribution");
//	}
//}
////!< Constructor - Gaussian Random (custom)
//Distribution::Distribution( DistributionType distributionType, double sigma,
//							double mean, double min, double max):
//							mNormal_dist_engn( mean, sigma ),
//							mUniform_dist_engn( min, max),
//							mGenerator(mDevice, mNormal_dist_engn){
//	mDistributionType = distributionType;
//	if( distributionType==DistributionType::GAUSSIAN_RANDOM ){
//		mSigma = sigma;
//		mMean = mean;
//		mMin = min;
//		mMax = max;
//	}else{
//		throw invalid_argument("No constructor for this type of Distribution");
//	}
//}
////!< Constructor - Cosine Distribution
//Distribution::Distribution( DistributionType distributionType, double min,
//		double max, std::vector<double> angles, std::vector<double> weights):
//									mNormal_dist_engn( 0.0, 0.0 ),
//									mUniform_dist_engn( 0.0, 0.0),
//									mGenerator(mDevice, mNormal_dist_engn),
//									mDiscrt_dist_engn(weights){
//	mDistributionType = distributionType;
//	rnd::discrete_distribution<>* tempDiscrtEngn = new rnd::discrete_distribution<>(weights);
//	if( distributionType==DistributionType::DISCRETE ){
//		mMin = min;
//		mMax = max;
//		mDiscrtAngle = angles;
//		mDiscrtWeights = weights;
//	}else{
//		throw invalid_argument("No constructor for this type of Distribution");
//	}
//}
//
////!< Constructor - No Distribution, returns fixed value
//Distribution::Distribution( DistributionType distributionType, double fixed):
//							mNormal_dist_engn( 0.0, 0.0 ),
//							mUniform_dist_engn( 0.0, 0.0),
//							mGenerator(mDevice, mNormal_dist_engn){
//	mDistributionType = distributionType;
//	if( distributionType==DistributionType::NONE ){
//		mMean = fixed;
//	}else{
//		throw invalid_argument("No constructor for this type of Distribution");
//	}
//}
////!< Normal Dist
//void Distribution::getDistribution( vector<double> &result ){
//	if( mDistributionType == DistributionType::NORMAL ){
//		genNormalDist(result);
//	}else{
//		throw invalid_argument("Wrong Distribution Type");
//	}
//}
////!< Gaussian Random custom, Gaussian Random boost
//double Distribution::getDistribution( bool reset ){
//	double number = NaN;
//	if( mDistributionType == DistributionType::UNIFORM_RANDOM ){
//		number = getUniformRand(reset);
//	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM ){
//		number = getGaussianRand(reset);
//	}else if( mDistributionType == DistributionType::DISCRETE ){
//		number = getDiscrtRand();
//	}else if( mDistributionType == DistributionType::NONE ){
//		number = mMean;
//	}else{
//		throw invalid_argument("Wrong Distribution Type");
//	}
//	return number;
//}
//void Distribution::genNormalDist(vector<double> &result){
//    for (auto it = result.begin(); it != result.end(); ++it){
//        *it = mGenerator();
//    }
//}
//
//double Distribution::getGaussianRand(bool reset){
//	if (reset) {
//	    	mNormal_dist_engn.reset();
//	}
//	double number = NaN;
//	// is_nan checking
//	if ( mMin != mMin && mMax != mMax ){   // without Min and Max Feature
//		number = mGenerator();
//	}else{
//		do {
//		   number = mGenerator();
//		} while (number < mMin || number > mMax);
//	}
//    return number;
//}
//
//double Distribution::getUniformRand(bool reset){
//    if (reset) {
//    	mUniform_dist_engn.reset();
//    }
//    return mUniform_dist_engn( mDevice );
//}
//
//double Distribution::getDiscrtRand(){
//	double number = NaN;
//	do {
//		number = mDiscrtAngle[ mDiscrt_dist_engn(mDevice) ];
//	} while (number < mMin || number > mMax);
//	return number;
//}
//
//string Distribution::toString() const{
//
//	stringstream out;
//	out << " Distribution: " << endl;
//	out.precision(4);
//	out << std::fixed;
//
//	string distype;
//	if( mDistributionType == DistributionType::NONE ){
//			distype = "no distribution";
//	}else if( mDistributionType == DistributionType::UNIFORM_RANDOM ){
//		distype = "uniform random";
//	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM && mMin!=mMin && mMax != mMax){
//		distype = "gaussian random without minm and maxm";
//	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM && mMin==mMin && mMax==mMax){
//		distype = "gaussian random with minm and maxm";
//	}else if( mDistributionType == DistributionType::NORMAL){
//		distype = "normal";
//	}else if( mDistributionType == DistributionType::DISCRETE){
//		distype = "discrete";
//	}
//	out << " type: " << std::fixed  << distype << endl;
//	out << " Sigma or Variance: " << mSigma << endl;
//	out << " Mean: " << mMean << endl;
//	out << " Minimum: " << mMin << endl;
//	out << " Maximum: " << mMax << endl;
//
//	return out.str();
//}


}}

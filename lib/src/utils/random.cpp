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

UniformRandom::UniformRandom (double min, double max) 
: max(max), min(min), distribution(min, max) { 
}


double UniformRandom::generate () {
    return distribution(generator);
}

void UniformRandom::reset() {
    distribution.reset();
}


//!< Constructor - Gaussian Random (boost) or Uniform Random or Normal Dist
Distribution::Distribution( DistributionType distributionType, double sigma_or_min,
		double mean_or_max) : mNormal_dist_engn( mean_or_max, sigma_or_min ),
			mUniform_dist_engn(  sigma_or_min, mean_or_max ),
			mGenerator(mDevice, mNormal_dist_engn){
	mDistributionType = distributionType;
	if( distributionType==DistributionType::GAUSSIAN_RANDOM ||
			distributionType==DistributionType::NORMAL){
		mSigma = sigma_or_min;
		mMean = mean_or_max;
	}else if( distributionType==DistributionType::UNIFORM_RANDOM ){
		mMin = sigma_or_min;
		mMax = mean_or_max;
	}else{
		throw invalid_argument("No constructor for this type of Distribution");
	}
}
//!< Constructor - Gaussian Random (custom)
Distribution::Distribution( DistributionType distributionType, double sigma,
							double mean, double min, double max):
							mNormal_dist_engn( mean, sigma ),
							mUniform_dist_engn( min, max),
							mGenerator(mDevice, mNormal_dist_engn){
	mDistributionType = distributionType;
	if( distributionType==DistributionType::GAUSSIAN_RANDOM ){
		mSigma = sigma;
		mMean = mean;
		mMin = min;
		mMax = max;
	}else{
		throw invalid_argument("No constructor for this type of Distribution");
	}
}
//!< Constructor - Cosine Distribution
Distribution::Distribution( DistributionType distributionType, double min,
		double max, std::vector<double> angles, std::vector<double> weights):
									mNormal_dist_engn( 0.0, 0.0 ),
									mUniform_dist_engn( 0.0, 0.0),
									mGenerator(mDevice, mNormal_dist_engn),
									mDiscrt_dist_engn(weights){
	mDistributionType = distributionType;
	rnd::discrete_distribution<>* tempDiscrtEngn = new rnd::discrete_distribution<>(weights);
	if( distributionType==DistributionType::DISCRETE ){
		mMin = min;
		mMax = max;
		mDiscrtAngle = angles;
		mDiscrtWeights = weights;
	}else{
		throw invalid_argument("No constructor for this type of Distribution");
	}
}

//!< Constructor - No Distribution, returns fixed value
Distribution::Distribution( DistributionType distributionType, double fixed):
							mNormal_dist_engn( 0.0, 0.0 ),
							mUniform_dist_engn( 0.0, 0.0),
							mGenerator(mDevice, mNormal_dist_engn){
	mDistributionType = distributionType;
	if( distributionType==DistributionType::NONE ){
		mMean = fixed;
	}else{
		throw invalid_argument("No constructor for this type of Distribution");
	}
}
//!< Normal Dist
void Distribution::getDistribution( vector<double> &result ){
	if( mDistributionType == DistributionType::NORMAL ){
		genNormalDist(result);
	}else{
		throw invalid_argument("Wrong Distribution Type");
	}
}
//!< Gaussian Random custom, Gaussian Random boost
double Distribution::getDistribution( bool reset ){
	double number = NaN;
	if( mDistributionType == DistributionType::UNIFORM_RANDOM ){
		number = getUniformRand(reset);
	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM ){
		number = getGaussianRand(reset);
	}else if( mDistributionType == DistributionType::DISCRETE ){
		number = getDiscrtRand();
	}else if( mDistributionType == DistributionType::NONE ){
		number = mMean;
	}else{
		throw invalid_argument("Wrong Distribution Type");
	}
	return number;
}
void Distribution::genNormalDist(vector<double> &result){
    for (auto it = result.begin(); it != result.end(); ++it){
        *it = mGenerator();
    }
}

double Distribution::getGaussianRand(bool reset){
	if (reset) {
	    	mNormal_dist_engn.reset();
	}
	double number = NaN;
	// is_nan checking
	if ( mMin != mMin && mMax != mMax ){   // without Min and Max Feature
		number = mGenerator();
	}else{
		do {
		   number = mGenerator();
		} while (number < mMin || number > mMax);
	}
    return number;
}

double Distribution::getUniformRand(bool reset){
    if (reset) {
    	mUniform_dist_engn.reset();
    }
    return mUniform_dist_engn( mDevice );
}

double Distribution::getDiscrtRand(){
	double number = NaN;
	do {
		number = mDiscrtAngle[ mDiscrt_dist_engn(mDevice) ];
	} while (number < mMin || number > mMax);
	return number;
}

string Distribution::toString() const{

	stringstream out;
	out << " Distribution: " << endl;
	out.precision(4);
	out << std::fixed;

	string distype;
	if( mDistributionType == DistributionType::NONE ){
			distype = "no distribution";
	}else if( mDistributionType == DistributionType::UNIFORM_RANDOM ){
		distype = "uniform random";
	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM && mMin!=mMin && mMax != mMax){
		distype = "gaussian random without minm and maxm";
	}else if( mDistributionType == DistributionType::GAUSSIAN_RANDOM && mMin==mMin && mMax==mMax){
		distype = "gaussian random with minm and maxm";
	}else if( mDistributionType == DistributionType::NORMAL){
		distype = "normal";
	}else if( mDistributionType == DistributionType::DISCRETE){
		distype = "discrete";
	}
	out << " type: " << std::fixed  << distype << endl;
	out << " Sigma or Variance: " << mSigma << endl;
	out << " Mean: " << mMean << endl;
	out << " Minimum: " << mMin << endl;
	out << " Maximum: " << mMax << endl;

	return out.str();
}


}}

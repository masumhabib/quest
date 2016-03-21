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
#include "maths/constants.h"
#include "utils/stringutils.h"
#include "utils/Printable.hpp"
#include "utils/serialize.hpp"
#include "utils/vout.h"
#include "utils/std.hpp"



using utils::Printable;
using utils::strings::trim;
using namespace utils::stds;
namespace utils{namespace random{
using std::vector;
using maths::constants::NaN;
namespace br = boost::random;
class Distribution {
public:
	//!< Constructor - Gaussian Random (boost) or Uniform Random or Normal Dist
	Distribution( int distributionType, double sigma_or_min,
			double mean_or_max);
	//!< Constructor - Gaussian Random (custom)
	Distribution( int distributionType, double sigma, double mean, double min,
			double max);
	//!< Normal Dist
	void getDistribution( vector<double> &result );
	//!< Gaussian Random custom, Gaussian Random boost, Uniform Random
	double getDistribution( bool reset=false );
	string toString() const;

public:
	// distribution types
    static constexpr int NORMAL= 1;
    static constexpr int GAUSSIAN_RANDOM = 2;
    static constexpr int UNIFORM_RANDOM = 3;

private:
    int mDistributionType;
    double mSigma;
    double mMean;
    double mMax = NaN;
    double mMin = NaN;

    //!< generator and engines
	br::random_device mDevice;
	br::normal_distribution<> mNormal_dist_engn;
	br::uniform_real_distribution<> mUniform_dist_engn;
	br::variate_generator< br::random_device&,
		br::normal_distribution<> > mGenerator;

    void genNormalDist(vector<double> &result);
    double getGaussianRand(bool reset = false);
    double getUniformRand(bool reset = false);
};
}}

#endif	/* RANDOM_H */


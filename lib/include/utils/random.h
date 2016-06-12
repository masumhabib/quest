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
namespace br = boost::random;

enum class DistributionType{
	NONE = 0,
	NORMAL = 1,
	GAUSSIAN_RANDOM = 2,
	UNIFORM_RANDOM = 3,
	DISCRETE = 4
};
class Distribution {
public:
	typedef shared_ptr<Distribution> ptr;
	//!< Constructor - Gaussian Random (boost) or Uniform Random or Normal Dist
	Distribution( DistributionType distributionType, double sigma_or_min,
			double mean_or_max);
	//!< Constructor - Cosine
	Distribution( DistributionType distributionType, double sigma_or_min,
			double mean_or_max, double deltaAngle);
	//!< Constructor - Gaussian Random (custom)
	Distribution( DistributionType distributionType, double sigma, double mean,
			double min, double max);
	//!< Constructor - Cosine
	Distribution( DistributionType distributionType, double min, double max,
			std::vector<double> angles, std::vector<double> weights);
	//!< Constructor - No Distribution, returns constant value
	Distribution( DistributionType distributionType, double fixed=0.0 );
	//!< Normal Dist
	void getDistribution( vector<double> &result );
	//!< Gaussian Random custom, Gaussian Random boost, Uniform Random
	double getDistribution( bool reset=false );
	string toString() const;

private:
    DistributionType mDistributionType;
    double mSigma;
    double mMean;
    double mMax = NaN;
    double mMin = NaN;

    //!< generator and engines
	br::random_device mDevice;
	br::normal_distribution<> mNormal_dist_engn; //!< Normal Engine
	br::uniform_real_distribution<> mUniform_dist_engn; //!< Uniform Engine
	br::variate_generator< br::random_device&,
	br::normal_distribution<> > mGenerator;
	br::discrete_distribution<> mDiscrt_dist_engn; //!< Discrete Engine

    void genNormalDist(vector<double> &result);
    double getGaussianRand(bool reset = false);
    double getUniformRand(bool reset = false);
    double getDiscrtRand();

    vec mDiscrtAngle; //!< Discrete Values
    vec mDiscrtWeights ; //!< Discrete Weights
};
}}

#endif	/* RANDOM_H */


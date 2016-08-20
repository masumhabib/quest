/** Test cases for Segment class.
 *
 */

#include "utils/random.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RandomTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>
#include <map>
#include <string>
#include <iomanip>

using namespace utils::random;
using namespace std;

BOOST_AUTO_TEST_CASE(UniformRandom_generator)
{
	int N = 1000000;
	double min = -8, max = 12, span;
	span = max - min;
	double expected_freq = N/span;

	map<int, int> hist;
    for (int i = min; i <= max; i += 1) {
        hist[i] = 0;
    }

	UniformRandom uniform_random (min, max);
	for (int n = 0; n < N; n += 1) {
	    hist[std::round(uniform_random.generate())] += 1;
    }

    double freq_error = 0;
    for (int i = min + 1; i < max; i += 1) {
        cout << "freq[" << i << "] = " << hist[i] << endl;
        freq_error += std::abs((hist[i] - expected_freq)/expected_freq*100);
    }
    double avg_freq_error = freq_error/span;

    cout << "freq_error = " << avg_freq_error << endl;
    BOOST_CHECK(avg_freq_error < 1);

    vector <int> keys;
    for (auto item : hist) {
        keys.push_back(item.first);
    }
    double min_obtained = *std::min_element(keys.begin(), keys.end());
    BOOST_CHECK_EQUAL(min_obtained, min);
    double max_obtained = *std::max_element(keys.begin(), keys.end());
    BOOST_CHECK_EQUAL(max_obtained, max);
}

//bool print_histogram = true;
//
//void printBannerTestCase(string testName)
//{
//	int N = 50;
//	cout << "\n";
//	for (int i = 1; i<=N; i++){
//		cout << "*" ;
//	}
//	cout << "\n" ;
//	int dashNo = ( N - testName.length() - 2 )  /  2;
//	for (int i = 1; i<=dashNo; i++){
//		cout << "-" ;
//	}
//	cout << " " << testName ;
//
//	for (int i = 1; i <= N-2*dashNo-testName.length()-1; i++){
//			cout << " " ;
//	}
//	for (int i = 1; i<=dashNo; i++){
//			cout << "-" ;
//	}
//	cout << "\n";
//	for (int i = 1; i<=N; i++){
//			cout << "*" ;
//	}
//	cout << "\n\n";
//}
//
//BOOST_AUTO_TEST_CASE(Generate_Normal_Distribution)
//{
//	printBannerTestCase( "Testing Normal Distribution" );
//    int N = 100;
//    double expmean = 0;
//    double sigma = 1;
//    vector<double> vs(N);
//	Distribution normaldist(DistributionType::NORMAL, sigma, expmean);
//	std::cout << normaldist.toString() << endl;
//	normaldist.getDistribution(vs);
//    //genNormalDist(vs, sigma, expmean);
//    double foundMean = 0;
//    for (auto v : vs) {
//    	foundMean += v;
//    }
//    foundMean = foundMean/N;
//    BOOST_CHECK(abs(foundMean - expmean ) < 0.1*abs(foundMean-sigma));
//    cout << "Expected Mean: " << expmean << endl;
//    cout << "Generated Mean: " << foundMean << endl;
//}
//
//BOOST_AUTO_TEST_CASE(Gaussian_Random_Distribution)
//{
//	printBannerTestCase( "Testing Gaussian Random Distribution" );
//	int N = 1000000;
//	std::map<int, int> hist;
//	double sigma, expmean, min, max, foundMean;
//	sigma = 2;
//	max = 12;
//	min = -2;
//	expmean = 5;
//	foundMean = 0;
//	bool reset = false;
//	Distribution gaussianRandDist( DistributionType::GAUSSIAN_RANDOM, sigma, expmean, min, max );
//	std::cout << gaussianRandDist.toString() << endl;
//	gaussianRandDist.getDistribution(true);
//	for( int n=0; n<N; n++ ){
//		double number = gaussianRandDist.getDistribution(reset);
//		foundMean += number;
//		++hist[std::round( number)];
//	}
//	foundMean = foundMean/double(N);
//	//checking mean
//	if ( expmean!=0 ){
//		BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//	}
//	BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*sigma);
//	//checking minimum
//	BOOST_CHECK(hist.begin()->first >= min);
//	//checking maximum
//	BOOST_CHECK(hist.rbegin()->first <= max );
//	if( print_histogram ){
//		cout << "Histogram : " << endl;
//		for(auto p : hist) {
//				std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//						  << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//		}
//		cout << "\n";
//	}
//	cout << "Expected Mean: " << expmean << endl;
//	cout << "Generated Mean: " << foundMean << endl;
//}
//
//BOOST_AUTO_TEST_CASE(Uniform_Random_Distribution)
//{
//	printBannerTestCase( "Testing Uniform Random Distribution" );
//	int N = 1000000;
//	std::map<int, int> hist;
//	double expmean, min, max, foundMean, dSigma;
//	min = -2;
//	max = 12;
//	expmean = (min + max)/2;
//
//	foundMean = 0;
//	bool reset = false;
//	Distribution uniformRandDist( DistributionType::UNIFORM_RANDOM, min, max);
//	std::cout << uniformRandDist.toString() << endl;
//	uniformRandDist.getDistribution(true);
//	for( int n=0; n<N; n++ ){
//		double number = uniformRandDist.getDistribution(reset);
//		foundMean += number;
//		++hist[std::round( number)];
//	}
//	foundMean = foundMean/double(N);
//	//checking mean
//	if ( expmean!=0 ){
//			BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//	}
//	BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*abs(abs(min)-abs(expmean)));
//	//checking minimum
//	BOOST_CHECK(hist.begin()->first >= min);
//	//checking maximum
//	BOOST_CHECK(hist.rbegin()->first <= max );
//	if( print_histogram ){
//		cout << "Histogram : " << endl;
//		for(auto p : hist) {
//				std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//						  << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//		}
//		cout << "\n";
//	}
//	cout << "Expected Mean: " << expmean << endl;
//	cout << "Generated Mean: " << foundMean << endl;
//}
//
//
//BOOST_AUTO_TEST_CASE(Cosine_Distribution)
//{
//	printBannerTestCase( "Testing Cosine Random Distribution" );
//	int N = 1000000;
//	std::map<int, int> hist;
//	double expmean, min, max, foundMean;
//	min = -pi/2;
//	max = pi/2;
//	expmean = (min + max)/2;
//	vec angles = linspace(min, max, (max - min)/100.0);
//	vec tempweights = cos( angles );
//	vec weights = tempweights/ accu( tempweights );
//	vector<double> anglesStdVec = arma::conv_to< vector<double> >::from( angles );
//	vector<double> weightsStdVec = arma::conv_to< vector<double> >::from( weights );
//	foundMean = 0;
//	bool reset = false;
//	Distribution cosineRandDist( DistributionType::DISCRETE, min, max, anglesStdVec, weightsStdVec);
//	std::cout << cosineRandDist.toString() << endl;
//	cosineRandDist.getDistribution(true);
//	for( int n=0; n<N; n++ ){
//		double number = cosineRandDist.getDistribution(reset);
//		foundMean += number;
//		++hist[std::round( number)];
//	}
//	foundMean = foundMean/double(N);
//	//checking mean
//	if ( expmean!=0 ){
//			BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//	}
//	BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*abs(abs(min)-abs(expmean)));
//	//checking minimum
//	BOOST_CHECK(hist.begin()->first >= min);
//	//checking maximum
//	BOOST_CHECK(hist.rbegin()->first <= max );
//	if( print_histogram ){
//		cout << "Histogram : " << endl;
//		for(auto p : hist) {
//				std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//						  << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//		}
//		cout << "\n";
//	}
//	cout << "Expected Mean: " << expmean << endl;
//	cout << "Generated Mean: " << foundMean << endl;
//}





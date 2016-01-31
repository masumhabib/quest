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

bool print_histogram = true;

void printBannerTestCase(string testName)
{
	int N = 50;
	cout << "\n";
	for (int i = 1; i<=N; i++){
		cout << "*" ;
	}
	cout << "\n" ;
	int dashNo = ( N - testName.length() - 2 )  /  2;
	for (int i = 1; i<=dashNo; i++){
		cout << "-" ;
	}
	cout << " " << testName ;

	for (int i = 1; i <= N-2*dashNo-testName.length()-1; i++){
			cout << " " ;
	}
	for (int i = 1; i<=dashNo; i++){
			cout << "-" ;
	}
	cout << "\n";
	for (int i = 1; i<=N; i++){
			cout << "*" ;
	}
	cout << "\n\n";
}

BOOST_AUTO_TEST_CASE(Generate_Normal_Distribution)
{
	printBannerTestCase( "Testing Normal Distribution" );
    int N = 1000000;
    double expmean = 0;
    double sigma = 1;
    vector<double> vs(N);
    genNormalDist(vs, sigma, expmean);
    
    double foundMean = 0;
    for (auto v : vs) {
    	foundMean += v;
    }
    foundMean = foundMean/N;
    BOOST_CHECK(abs(foundMean - expmean ) < 0.1*abs(foundMean-sigma));
    cout << "Expected Mean: " << expmean << endl;
    cout << "Generated Mean: " << foundMean << endl;
}

BOOST_AUTO_TEST_CASE(Gaussian_Random_Distribution)
{
	printBannerTestCase( "Testing Gaussian Random Distribution" );
	int N = 1000000;
	std::map<int, int> hist;
	double sigma, expmean, min, max, foundMean;
	sigma = 2;
	max = 12;
	min = -2;
	expmean = 5;
	foundMean = 0;
	bool reset = false;

	for( int n=0; n<N; n++ ){
		double number = getGaussianRand(sigma, expmean, min, max, reset );
		foundMean += number;
		++hist[std::round( number)];
	}
	foundMean = foundMean/double(N);
	//checking mean
	if ( expmean!=0 ){
		BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
	}
	BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*sigma);
	//checking minimum
	BOOST_CHECK(hist.begin()->first >= min);
	//checking maximum
	BOOST_CHECK(hist.rbegin()->first <= max );
	if( print_histogram ){
		cout << "Histogram : " << endl;
		for(auto p : hist) {
				std::cout << std::fixed << std::setprecision(1) << std::setw(2)
						  << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
		}
		cout << "\n";
	}
	cout << "Expected Mean: " << expmean << endl;
	cout << "Generated Mean: " << foundMean << endl;
}

BOOST_AUTO_TEST_CASE(Uniform_Random_Distribution)
{
	printBannerTestCase( "Testing Uniform Random Distribution" );
	int N = 1000000;
	std::map<int, int> hist;
	double expmean, min, max, foundMean;
	min = -2;
	max = 12;
	expmean = (min + max)/2;
	foundMean = 0;
	bool reset = false;

	for( int n=0; n<N; n++ ){
		double number = getUniformRand(min, max, reset );
		foundMean += number;
		++hist[std::round( number)];
	}
	foundMean = foundMean/double(N);
	//checking mean
	if ( expmean!=0 ){
			BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
	}
	BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*abs(abs(min)-abs(expmean)));
	//checking minimum
	BOOST_CHECK(hist.begin()->first >= min);
	//checking maximum
	BOOST_CHECK(hist.rbegin()->first <= max );
	if( print_histogram ){
		cout << "Histogram : " << endl;
		for(auto p : hist) {
				std::cout << std::fixed << std::setprecision(1) << std::setw(2)
						  << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
		}
		cout << "\n";
	}
	cout << "Expected Mean: " << expmean << endl;
	cout << "Generated Mean: " << foundMean << endl;
}





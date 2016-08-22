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


class RandomAnalyzer {
public:
    RandomAnalyzer (const std::vector<double>& random_numbers) : 
    m_random_numbers(random_numbers) {
        analyze();
    }

public:
    double max () { return m_max; }
    double min () { return m_min; }
    double mean () { return m_mean; }
    double std_dev () { return m_std_dev; }
    size_t N () { return m_N; }

    double distribution_error(double expected_frequency) {
        std::vector<double> freqs (m_hist.size(), expected_frequency);
        return distribution_error(freqs);
    }

    double distribution_error(const std::vector<double>& expected_frequency) {
        assert (m_hist.size() == expected_frequency.size());
        double dist_error = 0;

        double span = m_max - m_min;
        for (int i = 1; i < m_hist.size() - 1; i += 1) {
            int x = i + m_min;
            cout << "freq[" << i << "] = " << m_hist[x] << endl;
            cout << "expected_freq[" << i << "] = " << expected_frequency[i] << endl;
            dist_error += std::abs((m_hist[x] - expected_frequency[i])
                          / expected_frequency[i]*100);
        }
        dist_error /= span;
        cout << "Dist error: " << dist_error << endl;

        return dist_error/span;
    }

private:
    void analyze() {
        m_N = m_random_numbers.size();
        m_max = utils::stds::max(m_random_numbers);
        m_min = utils::stds::min(m_random_numbers);
        m_mean = utils::stds::mean(m_random_numbers);
        m_std_dev = utils::stds::std_dev(m_random_numbers);

        for (int i = m_min; i <= m_max; i += 1) {
            m_hist[i] = 0;
        }

        for (int n = 0; n < m_N; n += 1) {
            m_hist[std::round(m_random_numbers[n])] += 1;
        }
    }

private:
    std::vector<double> m_random_numbers;
    size_t m_N = 0;
    double m_min = std::nan("NaN");
    double m_max = std::nan("NaN");
    double m_mean = std::nan("NaN");
    double m_std_dev = std::nan("NaN");
    std::map<int, int> m_hist;
};

BOOST_AUTO_TEST_CASE(UniformRandom_generates_uniform_distribution)
{
    int N = 1000000;
    double min = -8, max = 12, span;
    span = max - min;
    double expected_freq = N/span;

    UniformRandom rnd_generator (min, max);
    std::vector<double> random_numbers;
    random_numbers.reserve(N);
    for (int i = 0; i < N; i += 1) {
        random_numbers.push_back(rnd_generator.generate());
    }

    RandomAnalyzer analyzer (random_numbers);
    BOOST_CHECK_EQUAL(std::round(analyzer.min()), min);
    BOOST_CHECK_EQUAL(std::round(analyzer.max()), max);
    BOOST_CHECK(analyzer.distribution_error(expected_freq) < 1);
}

BOOST_AUTO_TEST_CASE(Random_returns_a_vector_of_numbers)
{
    size_t N = 1000000;
    double min = -8, max = 12, span;
    span = max - min;

    UniformRandom uniform_random (min, max);
    auto random_numbers = uniform_random.generate(N);
    RandomAnalyzer analyzer (random_numbers);

    BOOST_REQUIRE_EQUAL(random_numbers.size(), N);

    double expected_average = (max+min)/2;
    BOOST_CHECK_CLOSE(expected_average, analyzer.mean(), 1);
}

BOOST_AUTO_TEST_CASE(NormalRandom_generates_random_numbers_with_correct_stddev_and_mean)
{
    size_t N = 1000000;
    double std_dev = 10, mean = 5;

    NormalRandom random (std_dev, mean);
    auto random_numbers = random.generate(N);
    RandomAnalyzer analyzer (random_numbers);

    BOOST_REQUIRE_EQUAL(random_numbers.size(), N);

    BOOST_CHECK_CLOSE(mean, analyzer.mean(), 1);
    BOOST_CHECK_CLOSE(std_dev, analyzer.std_dev(), 1);
}

BOOST_AUTO_TEST_CASE(NormalRandom_generates_random_numbers_with_correct_min_max)
{
    size_t N = 1000000;
    double std_dev = 10, mean = 0;
    double min = -10, max = 20;

    NormalRandom random (std_dev, mean, min, max);
    auto random_numbers = random.generate(N);
    RandomAnalyzer analyzer (random_numbers);

    BOOST_REQUIRE_EQUAL(random_numbers.size(), N);

    BOOST_CHECK_CLOSE(analyzer.min(), min, 1);
    BOOST_CHECK_CLOSE(analyzer.max(), max, 1);
}

BOOST_AUTO_TEST_CASE(DiscreteRandom_generates_uniform_distribution_with_given_weights)
{
    int N = 1000000;
    double min = -8, max = 12, span;
    span = max - min;
    double expected_freq = N/span;

    std::vector<double> domain = utils::stds::linspace(min, max, size_t(span+1));
    std::vector<double> weights (span+1, expected_freq);

    DiscreteRandom rnd_generator (domain, weights);
    std::vector<double> random_numbers = rnd_generator.generate(N);

    RandomAnalyzer analyzer (random_numbers);
    BOOST_CHECK_EQUAL(std::round(analyzer.min()), min);
    BOOST_CHECK_EQUAL(std::round(analyzer.max()), max);
    BOOST_CHECK(analyzer.distribution_error(expected_freq) < 1);
}

BOOST_AUTO_TEST_CASE(DiscreteRandom_generates_uniform_distribution_with_given_weight_function)
{
    struct UniformWeight : public DiscreteRandom::WeightFunction {
        UniformWeight (double weight) : m_weight(weight) {}
        double operator () (const double& value) const {
            return m_weight;
        }

        double m_weight;
    };

    int N = 1000000;
    double min = -8, max = 12, span;
    span = max - min;
    double expected_freq = N/span;

    std::vector<double> domain = utils::stds::linspace(min, max, size_t(span+1));

    DiscreteRandom rnd_generator (domain, UniformWeight(expected_freq));
    std::vector<double> random_numbers = rnd_generator.generate(N);

    RandomAnalyzer analyzer (random_numbers);
    BOOST_CHECK_EQUAL(std::round(analyzer.min()), min);
    BOOST_CHECK_EQUAL(std::round(analyzer.max()), max);
    BOOST_CHECK(analyzer.distribution_error(expected_freq) < 1);
}

BOOST_AUTO_TEST_CASE(DiscreteRandom_generates_uniform_distribution_with_given_min_max_count)
{
    struct UniformWeight : public DiscreteRandom::WeightFunction {
        UniformWeight (double weight) : m_weight(weight) {}
        double operator () (const double& value) const {
            return m_weight;
        }

        double m_weight;
    };

    int N = 1000000;
    double min = -8, max = 12, span;
    span = max - min;
    double expected_freq = N/span;

    DiscreteRandom rnd_generator (size_t(span+1), min, max, UniformWeight(expected_freq));
    std::vector<double> random_numbers = rnd_generator.generate(N);

    RandomAnalyzer analyzer (random_numbers);
    BOOST_CHECK_EQUAL(std::round(analyzer.min()), min);
    BOOST_CHECK_EQUAL(std::round(analyzer.max()), max);
    BOOST_CHECK(analyzer.distribution_error(expected_freq) < 1);
}

//bool print_histogram = true;
//
//void printBannerTestCase(string testName)
//{
//  int N = 50;
//  cout << "\n";
//  for (int i = 1; i<=N; i++){
//      cout << "*" ;
//  }
//  cout << "\n" ;
//  int dashNo = ( N - testName.length() - 2 )  /  2;
//  for (int i = 1; i<=dashNo; i++){
//      cout << "-" ;
//  }
//  cout << " " << testName ;
//
//  for (int i = 1; i <= N-2*dashNo-testName.length()-1; i++){
//          cout << " " ;
//  }
//  for (int i = 1; i<=dashNo; i++){
//          cout << "-" ;
//  }
//  cout << "\n";
//  for (int i = 1; i<=N; i++){
//          cout << "*" ;
//  }
//  cout << "\n\n";
//}
//
//BOOST_AUTO_TEST_CASE(Generate_Normal_Distribution)
//{
//  printBannerTestCase( "Testing Normal Distribution" );
//    int N = 100;
//    double expmean = 0;
//    double sigma = 1;
//    vector<double> vs(N);
//  Distribution normaldist(DistributionType::NORMAL, sigma, expmean);
//  std::cout << normaldist.toString() << endl;
//  normaldist.getDistribution(vs);
//    //genNormalDist(vs, sigma, expmean);
//    double foundMean = 0;
//    for (auto v : vs) {
//      foundMean += v;
//    }
//    foundMean = foundMean/N;
//    BOOST_CHECK(abs(foundMean - expmean ) < 0.1*abs(foundMean-sigma));
//    cout << "Expected Mean: " << expmean << endl;
//    cout << "Generated Mean: " << foundMean << endl;
//}
//
//BOOST_AUTO_TEST_CASE(Gaussian_Random_Distribution)
//{
//  printBannerTestCase( "Testing Gaussian Random Distribution" );
//  int N = 1000000;
//  std::map<int, int> hist;
//  double sigma, expmean, min, max, foundMean;
//  sigma = 2;
//  max = 12;
//  min = -2;
//  expmean = 5;
//  foundMean = 0;
//  bool reset = false;
//  Distribution gaussianRandDist( DistributionType::GAUSSIAN_RANDOM, sigma, expmean, min, max );
//  std::cout << gaussianRandDist.toString() << endl;
//  gaussianRandDist.getDistribution(true);
//  for( int n=0; n<N; n++ ){
//      double number = gaussianRandDist.getDistribution(reset);
//      foundMean += number;
//      ++hist[std::round( number)];
//  }
//  foundMean = foundMean/double(N);
//  //checking mean
//  if ( expmean!=0 ){
//      BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//  }
//  BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*sigma);
//  //checking minimum
//  BOOST_CHECK(hist.begin()->first >= min);
//  //checking maximum
//  BOOST_CHECK(hist.rbegin()->first <= max );
//  if( print_histogram ){
//      cout << "Histogram : " << endl;
//      for(auto p : hist) {
//              std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//                        << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//      }
//      cout << "\n";
//  }
//  cout << "Expected Mean: " << expmean << endl;
//  cout << "Generated Mean: " << foundMean << endl;
//}
//
//BOOST_AUTO_TEST_CASE(Uniform_Random_Distribution)
//{
//  printBannerTestCase( "Testing Uniform Random Distribution" );
//  int N = 1000000;
//  std::map<int, int> hist;
//  double expmean, min, max, foundMean, dSigma;
//  min = -2;
//  max = 12;
//  expmean = (min + max)/2;
//
//  foundMean = 0;
//  bool reset = false;
//  Distribution uniformRandDist( DistributionType::UNIFORM_RANDOM, min, max);
//  std::cout << uniformRandDist.toString() << endl;
//  uniformRandDist.getDistribution(true);
//  for( int n=0; n<N; n++ ){
//      double number = uniformRandDist.getDistribution(reset);
//      foundMean += number;
//      ++hist[std::round( number)];
//  }
//  foundMean = foundMean/double(N);
//  //checking mean
//  if ( expmean!=0 ){
//          BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//  }
//  BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*abs(abs(min)-abs(expmean)));
//  //checking minimum
//  BOOST_CHECK(hist.begin()->first >= min);
//  //checking maximum
//  BOOST_CHECK(hist.rbegin()->first <= max );
//  if( print_histogram ){
//      cout << "Histogram : " << endl;
//      for(auto p : hist) {
//              std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//                        << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//      }
//      cout << "\n";
//  }
//  cout << "Expected Mean: " << expmean << endl;
//  cout << "Generated Mean: " << foundMean << endl;
//}
//
//
//BOOST_AUTO_TEST_CASE(Cosine_Distribution)
//{
//  printBannerTestCase( "Testing Cosine Random Distribution" );
//  int N = 1000000;
//  std::map<int, int> hist;
//  double expmean, min, max, foundMean;
//  min = -pi/2;
//  max = pi/2;
//  expmean = (min + max)/2;
//  vec angles = linspace(min, max, (max - min)/100.0);
//  vec tempweights = cos( angles );
//  vec weights = tempweights/ accu( tempweights );
//  vector<double> anglesStdVec = arma::conv_to< vector<double> >::from( angles );
//  vector<double> weightsStdVec = arma::conv_to< vector<double> >::from( weights );
//  foundMean = 0;
//  bool reset = false;
//  Distribution cosineRandDist( DistributionType::DISCRETE, min, max, anglesStdVec, weightsStdVec);
//  std::cout << cosineRandDist.toString() << endl;
//  cosineRandDist.getDistribution(true);
//  for( int n=0; n<N; n++ ){
//      double number = cosineRandDist.getDistribution(reset);
//      foundMean += number;
//      ++hist[std::round( number)];
//  }
//  foundMean = foundMean/double(N);
//  //checking mean
//  if ( expmean!=0 ){
//          BOOST_CHECK_CLOSE( foundMean,  expmean, 1);
//  }
//  BOOST_CHECK(abs( expmean - abs(foundMean) ) < 0.1*abs(abs(min)-abs(expmean)));
//  //checking minimum
//  BOOST_CHECK(hist.begin()->first >= min);
//  //checking maximum
//  BOOST_CHECK(hist.rbegin()->first <= max );
//  if( print_histogram ){
//      cout << "Histogram : " << endl;
//      for(auto p : hist) {
//              std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//                        << p.first << ' ' << std::string(p.second/4000, '*') << '\n';
//      }
//      cout << "\n";
//  }
//  cout << "Expected Mean: " << expmean << endl;
//  cout << "Generated Mean: " << foundMean << endl;
//}





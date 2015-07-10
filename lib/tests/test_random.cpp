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

using namespace utils::random;
using namespace std;

BOOST_AUTO_TEST_CASE(basics)
{
    int N = 1000000;
    vector<double> vs(N);
    genNormalDist(vs, 1, 0);
    
    double mean = 0;
    for (auto v : vs) {
        mean += v;
    }

    mean = mean/N;
    BOOST_CHECK_CLOSE(0.0, abs(mean), 1);

    cout << "Mean: " << mean << endl;
}



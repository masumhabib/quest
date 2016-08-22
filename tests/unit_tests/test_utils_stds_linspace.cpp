/** 
 * Test cases for functions in utils::stds.
 *
 */

#include "utils/random.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LinspaceTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>
#include <map>
#include <string>
#include <iomanip>

using namespace utils::random;
using namespace std;



BOOST_AUTO_TEST_CASE(linspace_with_count_generates_int_vector)
{
    size_t N = 0;
    double min = 0, max = 9;

    auto vec = utils::stds::linspace(min, max, N);
    BOOST_REQUIRE_EQUAL(vec.size(), N);


    N = 1;
    vec = utils::stds::linspace(min, max, N);
    BOOST_REQUIRE_EQUAL(vec.size(), N);
    BOOST_CHECK_EQUAL(vec[0], min);

    N = 2;
    vec = utils::stds::linspace(min, max, N);
    BOOST_REQUIRE_EQUAL(vec.size(), N);
    BOOST_CHECK_EQUAL(vec[0], min);
    BOOST_CHECK_EQUAL(vec[1], max);

    N = 10;
    vec = utils::stds::linspace(min, max, N);
    BOOST_REQUIRE_EQUAL(vec.size(), N);
    BOOST_CHECK_EQUAL(vec[0], min);
    BOOST_CHECK_EQUAL(vec[4], 4);
    BOOST_CHECK_EQUAL(vec[5], 5);
    BOOST_CHECK_EQUAL(vec[N-1], max);

}


BOOST_AUTO_TEST_CASE(linspace_with_delta_generates_int_vector)
{
    double min, max, delta;

    min = 0, max = 0, delta = 0;
    auto vec = utils::stds::linspace(min, max, delta);
    BOOST_REQUIRE_EQUAL(vec.size(), 0);


    min = 1, max = 1, delta = 1;
    vec = utils::stds::linspace(min, max, delta);
    BOOST_REQUIRE_EQUAL(vec.size(), 1);
    BOOST_CHECK_EQUAL(vec[0], min);

    min = 0, max = 1, delta = 1;
    vec = utils::stds::linspace(min, max, delta);
    BOOST_REQUIRE_EQUAL(vec.size(), 2);
    BOOST_CHECK_EQUAL(vec[0], min);
    BOOST_CHECK_EQUAL(vec[1], max);

    min = 0, max = 9, delta = 1;
    vec = utils::stds::linspace(min, max, delta);
    BOOST_REQUIRE_EQUAL(vec.size(), 10);
    BOOST_CHECK_EQUAL(vec[0], min);
    BOOST_CHECK_EQUAL(vec[4], 4);
    BOOST_CHECK_EQUAL(vec[5], 5);
    BOOST_CHECK_EQUAL(vec[9], max);

}




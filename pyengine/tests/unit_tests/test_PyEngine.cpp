/** Test cases for Segment class.
 *
 */

#include "PyEngine.hpp"

#ifndef LINK_STATIC
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE PyEngineTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>
#include <map>
#include <string>
#include <iomanip>

using namespace std;
using namespace quester;

BOOST_AUTO_TEST_SUITE (PyEngine_constructor);

BOOST_AUTO_TEST_CASE(PyEngine_default_construction_initializes_python_engine)
{
    PyEngine py;
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE (PyEngine_eval);

BOOST_AUTO_TEST_CASE(PyEngine_run_hello_world_succeeds)
{
    PyEngine py;
    auto status = py.eval ("print('Hello World!')");
    //BOOST_REQUIRE (status == PyEngine::CommandStatus::SUCCESS);

    //auto output = py.get_output_msg ();
    //BOOST_CHECK_EQUAL (output, "Hello World!");
}

BOOST_AUTO_TEST_CASE(PyEngine_run_simple_variable_succeeds)
{
    PyEngine py;
    //auto status = py.eval ("x = 10; print (x)");
    //BOOST_REQUIRE (status == PyEngine::CommandStatus::SUCCESS);

    //auto output = py.get_output_msg ();
    //BOOST_CHECK_EQUAL (output, "10");
}

BOOST_AUTO_TEST_SUITE_END();


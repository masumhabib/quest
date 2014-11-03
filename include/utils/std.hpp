/* 
 * File:   std.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 9, 2014, 12:46 AM
 */

#ifndef STD_HPP
#define	STD_HPP

#include <ios>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <cmath>
#include <memory>

namespace utils{
namespace stds{

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

using std::ios;
using std::ios_base;
using std::basic_ostream;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::ostringstream;
using std::istringstream;

using std::string;
using std::vector;
using std::list;
using std::pair;
using std::map;

using std::abs;

using std::swap;

using std::invalid_argument;
using std::runtime_error;

using std::binary_function;

using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;
using std::static_pointer_cast;

}
}


#endif	/* STD_HPP */


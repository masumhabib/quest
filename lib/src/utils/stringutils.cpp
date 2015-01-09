/* 
 * File:   stringutils.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 7, 2013, 10:09 PM
 * 
 * Description: String utilities.
 * 
 */

#include "utils/stringutils.h"

namespace utils{namespace strings{

/**
 * trim:
 * Remove surrounding whitespace from a std::string.
 * @param s The string to be modified.
 * @param t The set of characters to delete from each end
 * of the string.
 * @return The same string passed in as a parameter reference.
 */

string trim(const string& ss, const char* t){
    string s = ss;
	s.erase(0, s.find_first_not_of(t));
	s.erase(s.find_last_not_of(t) + 1);
	return s;
}

}}


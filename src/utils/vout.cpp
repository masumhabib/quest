/* 
 * File:   vout.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 16, 2014, 12:31 PM
 */

#include "vout.h"

namespace utils{
namespace stds{

/**
 * Global variables.
 */

const verbosity  vquiet(0);
const verbosity  vnormal(1);
const verbosity  vdebug(10);
const verbosity  vdump(100);

vostream         dout(cerr);
vostream         vout(cout);

const string     dbg = "  DBG: ";

}
}
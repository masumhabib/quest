/* 
 * File:   RgfResult.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 13, 2014, 3:41 PM
 */

#ifndef RGFRESULT_H
#define	RGFRESULT_H

#include "maths/arma.hpp"
#include "utils/std.hpp"
#include "utils/serialize.hpp"

namespace quest{
namespace negf{

using namespace utils::stds;
using namespace maths::armadillo;

/**
 * Result as a function of energy.
 */

struct RgfResult{
public:    
    typedef std::list<cxmat>::iterator iter;
    
    std::list<cxmat>R;           //!< List of results.
    int             ib;          //!< Block i.
    int             jb;          //!< Block j.
    uint            N;           //!< Size of the cxmat matrix. N == 0 means do not calculate.
    string          tag;         //!< Suffix added to the output filename.
    
    RgfResult(string tag = "", uint N = 0, int ib = -1, int jb = -1);
    
    bool    isEnabled() { return N > 0; };  //!< Is calculation enabled? N = 0: no.
    void    save(ostream &out, bool isText);
};

}
}

#endif	/* RGFRESULT_H */


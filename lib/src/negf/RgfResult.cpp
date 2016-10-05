/* 
 * File:   RgfResult.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 13, 2014, 3:41 PM
 */

#include "negf/RgfResult.h"

namespace quest{
namespace negf{

RgfResult::RgfResult(string tag, uint N, int ib, int jb):
        tag(tag), N(N), ib(ib), jb(jb)
{
}

void RgfResult::save(ostream &out, bool isText){
    if (isText){
        out << tag << endl;
        out << R.size() << endl;
        out << ib << " " << jb << endl;
        out << N << endl;
        iter it;
        for (it = R.begin(); it != R.end(); ++it){
            out << *it << endl;
        }
    }else{
        
    }
}
 
}
}

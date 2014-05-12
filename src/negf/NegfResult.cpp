/* 
 * File:   NegfResult.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 13, 2014, 3:41 PM
 */

#include "NegfResult.h"

namespace qmicad{
namespace negf{

NegfResultList::NegfResultList(string tag, uint N, int ib, int jb):
        tag(tag), N(N), ib(ib), jb(jb)
{
}

void NegfResultList::sort(){
    R.sort(ResultComparator());
}

void NegfResultList::merge(NegfResultList &second){
    R.merge(second.R, ResultComparator());
}

void NegfResultList::save(ostream &out){
    out << tag << endl;
    out << ib << " " << jb << endl;
    out << R.size() << endl;
    out << N << endl;
    iter it;
    for (it = R.begin(); it != R.end(); ++it){
        out << it->E << endl << it->M << endl;
    }
}
 
}
}
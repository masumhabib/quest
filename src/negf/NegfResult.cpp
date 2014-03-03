/* 
 * File:   NegfResult.cpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 * 
 * Created on February 13, 2014, 3:41 PM
 */

#include "NegfResult.h"

namespace qmicad{

NegfResultList::NegfResultList(string tag, uint N):
        tag(tag), N(N)
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
    out << R.size() << endl;
    out << N << endl;
    list<negf_result>::iterator it;
    for (it = R.begin(); it != R.end(); ++it){
        out << it->first << endl << it->second << endl;
    }
}
 
}
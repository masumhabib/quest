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
    list<negf_result>::iterator it;
    for (it = R.begin(); it != R.end(); ++it){
        if (N == 1){
            row rvec(2); 
            rvec(0) = it->first;                   // energy
            rvec(1) = real((it->second)(0,0));     // result
            out << rvec;  
        }else{
            out << it->first << endl << it->second << endl;
        }
    }
}
 
}
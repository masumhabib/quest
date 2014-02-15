/* 
 * File:   NegfResult.cpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 * 
 * Created on February 13, 2014, 3:41 PM
 */

#include "NegfResult.h"

namespace qmicad{

NegfResultList::NegfResultList(string suffix, uint N, bool saveAscii):
        suffix(suffix), N(N), saveAscii(saveAscii)
{
}

void NegfResultList::sort(){
    R.sort(ResultComparator());
}

void NegfResultList::merge(NegfResultList &second){
    R.merge(second.R, ResultComparator());
}

void NegfResultList::save(string fileName){
    fileName = fileName + suffix + ".dat";
    // save to a file
    ofstream outFile;
    if(saveAscii){
        outFile.open(fileName.c_str(), ostream::binary);
    }else{
        outFile.open(fileName.c_str());
    }
    if (!outFile.is_open()){
        throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                + fileName + ".");
    }

    list<negfresult>::iterator it;
    for (it = R.begin(); it != R.end(); ++it){

        if (N == 1){
            row rvec(2); 
            rvec(0) = it->first;                   // energy
            rvec(1) = real((it->second)(0,0));     // result
            outFile << rvec;  
        }else{
            outFile << it->first << it->second << endl;
        }
    }
    outFile.close();
}
 
}
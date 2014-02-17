/* 
 * File:   parloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 2:41 PM
 */

#include "parloop.h"

namespace utils{
    

ParLoop::ParLoop(const Workers &workers, int N):
        mWorkers(workers)
{
    assign(N);
}    

void ParLoop::assign(int N){
    mN = N;

    int ncpu = mWorkers.N();
    int myId = mWorkers.MyId();
    
    int quo = mN/ncpu ;
    int rem = mN%ncpu;
    mMyStart = myId*quo + (myId < rem ? myId:rem);
    mMyEnd = mMyStart + quo + (myId < rem ? 1:0) - 1;        
}
            
void ParLoop::run(){
    prepare();
    for(int il = mMyStart; il <= mMyEnd; ++il){
        preCompute(il);
        compute(il);
        postCompute(il);
    }    
    collect();
}
       
}
    
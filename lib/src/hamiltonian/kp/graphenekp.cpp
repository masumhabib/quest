/* 
 * File:   graphenekp.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */


#include <maths/svec.h>

#include "hamiltonian/kp/graphenekp.h"

namespace qmicad{
namespace hamiltonian{

//GrapheneKpParams::GrapheneKpParams(const string &prefix) //:DiracKpParams(prefix)
//{
//    mTitle  = "Graphene k.p parameters";
//}

string GrapheneKpParams::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ":" << endl;
    ss << mPrefix << " a     = " << ma << endl;
    ss << mPrefix << " K     = " << mK  << endl;
    ss << mPrefix << " gamma = " << mgamma << endl;
    ss << mPrefix << " dtol  = " << mdtol << endl;
    ss << mPrefix << " Bz    = " << mBz << endl;
    ss << mPrefix << " A     = " << (mBzGauge == coord::X?" (-Bz*y, 0, 0)":" (0, Bz*x, 0)") << endl;
    ss << mPrefix << " eps: " << endl << meps;
    ss << mPrefix << " t01x: " << endl << mt01x;
    ss << mPrefix << " t01y: " << endl << mt01y;

    return ss.str(); 
};

GrapheneOneValleyKpParams::GrapheneOneValleyKpParams(const string &prefix)
    :GrapheneKpParams(prefix)
{
    mI = eye<cxmat>(2,2);    
    setDefaultParams();
    update();
}

void GrapheneOneValleyKpParams::update(){
    if(!(is_finite(mdtol) && is_finite(mgamma) && is_finite(ma) 
            && is_finite(mK))){
        throw runtime_error("GrapheneOneValleyKpParams: invalid k.p parameters.");    
    }
    double Kx = mK, Ky = mK, ax = ma, ay = ma;
    meps = -mgamma*(Kx/ax + Ky/ay)*sz();
    mt01x = (mgamma/(2*ax))*sx()*i + (Kx*mgamma/(2*ax))*sz();
    mt10x = trans(mt01x);
    mt01y = (mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
    mt10y = trans(mt01y);                
}

GrapheneTwoValleyKpParams::GrapheneTwoValleyKpParams(const string &prefix)
    :GrapheneKpParams(prefix)
{
    mI = eye<cxmat>(4, 4);    
    setDefaultParams();
    update();
}

void GrapheneTwoValleyKpParams::update(){
    if(!(is_finite(mdtol) && is_finite(mgamma) && is_finite(ma) 
            && is_finite(mK))){
        throw runtime_error("GrapheneTwoValleyKpParams: invalid k.p parameters.");    
    }
    
    double Kx = mK, Ky = mK, ax = ma, ay = ma;
    meps = zeros<cxmat>(4,4);
    cxmat eps22 = -mgamma*(Kx/ax + Ky/ay)*sz();
    meps(span(0,1), span(0,1)) = eps22;
    meps(span(2,3), span(2,3)) = eps22;

    mt01x = zeros<cxmat>(4,4);
    cxmat t01x_K = (mgamma/(2*ax))*sx()*i + (Kx*mgamma/(2*ax))*sz();
    cxmat t01x_Kp = (-mgamma/(2*ax))*sx()*i + (Kx*mgamma/(2*ax))*sz();
    mt01x(span(0,1), span(0,1)) = t01x_K;
    mt01x(span(2,3), span(2,3)) = t01x_Kp;
    mt10x = trans(mt01x);

    mt01y = zeros<cxmat>(4,4);
    cxmat t01y_K = (mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
    // cxmat t01y_Kp = (-mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
    cxmat t01y_Kp = (mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
    mt01y(span(0,1), span(0,1)) = t01y_K;
    mt01y(span(2,3), span(2,3)) = t01y_Kp;
    mt10y = trans(mt01y);
}

}
}



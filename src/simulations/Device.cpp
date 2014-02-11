/* 
 * File:   Device.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 * 
 * Created on January 25, 2014, 7:37 PM
 * 
 */

#include "Device.h"
namespace qmicad{

Device::Device(const communicator &workers, const DeviceParams &p):
interval(0), mWorkers(workers), p(p)
{
    
}

void Device::tic(){
    timer.tic();
}

void Device::toc(){
    interval = timer.toc();
}

string Device::time(){
    return " Runtime: " + ttos(interval) + ".";
}

void Device::setDeviceParams(const DeviceParams &dp){
    p = dp;
}

void Device::addSource(const squadrilateral &sq){
    // source
    source = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    source.Title("Source");
}

void Device::addDrain(const squadrilateral &sq) {
    // drain
    drain = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    drain.Title("Drain");
}

void Device::addGate(const squadrilateral& sq){
    gates.push_back(gate(sq.lb, sq.rb, sq.rt, sq.lt)); 
    int it = gates.size() - 1;
    gates[it].Title("Gate # " + itos(it+1));
}

void Device::addLinearRegion(const squadrilateral& sq){
    linrs.push_back(linear_region(sq.lb, sq.rb, sq.rt, sq.lt));    
    int it = linrs.size() - 1;
    linrs[it].Title("Linear Region # " + itos(it+1));
}

string Device::toString() const{
    stringstream ss;
    ss << "Device Information: " << endl;
    ss  << source << endl;
    ss  << drain << endl;
    for(int it = 0; it < gates.size(); ++it){
        ss << gates[it] << endl;
    }
    for(int it = 0; it < linrs.size(); ++it){
        ss << linrs[it] << endl;
    }
    
    return ss.str();
}

void Device::prepare(){
/*
    HamMethod hMethod;      // k.p, tight binding, EHT
    BasisType bType;        // orthogonal/non-orthogonal
    DeviceType dType;       // uniform/non-uniform

    // Prepare parameters required for Hamiltonina generation
    shared_ptr<HamParams> hp;       // Hamiltonian parameter
    if( p.HamiltonianType == "GrapheneKp"){
            shared_ptr<GrapheneKpParams> tmp(new GrapheneKpParams());
            tmp->dtol = p.dtol;
            tmp->ax = p.ax;
            tmp->ay = p.ay;
            tmp->K = p.K;
            tmp->gamma = p.gamma;
            tmp->update();
            hp = tmp;
            
            hMethod = kp;
            bType = Orthogonal;
            dType = Uniform;
        }else if(p.HamiltonianType == "TISurfKp"){
            shared_ptr<TISurfKpParams> tmp(new TISurfKpParams());
            tmp->dtol = p.dtol;
            tmp->ax = p.ax;
            tmp->ay = p.ay;
            tmp->K = p.K;
            tmp->A2 = p.A2;
            tmp->C = p.C;
            tmp->update();
            hp = tmp;

            hMethod = kp;
            bType = Orthogonal;
            dType = Uniform;
        }else{
            string msg = "Invalid HamiltoninaType. Allowed types are: \n";
            msg += "GrapheneKp\n";
            msg += "TISurfKp\n";
            throw runtime_error(msg);
    }    

    // Generate atomistic geometry
    // For k.p method
    if (hMethod == kp){
        double ax = p.ax;
        double ay = p.ay;
        int nw = p.nw;
        int nl = p.nl;
        w = (nw-1)*ax;                  // width (in A)
        l = (nl-1)*ay;                  // length  (in A)
        // Create fake atoms for discretized k.p model
        da = genKpAtoms(l, w, ax, ay, hp->PeriodicTable); 
        
    }
    na = da.NumOfAtoms();           // number of atoms
    no = da.NumOfOrbitals();        // number of orbitals   
    xMin = min(da.X());
    xMax = max(da.X());
    yMin = min(da.Y());
    yMax = max(da.Y());

    // Generate H and S
    // Hamiltonian and overlap matrix: 
    // for uniform device, we just need to store one H_i,i and 
    // one H_i+1,i
    if (dType == Uniform){
        AtomicStruct lyr0 = da(span(0, p.nw-1));        // take block 0
        AtomicStruct lyr1 = da(span(p.nw, 2*p.nw-1));   // take block 1
        H0.set_size(1);
        Hl.set_size(1);
        if(typeid(*hp) == typeid(GrapheneKpParams)){
            shared_ptr<GrapheneKpParams> grkpp = boost::static_pointer_cast<GrapheneKpParams>(hp);
            GrapheneKpHamGen generator(*grkpp);
            H0(0) = genHam<cxmat>(lyr0, lyr0, generator);   // H_0,0
            Hl(0) = genHam<cxmat>(lyr1, lyr0, generator);   // H_0,-1 = H_1,0
        }else if(typeid(*hp) == typeid(TISurfKpParams)){
            shared_ptr<TISurfKpParams> tikpp = boost::static_pointer_cast<TISurfKpParams>(hp);
            TISurfHamGen generator(*tikpp);
            H0(0) = genHam<cxmat>(lyr0, lyr0, generator);   // H_0,0
            Hl(0) = genHam<cxmat>(lyr1, lyr0, generator);   // H_0,-1 = H_1,0
            
        }
        // Overlap
        S0.set_size(1);
        S0(0).copy_size(H0(0));                          // S_0,0 = I        
        if (bType == Orthogonal){
            S0(0).eye();
        }else{
            Sl.set_size(1);
            // calculatate S0 and Sl
        }
    }else{
        // non-uniform
    }
      
    // Generate block Hamiltonian for NEGF parameters
    int nl = p.nl;    
    // create H and S arrays for NEGF
    np.H0.set_size(nl);
    np.S0.set_size(nl);
    np.Hl.set_size(nl+1);
    if (bType == Nonorthogonal){
        np.Sl.set_size(nl+1);
    }
    for (int ib = 0; ib <= nl; ++ib){
        int ibl;
        if (dType == Uniform){
            ibl = 0;
        }else{
            ibl = ib;
        }
        if (ib != nl){
            np.H0(ib) = &H0(ibl);
            np.S0(ib) = &S0(ibl);
        }
        // allow N+1 in H_i,i-1 and S_i,i-1
        np.Hl(ib) = &Hl(ibl);
        if (bType == Nonorthogonal){
            np.Sl(ib) = &Sl(ibl);
        }
    }
      
    // Other NEGF parameters
    np.nb = nl;
    np.isOrthogonal = bType == Orthogonal;
    np.V.set_size(p.nl);
    np.kT = p.kT;
    np.ieta = dcmplx(0,p.eta);        // imaginary potential
   
    // NEGF setup
    np.grcCache = myenums::Enabled; // enable grc cache
    np.DCache = myenums::Enabled;   // enable D cache

    // biases
    VDDs = BiasGrid(p.VDDmn, p.VDDmx, p.dVDD);
    VGGs = BiasGrid(p.VGmin, p.VGmax, p.dVG);  
 */  
}

void Device::VDS(double VD, double VS){
    np.muD = p.mu - VD;
    np.muS = p.mu - VS;
}

void Device::VG(int ig, double VG){
    gates[ig].V = VG;
}

void Device::VLR(int ilr, double Vl, double Vr, double Vt, double Vb){
    linrs[ilr].Vl = Vl;
    linrs[ilr].Vr = Vr;
    linrs[ilr].Vt = Vt;
    linrs[ilr].Vb = Vb;
}


void Device::computePotential(){
/*
    // block potential vector
    V.set_size(p.nl);
    if (p.PotentialSolver == "Linear"){
        //set drain and source potential
        source.V = gates[0].V;
        drain.V = gates[gates.size()-1].V;
        
        // atomistic linear potential
        LinearPot Va(da, source, drain, gates, linrs);
        Va.compute();
        // convert it to block potential
        for (int ib = 0; ib < p.nl; ++ib){
            V(ib) = Va.toOrbPot(span(p.nw*ib, p.nw*(ib+1)-1));
            np.V(ib) = &V(ib);
        }
        //Va.exportPotential("pot.dat");
    }else{
        throw runtime_error("ERROR: Forgot to set potential solver type?");
    }
 */
}

void Device::runNegfEloop(){
    /*
    VecGrid EE;
    // prepare energy grid
    if(p.AutoGenE){
        EE = VecGrid(np.muD - 10*p.kT, np.muS + 10*p.kT, p.dE);
    }else{
        EE = VecGrid(p.Emin, p.Emax, p.dE);
    }    
    NegfEloop Eloop(EE, np, mWorkers);
    Eloop.run();
    Eloop.saveTE(p.OutFileName + "TE.dat");
   */
}

NegfParams Device::NegfParam() { 
    return np; 
};

}


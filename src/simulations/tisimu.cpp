/* 
 * File:   tisimu.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 26, 2014, 2:57 PM
 * 
 * Topological Insulator device simulation file
 * 
 */

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "../include/qmad.hpp"
#include "tisimu.h"

using namespace std;

#ifndef NO_MPI
namespace mpi = boost::mpi;
#endif

using constants::pi;

#ifndef NO_MPI
void tiTrans(mpi::communicator &world){
#else
void tiTrans(){
#endif
    
    double dVD = 0.05;
    double VDmax = 0.0;
    double VDmin = 0.0;
    
    double dVG = 0.05;
    double VGmax = 1.0;
    double VGmin = 1.0;
      
    // setup
    double kT = 0.0259;      // Temperature in eV 
    double eta = 1E-3;       // imaginary potential    
    
    double mu = 0.0;         // Fermi level
    double Ew = 3;           // energy window
    double dE = 0.01;       // energy grid spcing
    
    int nl = 41;             // number of blocks in the device including contact
    int nw = 33;            // number of layers along the width
    double ax = 4;           // discretization length in x direction
    double ay = 4;           // discretization length in y direction
    double w = (nw-1)*ax;    // width (in A)
    double l = (nl-1)*ay;    // length  (in A)
    double Lx = l + ax;      // lattice vector ( along x)
    int dg = 20;             // distance between gate1 and gate2
    double thg = 10*pi/180;  // tilt angle

    int na = nl * nw;        // number of atoms
    int no = 2*nl*nw;        // number of orbitals
    
    // k.p params for Be2Se3
    double gamma = 3.16*1.42*3/2;   // h_bar*v_F
    double dtol = 1E-2;             // distance tolerance. 
    double K = ax;                  // avoid fermion doubling  
    
    string TEFileName = "gnr";
    
    // TI parameters
    GrapheneKpParams grb;
    grb.dtol = dtol;
    grb.gamma = gamma;
    grb.K = K;
    grb.ax = ax;
    grb.ay = ay;
    grb.update();
        
    // Create fake atoms for discretized k.p model
    Atoms atoms = genKpAtoms(l, w, ax, ay, grb.PeriodicTable);    
    
    // Define the terminals
    double minx = min(atoms.X()) - ax/2;
    double maxx = max(atoms.X()) + ax/2;
    double miny = min(atoms.Y()) - ay/2;
    double maxy = max(atoms.Y()) + ay/2;
    double ds = (dg-1)*ax;
    // source
    contact source(point(minx, miny), point(minx+ax, miny), 
            point(minx+ax, maxy),  point(minx, maxy));
    // drain
    contact drain(point(maxx-ax, miny), point(maxx, miny), 
            point(maxx, maxy), point(maxx-ax, maxy));
    // gates
    vector<gate> gates;
    // gate # 1
    gates.push_back(gate(point(minx+ax, miny), point(-ds/2, miny), 
            point(-ds/2, maxy), point(minx+ax, maxy)));
    // gate # 2
    gates.push_back(gate(point(ds/2, miny), point(maxx-ax, miny),
            point(maxx-ax, maxy), point(ds/2, maxy)));
    // linear region
    vector<linear_region> split;
    split.push_back(linear_region(point(-ds/2, miny), point(ds/2, miny), 
            point(ds/2, maxy), point(-ds/2, maxy)));
    
    // Hamiltonian: for uniform device, we just need to store one H_i,i and 
    // one H_i+1,i
    cx_mat H00, S00, H10, S10;                  // H and S
    NegfParams np;                              // NEGF parameters
    // Generate TI block Hamiltonian
    {
        Atoms lyr0 = atoms(span(0, nw-1));      // take block 0
        Atoms lyr1 = atoms(span(nw, 2*nw-1));   // take block 1
        
        // Generate H and S for only one block
        GrapheneKpHamGen generator(grb);
        H00 = genHam<cx_mat>(lyr0, lyr0, generator);                // H_0,0
        H10 = genHam<cx_mat>(lyr1, lyr0, generator);                // H_0,-1 = H_1,0
        S00.copy_size(H00);                      // S_0,0 = I
        S00.eye();
        
        // create H and S arrays for NEGF
        np.nb = nl;
        np.H0.set_size(nl);
        np.S0.set_size(nl);
        np.Hl.set_size(nl+1);
        for (int ib = 0; ib < nl; ++ib){
            np.H0(ib) = &H00;
            np.S0(ib) = &S00;
            np.Hl(ib) = &H10;
        }
        np.Hl(nl) = &H10;    
    }
    np.V.set_size(nl);
    np.kT = kT;
    np.calcTE = true;
    
    // NEGF setup
    np.grcCache = myenums::Enabled; // enable grc cache
    np.DCache = myenums::Enabled;   // enable D cache
    np.isOrthogonal = true;         // orthogonal basis set
    np.ieta = dcmplx(0,eta);        // imaginary potential
    
    // biases
    BiasGrid VD(VDmin, VDmax, dVD);
    BiasGrid VG(VGmin, VGmax, dVG);
    
    // bias loops
    for(int ivd = 0; ivd < VD.N; ++ivd){
        alculate();
        np.muD = mu - 0.5*VD(ivd);     // drain voltage
        for(int ivg = 0; ivg < VG.N; ++ivg){
            // setup terminal voltages
            double VG1 = +0.5*VG(ivg);         // Gate1 voltage
            double VG2 = -0.5*VG(ivg);         // Gate2 voltage
            source.V = VG1;                    //
            drain.V = VG2;
            gates[0].V = VG1;
            gates[1].V = VG2;
            split[0].Vl = VG1;
            split[0].Vr = VG2;
            // construct linear potential profile
            field<vec> V(nl);
            {
                // atomistic potential
                LinearPot Va(atoms, source, drain, gates, split);
                Va.calculate();
                // convert it to block potential
                for (int ib = 0; ib < nl; ++ib){
                    V(ib) = Va.toOrbPot(span(nw*ib, nw*(ib+1)-1));
                    np.V(ib) = &V(ib);
                }
                Va.exportPotential("pot.dat");
            }
            
            // Do NEGF
#ifndef NO_MPI          
            VecGrid EE(-Ew/2, Ew/2, dE);
            NegfLoop negfLoop(EE, np, world);
#else
            NegfLoop negfLoop(VecGrid(-Ew/2, Ew/2, dE), np);
#endif
            negfLoop.run();
            negfLoop.saveTE(TEFileName + ".dat");
            
        }
    }
     
}

/* 
 * File:   Device.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 25, 2014, 7:37 PM
 */

#ifndef DEVICE_H
#define	DEVICE_H

#include <typeinfo> 
#include <armadillo>
#include <boost/smart_ptr.hpp>

#include "../include/qmicad.hpp"

namespace qmicad{
using namespace arma;
using boost::shared_ptr;
using boost::static_pointer_cast;
using boost::mpi::communicator;

struct DeviceParams{
    // Geomtetry
    int             nl;          // number of blocks in the device including contact
    int             nw;          // number of layers along the width
    // Hamiltonina parameters
    double          dtol;        // distance tolerance. 
    // k.p parameters
    double          ax;          // discretization length in x direction
    double          ay;          // discretization length in y direction
    double          K;           // avoid fermion doubling  
    // for graphene
    double          gamma;       // h_bar*v_F
    // for TI
    double          A2;          // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double          C;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)  
    // Bias
    double          dVDD;         // drain bias interval
    double          VDDmx;       // maximum drain bias 
    double          VDDmn;       // minimum drain bias
    double          dVG;         // Gate bias interval
    double          VGmax;       // maximum drain bias
    double          VGmin;       // minimum drain bias
    string          PotSolveType;// potential solver type
    // NEGF parameters
    double          kT;          // Temperature in eV 
    double          eta;         // imaginary potential    
    double          mu;          // Fermi level
    bool            AutoGenE;    // Generate auto energy grid
    double          Emin;        // Min energy
    double          Emax;        // Max energy
    double          dE;          // Energy grid spacing
    
    
    string  OutFileName;         // Filename
    string  gjfFileName;         // GJF file name
    string  HamiltonianType;     // graphene/TI?
    string  PotentialSolver;     // Linear/Poisson/Laplace?
   
    bool            CalcTE;      // Calculate TE?

};

class Device { 
public:
    enum            BasisType{ Orthogonal, Nonorthogonal };
    enum            HamMethod{ kp, TightBinding, EHT };
    enum            DeviceType{ Uniform, NonUniform };
    
protected:    
    const communicator &mWorkers;  // MPI communicator
    DeviceParams    p;             // Device parameters
 
    // Atomistic geometry
    Atoms           da;            // Atomistic geometry of the device containing two contact leads
    double          w;             // width of the device (in A)
    double          l;             // length of the device (in A)
    double          xMin;          // x min of the atoms
    double          xMax;          // x max of the atoms
    double          yMin;          // y min of the atoms
    double          yMax;          // y max of the atoms
    double          na;            // total number of atoms in the device
    double          no;            // total number of orbitals in the device    
        
    // NEGF
    field<cx_mat>   H0;            // Diagonal blocks of Hamiltonian: H0(i) = [H]_i,i
    field<cx_mat>   Hl;            // Lower diagonal blocks of Hamiltonian: Hl(i) = [H]_i,i-1
    field<cx_mat>   S0;            // Diagonal blocks of overlap matrix: S0(i) = [S]_i,i
    field<cx_mat>   Sl;            // Lower diagonal blocks of overlap matrix: 
    NegfParams      np;            // NEGF parameters
        
    // Bias
    BiasGrid        VDDs;          // Drain bias
    BiasGrid        VGGs;          // Gate bias
    contact         source;        // source contact
    contact         drain;         // source contact
    vector<gate>    gates;         // gates array
    vector<linear_region> linrs;   // gates array
    field<vec>      V;             // block potential matrixes
    
    // utility
    wall_clock      timer;         // for tic toc: time measurement
    double          interval;      // total time elapsed

public:
    Device(const communicator &worker, const DeviceParams &p);
    
    void        tic();
    void        toc();
    string      time();
    
    void        setDeviceParams(const DeviceParams &dp);
    void        prepare();
    void        computePotential();
    void        runNegfEloop();
    
    void        addSource(const squadrilateral &sq);
    void        addDrain(const squadrilateral &sq);
    void        addGate(const squadrilateral &sq);
    void        addLinearRegion(const squadrilateral &sq);
    
    double      xmin() {return xMin; };
    double      xmax() {return xMax; };
    double      ymin() {return yMin; };
    double      ymax() {return yMax; };
    
    double      VDD(int ivd) { return VDDs(ivd); };
    int         NVDD() { return VDDs.N(); };
    double      VGG(int ivg) { return VGGs(ivg); };
    int         NVGG() { return VGGs.N(); };
    
    void        VDS(double VD, double VS = 0);
    void        VG(int ig, double VG);
    void        VLR(int ilr, double Vl, double Vr, double Vt = 0, double Vb = 0);
    
    NegfParams NegfParam();
    
    virtual string toString() const;

protected:    
private:

};

}
#endif	/* DEVICE_H */


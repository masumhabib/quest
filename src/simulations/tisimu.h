/* 
 * File:   tisimu.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 26, 2014, 2:57 PM
 */

#ifndef TISIMU_H
#define	TISIMU_H


#ifndef NO_MPI
#include <boost/mpi.hpp>
void tiTrans(boost::mpi::communicator &world);
#else
void tiTrans();
#endif

#endif	/* TISIMU_H */


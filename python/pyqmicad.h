/* 
 * File:   pyqmicad.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python wrapper for QMICAD package
 * 
 */

#ifndef PYQMICAD_H
#define	PYQMICAD_H

namespace qmicad{
namespace python{

void export_Option();

void export_cxmat();
void export_mat();
void export_vec();

void export_VecGrid();

void export_Printable();

void export_svec();
void export_pvec();
void export_lvec();
void export_lcoord();
void export_Atom();
void export_PeriodicTable();
void export_AtomicStruct();

void export_Workers();

void export_point();
void export_quadrilateral();

void export_Timer();

void export_cxhamparams();
void export_GrapheneTbParams();
void export_GrapheneKpParams();
void export_TISurfKpParams();
void export_TISurfKpParams4();

void export_Potential();
void export_LinearPot();

void export_CohRgfLoop();

void export_KPoints();

void export_BandStructParams();
void export_BandStruct();

void export_npyarma();


}
}
#endif	/* QMICAD_H */


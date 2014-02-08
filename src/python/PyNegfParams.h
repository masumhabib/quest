/* 
 * File:   PyNegfParams.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 8, 2014, 2:28 AM
 */

#ifndef PYNEGFPARAMS_H
#define	PYNEGFPARAMS_H

#include "../negf/NEGF.h"

namespace qmicad{
namespace python{

class PyNegfParams:public NegfParams{
public:
    PyNegfParams():NegfParams(){};
};

}
}
#endif	/* PYNEGFPARAMS_H */


/* 
 * File:   NegfResult.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on February 13, 2014, 3:41 PM
 */

#ifndef NEGFRESULT_H
#define	NEGFRESULT_H

#include "../maths/arma.hpp"
#include "../utils/std.hpp"

namespace qmicad{
using namespace utils::stds;
using namespace maths::armadillo;

typedef pair<double, cxmat> negfresult;     //!< Result as a function of energy.

struct ResultComparator : public binary_function<negfresult&, negfresult&, bool>{
    ResultComparator( double tol = 1e-7 ) : epsilon(tol) {}
    bool operator()( const negfresult &left, const negfresult &right  ) const{
        return (abs(left.first - right.first) > epsilon) && (left.first < right.first);
    }   
    double epsilon;
};

struct NegfResultList{
public:    
    std::list<negfresult>R;              //!< List of results.
    uint            N;              //!< Size of the cxmat matrix. N == 0 means do not calculate.
    bool            saveAscii;      //!< Save as ascii/binary?
    string          suffix;         //!< Suffix added to the output filename.
    
    NegfResultList(string suffix, uint N = 0, bool saveAscii = true);
    
    bool    isEnabled() { return N > 0; };  //!< Is calculation enabled? N = 0: no.
    void    sort();
    void    merge(NegfResultList &second);
    void    save(string fileName);
};

}
#endif	/* NEGFRESULT_H */


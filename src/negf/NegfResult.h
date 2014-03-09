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
namespace negf{

using namespace utils::stds;
using namespace maths::armadillo;

typedef pair<double, cxmat> negf_result;     //!< Result as a function of energy.

struct ResultComparator : public binary_function<negf_result&, negf_result&, bool>{
    ResultComparator( double tol = 1e-7 ) : epsilon(tol) {}
    bool operator()( const negf_result &left, const negf_result &right  ) const{
        return (abs(left.first - right.first) > epsilon) && (left.first < right.first);
    }   
    double epsilon;
};

struct NegfResultList{
public:    
    std::list<negf_result>R;              //!< List of results.
    uint            N;              //!< Size of the cxmat matrix. N == 0 means do not calculate.
    string          tag;         //!< Suffix added to the output filename.
    
    NegfResultList(string tag = "", uint N = 0);
    
    bool    isEnabled() { return N > 0; };  //!< Is calculation enabled? N = 0: no.
    void    sort();
    void    merge(NegfResultList &second);
    void    save(ostream &out);
};

}
}

#endif	/* NEGFRESULT_H */


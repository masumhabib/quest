/* 
 * File:   NegfResult.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 13, 2014, 3:41 PM
 */

#ifndef NEGFRESULT_H
#define	NEGFRESULT_H

#include "maths/arma.hpp"
#include "utils/std.hpp"
#include "utils/serialize.hpp"

namespace qmicad{
namespace negf{

using namespace utils::stds;
using namespace maths::armadillo;

/**
 * Result as a function of energy.
 */
typedef struct{
    double E;
    cxmat M;
    
private:
    //! For MPI send/receive.
    friend class boost::serialization::access;

    //! For MPI send/receive.
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & E;
        ar & M;
    }
    
} negf_result;

struct ResultComparator : public binary_function<negf_result&, negf_result&, bool>{
    ResultComparator( double tol = 1e-7 ) : epsilon(tol) {}
    bool operator()( const negf_result &left, const negf_result &right  ) const{
        return (abs(left.E - right.E) > epsilon) && (left.E < right.E);
    }   
    double epsilon;
};

struct NegfResultList{
public:    
    typedef std::list<negf_result>::iterator iter;
    
    std::list<negf_result>R;     //!< List of results.
    int             ib;          //!< Block i.
    int             jb;          //!< Block j.
    uint            N;           //!< Size of the cxmat matrix. N == 0 means do not calculate.
    string          tag;         //!< Suffix added to the output filename.
    
    NegfResultList(string tag = "", uint N = 0, int ib = -1, int jb = -1);
    
    bool    isEnabled() { return N > 0; };  //!< Is calculation enabled? N = 0: no.
    void    sort();
    void    merge(NegfResultList &second);
    void    save(ostream &out, bool isText = true);
};

}
}

#endif	/* NEGFRESULT_H */


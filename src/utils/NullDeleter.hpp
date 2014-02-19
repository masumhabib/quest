/* 
 * File:   nullDeleter.hpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on February 18, 2014, 4:59 PM
 */

#ifndef NULLDELETER_HPP
#define	NULLDELETER_HPP

namespace utils{

struct NullDeleter{
    void operator()(void const *) const
    {
    }
};

}

#endif	/* NULLDELETER_HPP */


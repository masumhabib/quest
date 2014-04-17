/* 
 * File:   nullDeleter.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
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


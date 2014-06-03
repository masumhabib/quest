/* 
 * File:   fermi.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 8:23 PM
 */

#ifndef FERMI_HPP
#define	FERMI_HPP

#include <armadillo>
using arma::exp;

template<typename T>
T fermi(const T& E, double mu, double kT){
    return (1/(1 + exp((E - mu)/kT)));
}

#endif	/* FERMI_HPP */


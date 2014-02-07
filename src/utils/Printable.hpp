/* 
 * File:   Printable.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 18, 2013, 4:27 PM
 * 
 * Description: Printable class. All class that want to print itself
 * extends this class and overloads the toString() method.
 * 
 */

#ifndef PRINTABLE_H
#define	PRINTABLE_H

#include <string>
#include <iostream>
using namespace std;

class Printable{
    protected:
        string mTitle;
        string mPrefix;
        
        Printable(const string &prefix = ""):mPrefix(prefix) {
            mTitle = "Printable";
        };
        virtual ~Printable(){};
        
    public:
        virtual string toString() const { return mTitle; };
        virtual bool isEmpty()const { return false; };

        /* Dump the data to the stream */
        friend ostream& operator << (ostream & out, const Printable &p){
            if(!p.isEmpty()){
                if (p.mTitle.length()){
                    out << p.mPrefix << p.mTitle << ":" << endl;
                }
                out << p.toString();
            }
            return out;
        };
        
        void Prefix(const string& prefix){
            mPrefix = prefix;  
        };
        
        string Prefix(){
            return mPrefix;
        };
        
        void Title(const string &title){
            mTitle = title;
        }
        
        string Title(){
            return mTitle;
        }
    };


#endif	/* PRINTABLE_H */


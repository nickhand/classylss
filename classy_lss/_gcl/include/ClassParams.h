#ifndef CLASS_PARAMS_H
#define CLASS_PARAMS_H

// CLASS
#include "class.h"
#include "output.h"

#include "Common.h"
#include "boost/any.hpp"
#include <utility>
#include <map>
#include <iostream>
#include <sstream>
#include <string>

// function to ouput an option
template<class T>
void output_option(const boost::any* v, std::stringstream& out) {
    out << *(boost::any_cast<T>(v));
}

/*----------------------------------------------------------------------------*/
/* class to encapsulate CLASS parameters from any type (numerical or string)  */
/*----------------------------------------------------------------------------*/
class ClassParams {

public:
    
    typedef std::map<std::string, std::string> param_vector;
    typedef param_vector::iterator iterator;
    typedef param_vector::const_iterator const_iterator;
    
    ClassParams();
    ClassParams(const std::string& param_file); 
    ClassParams(const ClassParams &other);
    
    void update(const ClassParams &other);

    // use this to add a CLASS variable
    int add(const std::string& key, const int& val); 
    int add(const std::string& key, const float& val); 
    int add(const std::string& key, const double& val); 
    int add(const std::string& key, const bool& val); 
    int add(const std::string& key, const std::string& val);
    int add(const std::string& key, const char* val); 
  
    std::vector<std::string> keys() const; 
    void print() const; 
    bool contains(const std::string& key) const;
       
    // accesors
    inline unsigned size() const {return pars.size();}
    
    // return the string representation
    inline const std::string& value(const std::string& key) const {return pars.at(key);}
    
    // iterate over the pars variable
    const_iterator begin() const { return pars.begin(); }
    const_iterator end() const { return pars.end(); }
  
    // overload the [] operator to return const reference
    const std::string& operator[](const std::string& key) const { return this->value(key); }

private:
    
    // a map of string keys and string values
    param_vector pars;
};

#endif
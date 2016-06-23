#include "ClassParams.h" 

#include <numeric>
#include <iomanip>
#include <algorithm>

using namespace std;

ClassParams::ClassParams() {}

// destructor
ClassParams::~ClassParams() {}

ClassParams::ClassParams(const string& param_file) {
      
      // for error messages
      ErrorMsg _errmsg; 
      
      // initialize an empty file_content to read the params file
      struct file_content fc_input;
      fc_input.size = 0;
      fc_input.filename=new char[1];
      
      // read the param file
      if (parser_read_file(const_cast<char*>(param_file.c_str()), &fc_input, _errmsg) == _FAILURE_){
          throw invalid_argument(_errmsg);
      }
      
      // set the pars (all strings)
      for (int i=0; i < fc_input.size; i++) {
          this->add(fc_input.name[i], std::string(fc_input.value[i]));
      }
      
      // free the input
      parser_free(&fc_input);
      
}

ClassParams::ClassParams(const ClassParams &other) 
{
    pars.clear();  
    for (auto it = other.pars.begin(); it != other.pars.end(); ++it)
        pars[it->first] = it->second;
}

void ClassParams::update(const ClassParams &other) 
{    
    for (auto it = other.pars.begin(); it != other.pars.end(); ++it)
        pars[it->first] = it->second;
}


std::vector<std::string> ClassParams::keys() const
{
    vector<std::string> v;
    for (const_iterator it = pars.begin(); it != pars.end(); ++it) 
        v.push_back(it->first);
    return v;
}

void ClassParams::print() const 
{  
    for (const_iterator it = pars.begin(); it != pars.end(); ++it) 
      cout << it->first << " = " << it->second << endl;
}

bool ClassParams::contains(const std::string& key) const
{
    if (pars.size() == 0)
        return false;
     
    return (pars.find(key) != pars.end());
}

int ClassParams::add(const std::string& key, const int val) 
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const long val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const long long val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const unsigned val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const unsigned long val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const unsigned long long val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const float val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const double val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const std::string& key, const long double val)
{
    pars[key] = std::to_string(val); 
    return pars.size();
}

int ClassParams::add(const string& key, const bool val) { 
       
    pars[key] = val ? "yes" : "no"; 
    return pars.size(); 
}

int ClassParams::add(const string& key, const string val) { 
    
    pars[key] = val; 
    return pars.size(); 
}

int ClassParams::add(const string& key, const char * val) { 
    
    pars[key] = string(val); 
    return pars.size(); 
}

string ClassParams::pop(const string& key) 
{
    // use-remove
    auto it = pars.find(key);
    if (it != pars.end()) 
    {
      // move-remove-use
      auto x = std::move(it->second);
      pars.erase(it);
      return x; 
      
    } else
        Common::throw_error("no such key to pop: " + key, __FILE__, __LINE__);

    return "";
}



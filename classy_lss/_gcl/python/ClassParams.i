%{
#include "ClassParams.h"
%}

namespace std {
    %template(VectorString) vector<std::string>;
};
%rename(_update) update(const ClassParams &other);

class ClassParams {
public:    
    
    ClassParams();
    ClassParams(const std::string& param_file);
    
    void update(const ClassParams &other);

    int add(const std::string& key, const int& val); 
    int add(const std::string& key, const float& val); 
    int add(const std::string& key, const double& val); 
    int add(const std::string& key, const bool& val); 
    int add(const std::string& key, const std::string& val);
    int add(const std::string& key, const char* val); 

    std::vector<std::string> keys() const;
    unsigned size() const;
    void print() const; 
    
    const std::string& value(const std::string& key) const;
    
    %pythoncode %{
        
        @classmethod
        def from_dict(cls, d):
            toret = cls()
            for k in d: toret.add(k, d[k])
            return toret
            
        def update(self, *args, **kwargs):
            if len(args):
                if len(args) != 1:
                    raise ValueError("only one positional argument, a dictionary")
                d = args[0]
                if not isinstance(d, dict):
                    raise TypeError("first argument must be a dictionary")
                kwargs.update(d)
                
            if not len(kwargs):
                raise ValueError("no parameters provided to update")
            pars = self.__class__.from_dict(kwargs)
            self._update(pars)
            
        def __getitem__(self, key):
            try: 
                return self.value(key)
            except Exception as e: 
                if "no such key" in str(e):
                    raise KeyError(str(e))
            raise ValueError("no function available to return value for '%s'" %key)
            
        def __setitem__(self, key, value):
            self.add(key, value)
            
        def __iter__(self):
            return iter(self.keys()) 
            
        def values(self): 
            return [self[k] for k in self.keys()]
            
        def items(self):
            return list(zip(self.keys(), self.values()))
            
        def __len__(self):
            return self.size()
            
        def __str__(self):
            return dict(self).__str__()
        
        def __repr__(self):
            return dict(self).__repr__()
    %}
  
};

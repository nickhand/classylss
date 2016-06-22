import sys
import os

valid_keys = ['OPTFLAG', 'OMPFLAG', 'LDFLAG', 'CCFLAG', 'CC', 'AR']
def get_var_from_line(line, d):
    """
    Extract lines from the open file, ignoring everything
    that isn't KEY = VALUE syntax
    """
    toret = None
    if '#' in line:
        line = line[:line.index('#')]
    if not len(line):
        return toret
        
    try:
        fields = [l.strip() for l in line.split("=")]
        if len(fields) == 2 and fields[0] in valid_keys:
            d[fields[0]] = fields[1]
            toret = fields[0]
    except:
        pass
        
    return toret
    
# the user config
with open("class.cfg", "r") as ff:
    
    usr_cfg = {}
    for line in ff.readlines():
        get_var_from_line(line, usr_cfg)
    
# first command line argument is the CLASS root  
classdir = sys.argv[1]

# second is the desired install directory
installdir = sys.argv[2];

make_lines = open(os.path.join(classdir, 'Makefile'), 'r').readlines()

# write a new makefile
with open(os.path.join(classdir, 'Makefile'), 'w') as new_makefile:
    
    make_cfg = {}
    for i, line in enumerate(make_lines):
        
        # update the classdir
        if '__CLASSDIR__' in line:
            make_lines[i] = "CCFLAG += -D__CLASSDIR__='\"%s\"'" %installdir
        
        # update env vars
        key = get_var_from_line(line, make_cfg)
        if key is not None and key in usr_cfg:
            make_lines[i] = "%s = %s\n" %(key, usr_cfg[key])
            
    new_makefile.write("".join(make_lines))
            
            
        
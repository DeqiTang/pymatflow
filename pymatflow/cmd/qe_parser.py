
def read_pwscf_in(filepath):
    """
    Note: read parameters from pwscf input template
    """
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    
    control = {}
    system = {}
    electrons = {}
    ions = {}   
    cell = {}
    
    for i in range(len(lines)):
        if lines[i].split()[0].lower() == "&control":       
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &control variable
                    contorl[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    control[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1
        if lines[i].split()[0].lower() == "&system":          
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &control variable
                    system[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    system[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1
        if lines[i].split()[0].lower() == "&electrons":         
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &control variable
                    electrons[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    electrons[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1                
        if lines[i].split()[0].lower() == "&ions":
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &control variable
                    ions[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    ions[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1                     
        if lines[i].split()[0].lower() == "&cell":
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &control variable
                    cell[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    cell[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1  
                
    return control, system, electrons, ions, cell
    
    
    
def read_neb_in(filepath):
    """
    Note: read parameters from neb.x input template
    """
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    
    path = {}
    
    for i in range(len(lines)):
        if lines[i].split()[0].lower() == "&path":       
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &PATH variable
                    path[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    path[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1    
    return path
    
def read_ph_in(filepath):
    """
    Note: read parameters from neb.x input template
    """
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    
    ph = {}
    
    for i in range(len(lines)):
        if lines[i].split()[0].lower() == "&inputph":       
            j = 1
            while lines[i+j].split()[0] != "/":
                if len(lines[i+j].split()) == 0:
                     pass
                if len(lines[i+j].split("\n")[0].split("#")[0].split("=")) == 2:
                    # in case of single value &INPUTPH variable
                    ph[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                else:
                    ph[lines[i+j].split("=")[0].split()[0]] = lines[i+j].split("\n")[0].split("#")[0].split("=")[1].split()
                j += 1    
    return ph
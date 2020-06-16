
def read_inp(file_path):
    """
    Note: read parameters from cp2k input file
    """
    out = {}
    args_all = [
        "GLOBAL-PRINT_LEVEL",                                                                            
        "FORCE_EVAL-SUBSYS-CELL-SYMMETRY",                                                               
        "FORCE_EVAL-DFT-LS_SCF",                                                                         
        "FORCE_EVAL-DFT-QS-METHOD",                                                                      
        "FORCE_EVAL-DFT-MGRID-CUTOFF",                                                                   
        "FORCE_EVAL-DFT-MGRID-REL_CUTOFF",                                                               
        "FORCE_EVAL-DFT-MGRID-NGRIDS",                                                                   
        "FORCE_EVAL-DFT-XC-XC_FUNCTIONAL",                                                               
        "FORCE_EVAL-DFT-SCF-EPS_SCF",                                                                    
        "FORCE_EVAL-DFT-SCF-ADDED_MOS",                                                                  
        "FORCE_EVAL-DFT-SCF-SMEAR",                                                                      
        "FORCE_EVAL-DFT-SCF-SMEAR-METHOD",                                                               
        "FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE",                                               
        "FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE",                                                          
        "FORCE_EVAL-DFT-SCF-DIAGONALIZATION",                                                            
        "FORCE_EVAL-DFT-SCF-OT",                                                                         
        "FORCE_EVAL-DFT-SCF-MIXING-ALPHA",                                                               
        "FORCE_EVAL-DFT-KPOINTS-SCHEME",                                                                 
        "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE",                                                
        "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE",                                           
        "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF",                                       
        "FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE",                                                          
        "FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE",                                                    
        "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE",                                                
        "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION",                                       
        "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST",                                            
        "MOTION-GEO_OPT-MAX_ITER",                                                                       
        "MOTION-GEO_OPT-OPTIMIZER",                                                                      
        "MOTION-GEO_OPT-TYPE",                                                                           
        "MOTION-GEO_OPT-MAX_DR",                                                                         
        "MOTION-GEO_OPT-MAX_FORCE",                                                                      
        "MOTION-GEO_OPT-RMS_DR",                                                                         
        "MOTION-GEO_OPT-RMS_FORCE",                                                                      
        "MOTION-CELL_OPT-MAX_ITER",                                                                      
        "MOTION-CELL_OPT-OPTIMIZER",                                                                     
        "MOTION-CELL_OPT-TYPE",                                                                          
        "MOTION-CELL_OPT-MAX_DR",                                                                        
        "MOTION-CELL_OPT-MAX_FORCE",                                                                     
        "MOTION-CELL_OPT-RMS_DR",                                                                        
        "MOTION-CELL_OPT-RMS_FORCE",                                                                     
        "MOTION-CELL_OPT-PRESSURE_TOLERANCE",                                                            
        "MOTION-CELL_OPT-KEEP_ANGLES",                                                                   
        "MOTION-CELL_OPT-KEEP_SYMMETRY",                                                                 
        "MOTION-BAND-BAND_TYPE",                                                                         
        "MOTION-BAND-NUMBER_OF_REPLICA",                                                                 
        "MOTION-BAND-ALIGN_FRAMES",                                                                      
        "MOTION-BAND-ROTATE-FRAMES",                                                                     
        "MOTION-BAND-K_SPRING",                                                                          
        "MOTION-MD-STEPS",                                                                              
        "MOTION-MD-TIMESTEP",                                                                            
        "MOTION-MD-ENSEMBLE",                                                                            
        "MOTION-MD-TEMPERATURE",                                                                         
        "MOTION-MD-TEMP_TOL",                                                                            
        "MOTION-PRINT-TRAJECTORY-FORMAT",                                                                
        "VIBRATIONAL_ANALYSIS-DX",                                                                       
        "VIBRATIONAL_ANALYSIS-FULLY_PERIODIC",                                                           
        "VIBRATIONAL_ANALYSIS-INTENSITIES",                                                              
        "VIBRATIONAL_ANALYSIS-TC_PRESSURE",                                                              
        "VIBRATIONAL_ANALYSIS-TC_TEMPERATURE",                                                           
        "VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY",
    ]
    with open(filepath) as fin:
        lines = fin.readlines()
    for arg in args_all:
        # len(arg.split("-")) == 2
        if len(arg.split("-")) == 2:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == arg.split("-")[1] or lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            out[arg] = lines[i+j].split()[1].upper()
                            break
                        j += 1
                    break
        # len(arg.split("-")) == 3
        if len(arg.split("-")) == 3:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == arg.split("-")[2] or lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    out[arg] = lines[i+j].split()[2].upper()
                                    break
                                k += 1
                            break
                        j += 1
                    break        
        # len(arg.split("-")) == 4
        if len(arg.split("-")) == 4:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == arg.split("-")[3] or lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            out[arg] = lines[i+j].split()[3].upper()
                                            break
                                        l += 1
                                    break
                                k += 1
                            break
                        j += 1
                    break     
        # len(arg.split("-")) == 5
        if len(arg.split("-")) == 5:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == arg.split("-")[4] or lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    out[arg] = lines[i+j].split()[4].upper()
                                                    break
                                                m += 1
                                            break
                                        l += 1
                                    break
                                k += 1
                            break
                        j += 1
                    break            
        # len(arg.split("-")) == 6
        if len(arg.split("-")) == 6:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    n = 1
                                                    while lines[i+j+k+l+m+n].split()[0].upper() != "&END %s" (arg.split("-")[4]):
                                                        if lines[i+j+k+l+m+n].split()[0].upper() == arg.split("-")[5] or lines[i+j+k+l+m+n].split()[0].upper() == ("&" + arg.split("-")[5]):
                                                            out[arg] = lines[i+j].split()[5].upper()
                                                            break
                                                        n += 1
                                                    break
                                                m += 1
                                            break
                                        l += 1
                                    break
                                k += 1
                            break
                        j += 1
                    break                             
        # len(arg.split("-")) == 7
        if len(arg.split("-")) == 7:
            for i in range(len(lines)):
                if lines[i].split()[0].upper() == arg.split("-")[0]:
                    j = 1
                    while lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    n = 1
                                                    while lines[i+j+k+l+m+n].split()[0].upper() != "&END %s" (arg.split("-")[4]):
                                                        if lines[i+j+k+l+m+n].split()[0].upper() == ("&" + arg.split("-")[5]):
                                                            o = 1
                                                            while lines[i+j+k+l+m+n+o].split()[0].upper() != "&END %s" % (arg.split("-")[5]):
                                                                if lines[i+j+k+l+m+n+o].split()[0].upper() == arg.split("-")[6] or lines[i+j+k+l+m+n+o].split()[0].upper() == ("&" + arg.split("-")[6]):
                                                                    out[arg] = lines[i+j].split()[6].upper()
                                                                    break
                                                                o += 1
                                                            break
                                                        n += 1
                                                    break
                                                m += 1
                                            break
                                        l += 1
                                    break
                                k += 1
                            break
                        j += 1
                    break                       

def read_inp_old(filepath):
    """
    Note: read parameters from cp2k input file
    """
    out = {}
    with open(filepath) as fin:
        lines = fin.readlines() 
    
    for i in range(len(lines)):
        
        # GLOBAL
        if lines[i].split()[i].upper() == "&GLOBAL":
            j = 1
            while lines[i+j].split()[0] != "&END GLOBAL"
                if lines[i].split()[0].upper() == "PRINT_LEVEL":
                    out["GLOBAL-PRINT_LEVEL"] = lines[i].split()[1]
                j += 1
            continue
        
        # FORCE_EVAL
        if lines[i].split()[i].upper() == "&FORCE_EVAL":
            # FORCE_EVAL-SUBSYS-CELL
            j = 1
            while lines[i+j].split()[0].upper() != "&END CELL":
                if lines[i+j].split()[0].upper() == "SYMMETRY":
                    out["FORCE_EVAL-SUBSYS-CELL-SYMMETRY"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-LS_SCF
            j = 1
            while lines[i+j].split()[0].upper() != "&END LS_SCF":
                if lines[i+j].split()[0].upper() == "&LS_SCF":
                    out["FORCE_EVAL-DFT-LS_SCF"] = lines[i+j].split()[1]
                    break
                j += 1
            continue
            # FORCE_EVAL_DFT-QS-METHOD
            j = 1 
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&QS":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END QS":
                        if lines[i+j+k].split()[0].upper() == "METHOD":
                            out["FORCE_EVAL-DFT-QS-METHOD"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-MGRID-CUTOFF
            j = 1 
            while lines[i+j].split()[0].upper() != "&END MGRID":
                if lines[i+j].split()[0].upper() == "CUTOFF":
                    out["FORCE_EVAL-DFT-MGRID-CUTOFF"] = lines[i+j].split()[1]
                    break
                j += 1 
            # FORCE_EVAL-DFT-MGRID-REL_CUTOFF 
            j = 1 
            while lines[i+j].split()[0].upper() != "&END MGRID":
                if lines[i+j].split()[0].upper() == "REL_CUTOFF":  
                    out["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-MGRID-NGRIDS
            j = 1 
            while lines[i+j].split()[0].upper() != "&END MGRID":
                if lines[i+j].split()[0].upper() == "NGRIDS":
                    out["FORCE_EVAL-DFT-MGRID-NGRIDS"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-XC-XC_FUNCITONAL
            j = 1
            while lines[i+j].split()[0].upper() != "&END XC":
                if lines[i+j].split()[0].upper() == "&XC_FUNCTIONAL":
                    out["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = lines[i+j].split()[1]
                    break
                j += 1 
            # FORCE_EVAL-DFT-SCF-EPS_SCF
            j = 1
            while lines[i+j].split()[0].upper() != "&END SCF":
                if lines[i+j].split()[0].upper() == "EPS_SCF":
                    out["FORCE_EVAL-DFT-SCF-EPS_SCF"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-ADDED_MOS
            j = 1 
            while lines[i+j].split()[0].upper() != "&END SCF":
                if lines[i+j].split()[0].upper() == "ADDED_MOS":
                    out["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = line[i+j].split()[1]
                    break
                j += 1 
            # FORCE_EVAL-DFT-SCF-SMEAR
            j = 1 
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&SMEAR":
                    out["FORCE_EVAL-DFT-SCF-SMEAR"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-SMEAR-METHOD
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&SMEAR":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END SMEAR":
                        if lines[i+j+k].split()[0].upper() == "METHOD":
                            out["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&SMEAR":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END SMEAR":
                        if lines[i+j+k].split()[0].upper() == "ELECTRONIC_TEMPERATURE":
                            out["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&SMEAR":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END SMEAR":
                        if lines[i+j+k].split()[0].upper() == "WINDOW_SIZE":
                            out["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-DIAGONALIZATION
            j = 1 
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&DIAGONALIZATION":
                    out["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-OT
            j = 1 
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&OT":
                    out["FORCE_EVAL-DFT-SCF-OT"] = lines[i+j].split()[1]
                    break
                j += 1
            # FORCE_EVAL-DFT-SCF-MIXING-ALPHA
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&MIXING":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MIXING":
                        if lines[i+j+k].split()[0].upper() == "ALPHA":
                            out["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-KPOINTS-SCHEME
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&KPOINS":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END KPOINTS":
                        if lines[i+j+k].split()[0].upper() == "SCHEME":
                            out["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1            
            # FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&VDW_POTENTIAL":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END VDW_POTENTIAL":
                        if lines[i+j+k].split()[0].upper() == "POTENTIAL_TYPE":
                            out["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&VDW_POTENTIAL":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END VDW_POTENTIAL":
                        if lines[i+j+k].split()[0].upper() == "&PAIR_POTENTIAL":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END PAIR_POTENTIAL":
                                if lines[i+j+k+m].split()[0].upper() == "TYPE":
                                    out["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = lines[i+j+k+m].split()[1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1            
            # FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&VDW_POTENTIAL":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END VDW_POTENTIAL":
                        if lines[i+j+k].split()[0].upper() == "&PAIR_POTENTIAL":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END PAIR_POTENTIAL":
                                if lines[i+j+k+m].split()[0].upper() == "R_CUTOFF"
                                    out["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF"] = lines[i+j+k+m].split()[1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1                       
            # FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&ELF_CUBE":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END ELF_CUBE":
                        if lines[i+j+k].split()[0].upper() == "STRIDE":
                            out["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1
            # FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&E_DENSITY_CUBE":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END E_DENSITY_CUBE":
                        if lines[i+j+k].split()[0].upper() == "STRIDE":
                            out["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1            
            # FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&RESP":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END RESP":
                        if lines[i+j+k].split()[0].upper() == "&SLAB_SAMPLING":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END SLAB_SAMPLING":
                                if lines[i+j+k+m].split()[0].upper() == "RANGE"
                                    out["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"] = lines[i+j+k+m].split()[1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1    
            # FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&RESP":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END RESP":
                        if lines[i+j+k].split()[0].upper() == "&SLAB_SAMPLING":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END SLAB_SAMPLING":
                                if lines[i+j+k+m].split()[0].upper() == "SURF_DIRECTION"
                                    out["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"] = lines[i+j+k+m].split()[1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1                    
            # FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST
            j = 1
            while lines[i+j].split()[0].upper() != "&END FORCE_EVAL":
                if lines[i+j].split()[0].upper() == "&RESP":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END RESP":
                        if lines[i+j+k].split()[0].upper() == "&SLAB_SAMPLING":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END SLAB_SAMPLING":
                                if lines[i+j+k+m].split()[0].upper() == "ATOM_LIST"
                                    out["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"] = lines[i+j+k+m].split("\n")[0].split()[1:-1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1    
            continue
        
        # MOTION
        if lines[i].split()[i].upper() == "&MOTION":
            # MOTION-GEO_OPT-MAXITER
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAXITER":
                            out["MOTION-GEO_OPT-MAX_ITER"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1     
            # MOTION-GEO_OPT-OPTIMIZER                
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "OPTIMIZER":
                            out["MOTION-GEO_OPT-OPTIMIZER"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                 
            # MOTION-GEO_OPT-TYPE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "TYPE":
                            out["MOTION-GEO_OPT-TYPE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                 
            # MOTION-GEO_OPT-MAX_DR
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAX_DR":
                            out["MOTION-GEO_OPT-MAX_DR"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1     
            # MOTION-GEO_OPT-MAX_FORCE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAX_FORCE":
                            out["MOTION-GEO_OPT-MAX_FORCE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1     
            # MOTION-GEO_OPT-RMS_DR
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "RMS_DR":
                            out["MOTION-GEO_OPT-RMS_DR"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                 
            # MOTION-GEO_OPT-RMS_FORCE 
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&GEO_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END GEO_OPT":
                        if lines[i+j+k].split()[0].upper() == "RMS_FORCE":
                            out["MOTION-GEO_OPT-RMS_FORCE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                 
            # MOTION-CELL_OPT-MAX_ITER
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAX_ITER":
                            out["MOTION-CELL_OPT-MAX_ITER"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                     
            # MOTION-CELL_OPT-OPTIMIZER
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "OPTIMIZER":
                            out["MOTION-CELL_OPT-OPTIMIZER"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-TYPE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "TYPE":
                            out["MOTION-CELL_OPT-TYPE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1             
            # MOTION-CELL_OPT-MAX_DR
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAC_DR":
                            out["MOTION-CELL_OPT-MAX_DR"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-MAX_FORCE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "MAX_ITER":
                            out["MOTION-CELL_OPT-MAX_ITER"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1             
            # MOTION-CELL_OPT-RMS_DR
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "RMS_DR":
                            out["MOTION-CELL_OPT-RMS_DR"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-RMS_FORCE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "RMS_FORCE":
                            out["MOTION-CELL_OPT-RMS_FORCE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-PRESSURE_TOLLERANCE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "PRESSURE_TOLLERANCE":
                            out["MOTION-CELL_OPT-PRESSURE_TOLERANCE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-KEEP_ANGLES
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "KEEP_ANGLES":
                            out["MOTION-CELL_OPT-KEEP_ANGLES"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-CELL_OPT-KEEP_SYMMETRY
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&CELL_OPT":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END CELL_OPT":
                        if lines[i+j+k].split()[0].upper() == "KEEP_SYMMETRY":
                            out["MOTION-CELL_OPT-KEEP_SYMMETRY"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-BAND-BAND_TYPE
            j = 1
            while lines[i+j].split()[0].upper() != "&END BAND":
                if lines[i+j].split()[0].upper() == "&BAND":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END BAND":
                        if lines[i+j+k].split()[0].upper() == "BAND_TYPE":
                            out["MOTION-BAND-BAND_TYPE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         
            # MOTION-BAND-NUMBER_OF_REPLICA
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&BAND":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END BAND":
                        if lines[i+j+k].split()[0].upper() == "NUMBER_OF_REPLICA":
                            out["MOTION-BAND-NUMBER_OF_REPLICA"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1   
            # MOTION-BAND-ALIGN_FRAMES
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&BAND":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END BAND":
                        if lines[i+j+k].split()[0].upper() == "ALIGN_FRAMES":
                            out["MOTION-BAND-ALIGN_FRAMES"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1               
            # MOTION-BAND-ROTATE-FRAMES
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&BAND":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END BAND":
                        if lines[i+j+k].split()[0].upper() == "&ROTATE":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END ROTATE":
                                if lines[i+j+k+m].split()[0].upper() == "FRAMES"
                                    out["MOTION-BAND-ROTATE-FRAMES"] = lines[i+j+k+m].split()[1]
                                    break
                                m += 1
                            break
                        k += 1
                    break
                j += 1                
            # MOTION-BAND-K_SPRING
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&BAND":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END BAND":
                        if lines[i+j+k].split()[0].upper() == "K_SPRING":
                            out["MOTION-BAND-K_SPRING"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                         

            # MOTION-MD-STEPS
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "STEPS":
                            out["MOTION-MD-STEPS"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1                   
            # MOTION-MD-TIMESTEP
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "TIMESTEP":
                            out["MOTION-MD-TIMESTEP"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1       
            # MOTION-MD-ENSEMBLE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "ENSEMBLE":
                            out["MOTION-MD-ENSEMBLE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1               
            # MOTION-MD-TEMPERATURE
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "TEMPERATURE":
                            out["MOTION-MD-TEMPERATURE"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1               
            # MOTION-MD-TEMP_TOL
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "TEMP_TOL":
                            out["MOTION-MD-TEMP_TOL"] = lines[i+j+k].split()[1]
                            break
                        k += 1
                    break
                j += 1               
            # MOTION-MD-PRINT-TRAJECTORY-FORMAT
            j = 1
            while lines[i+j].split()[0].upper() != "&END MOTION":
                if lines[i+j].split()[0].upper() == "&MD":
                    k = 1 
                    while lines[i+j+k].split()[0].upper() != "&END MD":
                        if lines[i+j+k].split()[0].upper() == "&PRINT":
                            m = 1
                            while lines[i+j+k+m].split()[0].upper() != "&END PRINT":
                                if lines[i+j+k+m].split()[0].upper() == "&TRAJECTORY":
                                    l = 1
                                    while lines[i+j+k+m+l].split()[0].upper() != "&END TRAJECTORY":
                                        if lines[i+j+k+m+l].split()[0].upper() == "FORMAT":
                                            out["MOTION-MD-PRINT-TRAJECTORY-FORMAT"] = lines[i+j+k+m].split()[1]
                                            break
                                        l += 1
                                m += 1
                            break
                        k += 1
                    break
                j += 1                    
                
            continue
        
        # VIBRATIONAL_ANALYSIS
        if lines[i].split()[0].upper() == "&VIBRATIONAL_ANALYSIS":
            # VIBRATIONAL_ANALYSIS-DX
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "DX":
                    out["VIBRATIONAL_ANALYSIS-DX"] = lines[i].split()[1]
                    break
                j += 1
            # VIBRATIONAL_ANALYSIS-FULLY_PERIODIC
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "FULLY_PERIODIC":
                    out["VIBRATIONAL_ANALYSIS-FULLY_PERIODIC"] = lines[i].split()[1]
                    break
                j += 1            
            # VIBRATIONAL_ANALYSIS-INTENSITIES
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "INTENSITIES":
                    out["VIBRATIONAL_ANALYSIS-INTENSITIES"] = lines[i].split()[1]
                    break
                j += 1                
            # VIBRATIONAL_ANALYSIS-TC_PRESSURE
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "TC_PRESSURE":
                    out["VIBRATIONAL_ANALYSIS-TC_PRESSURE"] = lines[i].split()[1]
                    break
                j += 1            
            # VIBRATIONAL_ANALYSIS-DX
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "TC_TEMPERATURE":
                    out["VIBRATIONAL_ANALYSIS-TC_TEMPERATURE"] = lines[i].split()[1]
                    break
                j += 1            
            # VIBRATIONAL_ANALYSIS-DX
            j = 1
            while lines[i+j].split()[0].upper() != "&END VIBRATIONAL_ANALYSIS":
                if lines[i].split()[0].upper() == "THERMOCHEMISTRY":
                    out["VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY"] = lines[i].split()[1]
                    break
                j += 1                
            continue
    return out
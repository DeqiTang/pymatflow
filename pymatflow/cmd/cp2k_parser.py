
def read_inp(filepath):
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
                if len(lines[i].split()) == 0:
                    continue
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == arg.split("-")[1] or lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            out[arg] = lines[i+j].split()[1].upper()
                            break
                        j += 1
                    break
        # len(arg.split("-")) == 3
        if len(arg.split("-")) == 3:
            for i in range(len(lines)):
                if len(lines[i].split()) == 0:
                    continue
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while len(lines[i+j+k].split()) > 0 and lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == arg.split("-")[2] or lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    out[arg] = lines[i+j+k].split()[1].upper()
                                    break
                                k += 1
                            break
                        j += 1
                    break        
        # len(arg.split("-")) == 4
        if len(arg.split("-")) == 4:
            for i in range(len(lines)):
                if len(lines[i].split()) == 0:
                    continue            
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while len(lines[i+j+k].split()) > 0 and lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while len(lines[i+j+k+l].split()) > 0 and lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == arg.split("-")[3] or lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            out[arg] = lines[i+j+k+l].split()[1].upper()
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
                if len(lines[i].split()) == 0:
                    continue            
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while len(lines[i+j+k].split()) > 0 and lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while len(lines[i+j+k+l].split()) > 0 and lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while len(lines[i+j+k+l+m]) > 0 and lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == arg.split("-")[4] or lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    out[arg] = lines[i+j+k+l+m].split()[1].upper()
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
                if len(lines[i].split()) == 0:
                    continue            
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while len(lines[i+j+k].split()) > 0 and lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while len(lines[i+j+k+l].split()) > 0 and lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while len(lines[i+j+k+l+m].split()) > 0 and lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    n = 1
                                                    while len(lines[i+j+k+l+m+n].split()) > 0 and lines[i+j+k+l+m+n].split()[0].upper() != "&END %s" (arg.split("-")[4]):
                                                        if lines[i+j+k+l+m+n].split()[0].upper() == arg.split("-")[5] or lines[i+j+k+l+m+n].split()[0].upper() == ("&" + arg.split("-")[5]):
                                                            out[arg] = lines[i+j+k+l+m+n].split()[1].upper()
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
                if len(lines[i].split()) == 0:
                    continue            
                if lines[i].split()[0].upper() == arg.split("-")[0] or lines[i].split()[0].upper() == "&" + arg.split("-")[0]:
                    j = 1
                    while len(lines[i+j].split()) > 0 and lines[i+j].split()[0].upper() != "&END %s" % (arg.split("-")[0]):
                        if lines[i+j].split()[0].upper() == ("&" + arg.split("-")[1]):
                            k = 1 
                            while len(lines[i+j+k].split()) > 0 and lines[i+j+k].split()[0].upper() != "&END %s" % (arg.split("-")[1]):
                                if lines[i+j+k].split()[0].upper() == ("&" + arg.split("-")[2]):
                                    l = 1
                                    while len(lines[i+j+k+l].split()) > 0 and lines[i+j+k+l].split()[0].upper() != "&END %s" % (arg.split("-")[2]):
                                        if lines[i+j+k+l].split()[0].upper() == ("&" + arg.split("-")[3]):
                                            m = 1
                                            while len(lines[i+j+k+l+m].split()) > 0 and lines[i+j+k+l+m].split()[0].upper() != "&END %s" % (arg.split("-")[3]):
                                                if lines[i+j+k+l+m].split()[0].upper() == ("&" + arg.split("-")[4]):
                                                    n = 1
                                                    while len(lines[i+j+k+l+m+n].split()) > 0 and lines[i+j+k+l+m+n].split()[0].upper() != "&END %s" (arg.split("-")[4]):
                                                        if lines[i+j+k+l+m+n].split()[0].upper() == ("&" + arg.split("-")[5]):
                                                            o = 1
                                                            while len(lines[i+j+k+l+m+n+o].split()) > 0 and lines[i+j+k+l+m+n+o].split()[0].upper() != "&END %s" % (arg.split("-")[5]):
                                                                if lines[i+j+k+l+m+n+o].split()[0].upper() == arg.split("-")[6] or lines[i+j+k+l+m+n+o].split()[0].upper() == ("&" + arg.split("-")[6]):
                                                                    out[arg] = lines[i+j+k+l+m+n+o].split()[1].upper()
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
    return out
# end read_inp

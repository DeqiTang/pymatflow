pwscf_character_var = [
    "prefix", "outdir", "pseudo_dir", "calculation",
    "occupations", "smearing", "vdw_corr", "U_projection_type", 
    "input_dft", "exxciv_teatment", "diagonalization", "ion_dynamics",
    "pot_extrapolation", "wfc_extrapolation", "ion_temperature", 
    "cell_dofree", ""
]

def qe_variable_to_string(variable):
    if False == variable.status:
        return ""
    
    if None == variable.value:
        return ""
    
    out = ""

    if 0 == len(variable.value):
        return out + variable.key

    special = ["starting_magnetization", "Hubbard_U", "Hubbard_J0", 
        "Hubbard_alpha", "Hubbard_beta"]
    if variable.key in special:
        for i in range(len(variable.as_val(t=float, dime=1))):
            out += "%s(%d) = %f\n" % (variable.key, i+1, variable.as_val(t=float, dim=1)[i])
        return out
    
    if variable.key in pwscf_character_var:
        out += "%s = '%s'\n" % (variable.key, variable.as_val(t=str, dim=0))
        return out


    if 1 == len(variable.value):
        if 1 == len(variable.value[0]):
            out += variable.key + " = " + variable.value[0][0]
        else:
            out += variable.key + " ="
            for item in variable.value[0]:
                out += " " + item
    else:
        out += variable.key + " ="
        for val in variable.value[0]:
            out += " " + val
        
        out += "\n"
        for row in range(1, len(variable.value)-1):
            for val in variable.value[row]:
                out += " " + val
            out += "\n"
        for val in variable.value[len(variable.value) - 1]:
            out += " " + val
    return out
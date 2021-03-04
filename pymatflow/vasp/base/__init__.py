
#
def vasp_variable_to_string(variable):
    if False == variable.status:
        return ""
    
    if None == variable.value:
        return ""
    
    out = ""

    if 0 == len(variable.value):
        return out + variable.key

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
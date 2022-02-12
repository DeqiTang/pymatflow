from .section import Cp2kSection

def gen_global():
    out = Cp2kSection("global")
    out.set_param("project", "cp2k_job")
    out.set_param('print_level', "low")
    out.set_param('run_type', "energy_force")
    return out 
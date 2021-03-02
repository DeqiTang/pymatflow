from .section import Cp2kSection


def dft():
    out = Cp2kSection("dft")
    
    out.set_param("basis_set_file_name", "basis_set")
    out.set_param("potential_file_name", "gth_poetntials")

    qs = out.add_subsection("qs")
    qs.set_param("eps_default", "1.0e-14")
    
    mgrid = out.add_subsection("mgrid")
    mgrid.set_param("ngrids", 4)
    mgrid.set_param("cutoff", 100)
    mgrid.set_param("rel_cutoff", 60)
    
    xc = out.add_subsection("xc")
    xc_functional = xc.add_subsection("xc_functional")
    xc_functional.section_parameter = "pbe"

    scf = out.add_subsection("scf")
    scf.set_param("scf_guess", "atomic")
    scf.set_param("eps_scf", "1.0e-7")
    scf.set_param("max_scf", 100)
    diag = scf.add_subsection("diagonalization")
    diag.set_param("algorithm", "standard")

    mixing = out.add_subsection("mixing")
    #mixing.value = "t"
    mixing.set_param("method", "broyden_mixing")
    mixing.set_param("alpha", 0.4)
    mixing.set_param("nbroyden", 8)

    return out

def force_eval():
    out = Cp2kSection("force_eval")

    out.set_param("method", "quickstep")

    out.add_subsection("dft", dft())

    _print = out.add_subsection("print")
    _print.add_subsection("forces").section_parameter = "on"

    return out

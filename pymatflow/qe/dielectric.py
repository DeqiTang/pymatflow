
import os
import re
import sys
import shutil


#from pymatflow.qe.base.control import qe_control
#from pymatflow.qe.base.system import qe_system
#from pymatflow.qe.base.electrons import qe_electrons
#from pymatflow.qe.base.arts import qe_arts
from pymatflow.qe.static import StaticRun



def dielectric_pw(xyz_f, directory="tmp-qe-static", mpi="", runopt="gen", auto=0,
        control={}, system={}, electrons={}, kpoints_option="automatic", kpoints_mp=[2, 2, 2, 0, 0, 0]):
    """
    directory: a place for all the generated files

    parameters:
        directory: the overall static calculation directory

        runopt: determine whether the calculation is executed.
            there are three values: 'gen', 'genrun', 'run'
            'gen': only generate the input files
            'genrun': generate input files and run
            'run': run from the previously generated input files
    Note:
        reference:
            https://github.com/QEF/q-e/tree/master/PW/examples/example10
        here we always use automatic type kpoints

        Berry Phase/electric fields only for insulators!
        so we should try occupations =  'tetrahedra'
    """


    qe = static_run()
    qe.get_xyz(xyz_f)

    if runopt == 'gen' or runopt == 'genrun':
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        #os.system("cp *.UPF %s/" % directory)
        #os.system("cp %s %s/" % (qe.arts.xyz.file, directory))

        # do not copy too many files at the same time or it will be slow
        # so we do not copy all UPF files in the directory but just copy
        # those used in the calculation.
        shutil.copyfile(qe.arts.xyz.file, os.path.join(directory, os.path.basename(qe.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        #

    qe.control.params["lelfield"] = True
    qe.electrons.params["startingwfc"] = 'random'


    # 0 electric field calculation
    qe.control.params["nberrycyc"] = 1
    qe.electrons.params["efield_cart(1)"] = 0
    qe.electrons.params["efield_cart(2)"] = 0
    qe.electrons.params["efield_cart(3)"] = 0
    qe.scf(directory=os.path.join(directory, "0-efield"), runopt=runopt, mpi=mpi, control=control, system=system, electrons=electrons, kpoints_option=kpoints_option, kpoints_mp=kpoints_mp)

    # finite electric field calculation
    qe.control.params["nberrycyc"] = 3
    qe.electrons.params["efield_cart(1)"] = 0
    qe.electrons.params["efield_cart(2)"] = 0
    qe.electrons.params["efield_cart(3)"] = 0.001

    qe.scf(directory=os.path.join(directory, "finite-efield"), runopt=runopt, mpi=mpi, control=control, system=system, electrons=electrons, kpoints_option=kpoints_option, kpoints_mp=kpoints_mp)


    # analysis

"""
Overall abstraction of CP2K
"""
import numpy as np
import sys
import os
import shutil
from pymatflow.remote.server import server_handle

from pymatflow.base.xyz import BaseXyz

from .section import Cp2kSection
from .gen_force_eval import gen_force_eval
from .gen_global import gen_global


"""
"""

class Cp2k:
    """
    Philosophy:
        Implementation of cp2k.base.xxx and cp2l.Cp2k is actually
        not very beautiful. So this is another version of Cp2k class agent
    
    Note:

    """
    def __init__(self):
        """
        """
        self.sections = {}
        self.xyz = BaseXyz()

        self.run_params = {}
        self.set_run()

        
    def to_string(self, indent="\t"):
        out = ""
        for item in self.sections:
            if True == self.sections[item].status:
                out += self.sections[item].to_string(indent=indent) + "\n"
        return out

    def get_xyz(self, xyzfile):
        """
        :param xyzfile:
            a modified xyz formatted file(the second line specifies the cell of the
            system).
        """
        self.xyz.get_xyz(xyzfile)
        self.set_subsys(xyz=self.xyz)

    def _with_bug_force_section_exists_recursive(self, section):
        """
        :param section: in format like this ->
            'force_eval-dft-xc'
        Note:
            check a subsection exist and if not add it
            Recursive version
        """
        def subsection_recursive_check(section_obj, subsection):
            item_list = subsection.lower().split("-")
            if 0 == len(item_list):
                return
            if not item_list[0] in section_obj.subsections:
                section_obj.subsections[item_list[0]] = Cp2kSection(name=item_list[0])
            subsection_recursive_check(section_obj.subsections[item_list[0]], "-".join(item_list[1:]))

        item_list = section.lower().split("-")
        if not item_list[0] in self.sections:
            self.sections[item_list[0]] = Cp2kSection(name=item_list[0])
        subsection_recursive_check(self.sections[item_list[0]], "-".join(item_list[1:]))

    def force_section_exists(self, section):
        """
        :param section: in format like this ->
            'force_eval-dft-xc'
        Note:
            check a subsection exist and if not add it

            we may implement this function recursively.
            however, we choose not to, since the length of number nested subsection is limited
            actually, so by doncitional check in limit range is acceptable.
        """
        item_list = section.lower().split("-")
        for i in range(len(item_list)):
            if 0 == i:
                if not item_list[i] in self.sections:
                    self.sections[item_list[i]] = Cp2kSection(name=item_list[i])
            elif 1 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections:
                    self.sections[item_list[0]].add_subsection(key=item_list[i])
            elif 2 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].add_subsection(key=item_list[i])
            elif 3 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].add_subsection(key=item_list[i])
            elif 4 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].add_subsection(key=item_list[i])
            elif 5 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].add_subsection(key=item_list[i])
            elif 6 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].add_subsection(key=item_list[i])
            elif 7 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].add_subsection(key=item_list[i])
            elif 8 == i:
                if not item_list[i] in self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections[item_list[7]].subsections:
                    self.sections[item_list[0]].subsections[item_list[1]].subsections[item_list[2]].subsections[item_list[3]].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections[item_list[7]].add_subsection(key=item_list[i])
            else:
                print("=========================================================\n")
                print("                     Warning\n")
                print("cp2k.cp2k_dev.Cp2k.force_section_exists: lenth of section string exceeding\n")
                print("limit -> 8\n")
                sys.exit(1)

    def set_section_parameters(self, params={}):
        """
        Note:
            however, we choose not to, since the length of number nested subsection is limited
            actually, so by doncitional check in limit range is acceptable.
        """
        for item in params:
            item_list = item.lower().split("-")
            self.force_section_exists(item)
            if 2 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].section_parameter = params[item]
            elif 3 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].section_parameter = params[item]
            elif 4 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].section_parameter = params[item]
            elif 5 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].section_parameter = params[item]
            elif 6 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].section_parameter = params[item]
            elif 7 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].section_parameter = params[item]
            elif 8 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections[item_list[7]].section_parameter = params[item]
            else:
                print("=========================================================\n")
                print("                     Warning\n")
                print("cp2k.cp2k_dev.Cp2k.set_param: lenth of item in params exceeding\n")
                print("limit -> 8\n")
                sys.exit(1)

    def set_section_vars(self, params={}):
        """
        Note:
            
            however, we choose not to, since the length of number nested subsection is limited
            actually, so by doncitional check in limit range is acceptable.
        """        
        for item in params:
            item_list = item.lower().split("-")
            self.force_section_exists(item)
            if 2 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].section_var.set(key="", value=params[item])
            elif 3 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].section_var.set(key="", value=params[item])
            elif 4 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].section_var.set(key="", value=params[item])
            elif 5 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].section_var.set(key="", value=params[item])
            elif 6 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].section_var.set(key="", value=params[item])
            elif 7 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].section_var.set(key="", value=params[item])
            elif 8 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections[item_list[7]].section_var.set(key="", value=params[item])
            else:
                print("=========================================================\n")
                print("                     Warning\n")
                print("cp2k.cp2k_dev.Cp2k.set_param: lenth of item in params exceeding\n")
                print("limit -> 8\n")
                sys.exit(1)

    def set_params(self, params={}):
        """
        Note:
            we should always use this function to set params in cp2k

            every item in params begin with a "XXX-" where XXX is the first level
            input section of cp2k, like "FORCE_EVAL", "ATOM", "MOTION".

            this function only set params in subsection, so the last item in item_list
            must be a param in subsection but not a subsection.

            however, we choose not to, since the length of number nested subsection is limited
            actually, so by doncitional check in limit range is acceptable.

        Important:
            if we have duplicate subsection name or variable key in a section or subsection(which
            is allowed in cp2k input file), we add a string to the end of the section.name or 
            variable.key "@%d", where "%d" means a digital number. This is used to distinguish 
            possible section name or variable key with same name.
        """
        
        for item in params:
            item_list = item.lower().split("-")
            self.force_section_exists("-".join(item_list[:-1]))
            if 2 == len(item_list):
                self.sections[item_list[0].lower()].set_param(item[-1].lower(), params[item])
            elif 3 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].set_param(item_list[-1].lower(), params[item])
            elif 4 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].set_param(item_list[-1].lower(), params[item])
            elif 5 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].set_param(item_list[-1].lower(), params[item])    
            elif 6 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].set_param(item_list[-1].lower(), params[item])
            elif 7 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].set_param(item_list[-1].lower(), params[item])
            elif 8 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].set_param(item_list[-1].lower(), params[item])
            else:
                print("=========================================================\n")
                print("                     Warning\n")
                print("cp2k.cp2k_dev.Cp2k.set_param: lenth of item in params exceeding\n")
                print("limit -> 8\n")
                sys.exit(1)

    def set_section_status(self, params):
        """
        Note:
            however, we choose not to, since the length of number nested subsection is limited
            actually, so by doncitional check in limit range is acceptable.
        """
        for item in params:
            item_list = item.lower().split("-")
            self.force_section_exists(item)
            if 2 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].status = params[item]
            elif 3 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].status = params[item]
            elif 4 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].status = params[item]
            elif 5 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].status = params[item]
            elif 6 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].status = params[item]
            elif 7 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].status = params[item]
            elif 8 == len(item_list):
                self.sections[item_list[0].lower()].subsections[item_list[1].lower()].subsections[item_list[2].lower()].subsections[item_list[3].lower()].subsections[item_list[4]].subsections[item_list[5]].subsections[item_list[6]].subsections[item_list[7]].status = params[item]
            else:
                print("=========================================================\n")
                print("                     Warning\n")
                print("cp2k.cp2k_dev.Cp2k.set_param: lenth of item in params exceeding\n")
                print("limit -> 8\n")
                sys.exit(1)

    def set_band(self, kpath):
        # See Version1 Cp2k.force_eval.dft.printout.band_structure.set_band(kpath=kpath)
        # to implement this function
        self.kpath = kpath
        for i in range(len(self.kpath) - 1):
            if self.kpath[i][4] != "|" and type(self.kpath[i][4] == int):
                self.set_params({
                    "force_eval-dft-print-band_structure-kpoint_set@%d-units" % (i): "B_VECTOR",
                    "force_eval-dft-print-band_structure-kpoint_set@%d-special_point@%d" % (i, 0): [self.kpath[i][3], str(self.kpath[i][0]), str(self.kpath[i][1]), str(self.kpath[i][2])],
                    "force_eval-dft-print-band_structure-kpoint_set@%d-special_point@%d" % (i, 1): [self.kpath[i+1][3], str(self.kpath[i+1][0]), str(self.kpath[i+1][1]), str(self.kpath[i+1][2])],
                    "force_eval-dft-print-band_structure-kpoint_set@%d-npoints" % (i): self.kpath[i][4] - 1,
                })

    def set_subsys(self, xyz):
        """
        :param xyz:
            an object of pymatflow.base.BaseXyz
        """
        self.xyz = xyz

        subsys = self.sections["force_eval"].add_subsection("subsys")
        cell = subsys.add_subsection("cell")
        cell.set_param("a", xyz.cell[0])
        cell.set_param("b", xyz.cell[1])
        cell.set_param('c', xyz.cell[2])
        
        coord = subsys.add_subsection("coord")
        matrix_str = []
        for atom in xyz.atoms:
            matrix_str.append([
                atom.name,
                str(atom.x),
                str(atom.y),
                str(atom.z)
            ])
        coord.section_var.set("", matrix_str)


    def set_pot_basis(self, kind_basis={}, kind_pot={}, basis_set_file="BASIS_MOLOPT", potential_file="GTH_POTENTIALS"):
        """
        :param kind_basis: {"Li": DZVP-MOLOPT-SR-GTH, ...}
            can be empty dict as all the elemnet has default setting
        :param kind_pot: {"Li": GTH-PBE, ...}
            can be empty dict as all the elemnet has default setting
        :param basis_set_file: PATH_TO_THE_BASIS_SET_FILE_NAME
        :param potential_file: PATH_TO_THE_POTENTIAL_FILE_NAME
        """
        self.set_params({
            "force_eval-dft-BASIS_SET_FILE_NAME": basis_set_file if basis_set_file is not None else "BASIS_MOLOPT", 
            "force_eval-dft-potential_file_name": potential_file if potential_file is not None else "GTH_POTENTIALS",
            })

        for item in kind_basis:
            self.basis_set[item] =  kind_basis[item]
        for item in kind_pot:
            self.potential[item] = kind_pot[item]


    def set_vdw(self, usevdw=False):
        if usevdw == True:
            self.sections["force_eval"].subsections["dft"].subsections["xc"].subsections["vdw_potential"].status = True
        else:
            self.sections["force_eval"].subsections["dft"].subsections["xc"].subsections["vdw_potential"].status = False


    def set_printout(self, option=[]):
        """
        Note:
            responsible for the parseing of the printout_option
        :param option:
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
        """
        self.set_section_status({
            "force_eval-dft-print": True,
            "force_eval-properties": True,
        })


        if 1 in option:
            #self.force_eval.dft.printout.pdos.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["pdos"].status = True
            self.set_section_status({
                "force_eval-dft-print-pdos": True
            })
        if 2 in option:
            #self.force_eval.dft.printout.band_structure.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["band_structure"].status = True
            self.set_section_status({
                "force_eval-dft-print-band_structure": True
            })
            # simply set status to True will not work !!!!
            # you have to also set kpath throught the following commented function
            #self.force_eval.dft.printout.band_structure.set_band(kpath=kpath)
        if 3 in option:
            #self.force_eval.dft.printout.e_density_cube.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["e_density_cube"].status = True
            self.set_section_status({
                "force_eval-dft-print-e_density_cube": True
            })            
        if 4 in option:
            #self.force_eval.dft.printout.elf_cube.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["elf_cube"].status = True
            self.set_section_status({
                "force_eval-dft-print-elf_cube": True
            })            
        if 5 in option:
            #self.force_eval.dft.printout.mo.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["mo"].status = True
            self.set_section_status({
                "force_eval-dft-print-mo": True
            })        
        if 6 in option:
            #self.force_eval.dft.printout.mo_cubes.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["mo_cubes"].status = True
            self.set_section_status({
                "force_eval-dft-print-mo_cubes": True
            })            
        if 7 in option:
            #self.force_eval.dft.printout.mulliken.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["mulliken"].status = True
            self.set_section_status({
                "force_eval-dft-print-mulliken": True
            })            
        if 8 in option:
            #self.force_eval.dft.printout.stm.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["stm"].status = True
            self.set_section_status({
                "force_eval-dft-print-stm": True
            })            
        if 9 in option:
            #self.force_eval.dft.printout.tot_density_cube.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["tot_density_cube"].status = True
            self.set_section_status({
                "force_eval-dft-print-tot_density_cube": True
            })            
        if 10 in option:
            #self.force_eval.dft.printout.v_hartree_cube.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["v_hartree_cube"].status = True
            self.set_section_status({
                "force_eval-dft-print-v_hartree_cube": True
            })            
        if 11 in option:
            #self.force_eval.dft.printout.v_xc_cube.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["v_xc_cube"].status = True
            self.set_section_status({
                "force_eval-dft-print-v_xc_cube": True
            })            
        if 12 in option:
            #self.force_eval.dft.printout.xray_diffraction_spectrum.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["xray_diffraction_spectrum"].status = True
            self.set_section_status({
                "force_eval-dft-print-xray_diffraction_spectrum": True
            })            
        if 13 in option:
            #self.force_eval.properties.resp.status = True
            #self.sections["force_eval"].subsections["properties"].subsections["resp"].status = True
            self.set_section_status({
                "force_eval-properties-resp": True
            })            
        if 14 in option:
            #self.force_eval.dft.printout.moments.status = True
            #self.sections["force_eval"].subsections["dft"].subsections["print"].subsections["moments"].status = True
            self.set_section_status({
                "force_eval-dft-print-moments": True
            })            

    def set_run(self, mpi="", server="pbs", jobname="cp2k", nodes=1, ppn=32, queue=None):
        """ used to set  the parameters controlling the running of the task
        :param mpi: you can specify the mpi command here, it only has effect on native running

        """
        self.run_params["server"] = server
        self.run_params["mpi"] = mpi
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn
        self.run_params["queue"] = queue
        #self.run_params["inpname"] = inpname
        #self.run_params["output"] = output

    def set_llhpc(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr


    def gen_llhpc(self, inpname, output, directory, cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        better pass in $PMF_CP2K
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".slurm"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("yhrun %s -in %s | tee %s\n" % (cmd, inpname, output))


    def gen_yh(self, inpname, output, directory, cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -in %s | tee %s\n" % (cmd, inpname, output))


    def gen_pbs(self, inpname, output, directory, cmd="cp2k.popt", jobname="cp2k", nodes=1, ppn=32, queue=None):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            if queue != None:
                fout.write("#PBS -q %s\n" % queue)
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s -in %s > %s\n" % (cmd, inpname, output))

    def set_cdcloud(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr


    def gen_cdcloud(self, inpname, output, directory, cmd="cp2k.psmp"):
        """
        generating yhbatch job script for calculation
        better pass in $PMF_CP2K
        """
        with open(os.path.join(directory, inpname.split(".inp")[0]+".slurm"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
            fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
            fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
            fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
            fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
            fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
            fout.write("#\n")
            fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
            fout.write("export FORT_BUFFERED=1\n")
            fout.write("srun --mpi=pmix_v3 %s -in %s | tee %s\n" % (cmd, inpname, output))
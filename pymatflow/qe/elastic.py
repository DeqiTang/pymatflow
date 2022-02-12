"""
calculate elastic properties of the system
"""
import re
import numpy as np

"""
calculate elastic constants of the system.
Reference:
    https://icme.hpc.msstate.edu/mediawiki/index.php/Calculation_of_elastic_constants#Elastic_constant_calculation_script
    http://muchong.com/html/201005/2098932.html
    https://wenku.baidu.com/view/0c6fecc458f5f61fb736668d.html
    Calculations of single-crystal elastic constants made simple”（Computer Physics Communications Volume 181, Issue 3, March 2010, Pages 671-675 ）

Note:
    there are two different ways to calculate elastic constant: stress-strain relationship, strain energy-strain relationship.

    strain_tensor = [
    [e1, 0.5*e6, 0.5*e5],
    [0.5*e6, e2, 0.5*e4],
    [0.5*e5, 0.5*e4, e3]
    ]

"""


def epsilon_stress_strain(xyz):
    """
    using stress-strain relationship to calculate elastica constant

    :param xyz:
        an instance of base_xyz()
    """
    delta = 0.05 # in unit of anstrom
    e = np.full(6, delta)
    strain_tensor = np.zeros([3, 3])
    strain_tensor[0, 0] = e[0]
    strain_tensor[1, 1] = e[1]
    strain_tensor[2, 2] = e[2]
    strain_tensor[0, 1] = 0.5 * e[5]
    strain_tensor[0, 2] = 0.5 * e[4]
    strain_tensor[1, 0] = 0.5 * e[5]
    strain_tensor[1, 2] = 0.5 * e[3]
    strain_tensor[2, 0] = 0.5 * e[4]
    strain_tensor[2, 1] =0.5 * e[3]

    #cell = np.array([xyz.cell[0:3], xyz.cell[3:6], xyz.cell[6:9]])
    cell = np.array(xyz.cell)
    new_cell = np.dot(cell, (np.eye(3) + strain_tensor))

    return new_cell


class ElasticRun:
    """
    """
    def __init__(self, xyz_f):
        self.file = xyz_f


    def elastic(self, directory="tmp-qe-elastic",
            mpi="", runopt="gen", control={}, system={}, electrons={}, ions={},
            kpoints_option="automatic", kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        :param directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)

        #os.system("cp *.UPF %s/" % directory)
        #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

        # do not copy too many files at the same time or it will be slow
        # so we do not copy all UPF files in the directory but just copy
        # those used in the calculation.
        shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
        #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
        all_file = os.listdir()
        for element in self.arts.xyz.specie_labels:
            for item in all_file:
                #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                    shutil.copyfile(item, os.path.join(directory, item))
                    break
        #


        os.chdir(directory)

        relax = opt_run(self.file)
        # need to get the stress tensor
        relax.control.params["tprnfor"] = True
        relax.control.params["tstress"] = True


        relax.relax(directory="tmp-qe-relax", control=control, system=system, electrons=electrons, ions=ions, kpoints_option=kpoints_option, kpoints_mp=kpoints_mp)

import os
import copy
import numpy as np
import matplotlib.pyplot as plt

from xml.etree.ElementTree import parse

class post_pdos:
    def __init__(self):
        pass

    def get_vasprun(self, vasprun="vasprun.xml"):
        self.vasprun = parse(vasprun)
        self.efermi = float(self.vasprun.getroot().find("calculation").find("dos").find("i").text)
        self.get_dos_partial()



    def get_dos_partial(self):
        if len(self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set").find("set")) == 1:
            self.magnetic_status = "spin-unpolarized"
        elif len(self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set").find("set")) == 2:
            self.magnetic_status = "spin-polarized"

        #
        if self.magnetic_status == "spin-unpolarized":
            """
            self.data:
                [
                    {
                        "ion": 'element label',
                        "spin_1": {"energy": [], "s": [], "py": [], .......}
                    },
                    ...
                ]
            a list contains partial dos data for each ion
            """
            self.data = []
            nion = len(self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set"))
            # field is in format of {"energy": [], "s": [], "py": [], .....}
            # used to initialize data of every ion
            field = {}
            for item in self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").findall("field"):
                field[item.text.split()[0]] = []
            #
            for ion in self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set"):
                iondata = {}
                iondata["spin_1"] = copy.deepcopy(field)
                for line in ion[0]:
                    i = 0
                    for key in iondata["spin_1"].keys():
                        iondata["spin_1"][key].append(float(line.text.split()[i]))
                        i = i + 1
                self.data.append(iondata)

        elif self.magnetic_status == "spin-polarized":
            """
            self.data:
                [
                    {
                        "ion": 'element label',
                        "spin_1": {"energy": [], "s": [], "py": [], .......},
                        "spin_2": {"energy": [], "s": [], "py": [], .......}
                    },
                    ...
                ]
            a list contains partial dos data for each ion
            """
            self.data = []
            nion = len(self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set"))
            # field is in format of {"energy": [], "s": [], "py": [], .....}
            # used to initialize data of every ion
            field = {}
            for item in self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").findall("field"):
                field[item.text.split()[0]] = []
            #
            for ion in self.vasprun.getroot().find("calculation").find("dos").find("partial").find("array").find("set"):
                iondata = {}
                iondata["spin_1"] = copy.deepcopy(field)
                iondata["spin_2"] = copy.deepcopy(field)
                for line in ion[0]:
                    i = 0
                    for key in iondata["spin_1"].keys():
                        iondata["spin_1"][key].append(float(line.text.split()[i]))
                        i = i + 1
                for line in ion[1]:
                    i = 0
                    for key in iondata["spin_2"].keys():
                        iondata["spin_2"][key].append(float(line.text.split()[i]))
                        i = i + 1
                self.data.append(iondata)
        # get element label for each ion in self.data
        # self.ion_list: ['H', 'H', 'He', 'Li', 'Li', 'Li', 'C', "N", 'N', ......]
        ion_list = []
        for ion in self.vasprun.getroot().find("atominfo").find("array").find("set"):
            ion_list.append(ion[0].text.split()[0])
        for i in range(len(ion_list)):
            self.data[i]["ion"] = ion_list[i]

        #
    def matplotlib_proj_elem_l_m(self, plotrange=[0, 1]):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])
        if self.magnetic_status == "spin-unpolarized":
            """
            data:{
                'element-label': {
                    "spin_1": {"energy": [], "s": [], "py": [], .......},
                },
                .....
            }
            """
            data = {}
            for ion in self.data:
                if ion["ion"] not in data:
                    data[ion["ion"]] = ion
                else:
                    for key in data[ion["ion"]]["spin_1"].keys():
                        if key != 'energy':
                            data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                        else:
                            continue
            # make the plot
            for element in data:
                for lm in data[element]["spin_1"]:
                    if lm != "energy":
                        plt.plot(data[element]["spin_1"]["energy"][begin:end], data[element]["spin_1"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm.png")
            plt.show()
            plt.close()

    def export(self, directory="tmp-vasp-static", plotrange=[0, 1], option="matplotlib"):
        """
        :parama option:
            gnuplot or matplotlib
        :param plotrange:
            a list of two values(between 0 and 1) defining the percentage
            of bands to plot.
            plotrange[0]: left boundary of the nth band to plot
            plotrange[1]: right boundary of the nth band to plot
            default is plotrange[0] = 0, plotrange[1], in which case
            all the band available will be plot.
            Be aware that the range if not for energy but for band number
        """
        os.system("mkdir -p %s/post-processing" % directory)
        os.chdir(os.path.join(directory, "post-processing"))
        #(option=option,  plotrange=plotrange)
        self.matplotlib_proj_elem_l_m(plotrange=plotrange)
        os.chdir("../../")

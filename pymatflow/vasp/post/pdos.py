import os
import copy
import numpy as np
import matplotlib.pyplot as plt

from xml.etree.ElementTree import parse

class post_pdos:
    def __init__(self):
        pass


    def get_efermi(self, vasprun="vasprun.xml"):
        """
        we set efermi in an individual function because we can choose to get efermi from nscf run
        or scf run in this way.
        if you want to get efermi from the scf run specify the vasprun.xml for the scf
        if you want to get efermi from the nscf run specify the vasprun.xml for the nscf
        """
        vasprun_xml = parse(vasprun)
        self.efermi = float(vasprun_xml.getroot().find("calculation").find("dos").find("i").text)

    def get_vasprun(self, vasprun="vasprun.xml"):
        self.vasprun = parse(vasprun)
        self.get_dos_partial()


    def get_dos_partial(self):
        # first check the magnetic_status
        for i in range(len(self.vasprun.getroot().find("parameters"))):
            if self.vasprun.getroot().find("parameters").getchildren()[i].attrib["name"] == "electronic":
                for item in self.vasprun.getroot().find("parameters")[i]:
                    if item.attrib["name"] == "electronic spin":
                        for spin_item in item:
                            if spin_item.attrib["name"] == "ISPIN":
                                ispin = int(spin_item.text.split()[0])
                            elif spin_item.attrib["name"] == "LNONCOLLINEAR":
                                lnoncollinear = spin_item.text.split()[0]
                            elif spin_item.attrib["name"] == "LSORBIT":
                                lsorbit = spin_item.text.split()[0]
        if lsorbit == "T":
            self.magnetic_status = "soc-ispin-%d" % ispin # soc-ispin-1 or soc-ispin-2
            # actually soc-ispin-2 never exists !!!
            # I found when soc and ISPIN=2 are set in the INCAR at the same time,
            # you can find in <incar> section of vasprun.xml ISPIN=2, but it will be turned to 1
            # in <parameters> -> <separator name="electronic"> - > <separator name="electronic spin"> 
            # So VASP will automatically set ISPIN to 1 even when you set ISPIN to 2 
            # if you are considering soc
            # I mean you can set soc and ISPIN=2 at the same time in INCAR, but 
            # VASP will turn ISPIN to 1, and this post script will read it to be ISPIN = 1
            # so never will therere be soc-ispin-2 even when you set it in INCAR.
            # in vasprun.xml, we read the acutally used ISPIN from <parameters>... rather than
            # the input from <incar>.
        else:
            self.magnetic_status = "non-soc-ispin-%d" % ispin # non-soc-ispin-1 or non-soc-ispin-2

        #
        if self.magnetic_status == "non-soc-ispin-1":
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
                        #print(line.text.split()[i])
                        iondata["spin_1"][key].append(float(line.text.split()[i]))
                        i = i + 1
                self.data.append(iondata)

        elif self.magnetic_status == "non-soc-ispin-2":
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
        elif self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            """
            self.data:
                [
                    {
                        "ion": 'element label',
                        "spin_1": {"energy": [], "s": [], "py": [], .......},
                        "spin_2": {"energy": [], "s": [], "py": [], .......},
                        "spin_3": {"energy": [], "s": [], "py": [], .......},
                        "spin_4": {"energy": [], "s": [], "py": [], .......}
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
                iondata["spin_3"] = copy.deepcopy(field)
                iondata["spin_4"] = copy.deepcopy(field)
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
                for line in ion[2]:
                    i = 0
                    for key in iondata["spin_3"].keys():
                        iondata["spin_3"][key].append(float(line.text.split()[i]))
                        i = i + 1
                for line in ion[3]:
                    i = 0
                    for key in iondata["spin_4"].keys():
                        iondata["spin_4"][key].append(float(line.text.split()[i]))
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

    def plot_proj_elem_l_m(self, plotrange=[0, 1], engine="matplotlib"):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])

        # construct data
        """
        in case of non-soc-ispin-1
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "py": [], .......},
            },
            .....
        }
        in case of non-soc-ispin-2
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "py": [], .......},
                "spin_2": {"energy": [], "s": [], "py": [], .......},                    
            },
            .....
        }
        in case of soc-ispin-1 or soc-ispin-2
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "py": [], .......},
                "spin_2": {"energy": [], "s": [], "py": [], .......},                    
                "spin_3": {"energy": [], "s": [], "py": [], .......},    
                "spin_4": {"energy": [], "s": [], "py": [], .......},                                            
            },
            .....
        }
        """        
        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            spin_n = 4



        data = {}
        for ion in self.data:
            if ion["ion"] not in data:
                data[ion["ion"]] = ion
            else:
                for i in range(spin_n):
                    for key in data[ion["ion"]]["spin_%d" % (i+1)].keys():
                        if key != 'energy':
                            for j in range(len(ion["spin_%d" % (i+1)][key])):
                                data[ion["ion"]]["spin_%d" % (i+1)][key][j] += ion["spin_%d" % (i+1)][key][j] 
                        else:
                            continue

        # export data in gnuplot format
        for element in data:
            for i in range(spin_n):
                with open("pdos-proj-elem(%s)-lm-%s-spin-%d.data" % (element, self.magnetic_status, i+1), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_%d" % (i+1)]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                        for item in data[element]["spin_%d" % (i+1)]:
                            fout.write("%f " % data[element]["spin_%d" % (i+1)][item][j])
                        fout.write("\n")
        
        if engine == "matplotlib":
            # make the plot
            for i in range(spin_n):
                for element in data:
                    for lm in data[element]["spin_%d" % (i+1)]:
                        if lm != "energy":
                            plt.plot([en - self.efermi for en in data[element]["spin_%d" % (i+1)]["energy"][begin:end]], data[element]["spin_%d" % (i+1)][lm][begin:end], label=element+"(%s)" % lm)
                plt.title("Partial Density of States (Spin %d)" % (i+1))
                plt.xlabel("Energy(ev)")
                plt.ylabel("PDOS")
                plt.legend()
                plt.savefig("pdos-proj-elem-lm-%s-spin-%d.png" % (self.magnetic_status, i+1))
                #plt.show()
                plt.close()
        elif engine == "gnuplot":
            # make the plot
            for i in range(spin_n):
                with open("pdos-proj-elem-lm-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1), 'w') as fout:
                    fout.write("set terminal gif\n")
                    fout.write("set output 'pdos-proj-elem-lm-%s-spin-%d.gif'\n" % (self.magnetic_status, i+1))
                    #fout.write("unset key\n")
                    fout.write("set parametric\n")
                    fout.write("set title 'Partial Density of States'\n")
                    fout.write("set xlabel 'Energy(eV)'\n")
                    fout.write("set ylabel 'PDOS'\n")
                    fout.write("set grid xtics ytics\n")
                    fout.write("set autoscale\n")
                    fout.write("plot ")
                    for element in data:
                        col = 2
                        for lm in data[element]["spin_%d" % (i+1)]:
                            if lm == "energy":
                                continue
                            fout.write("'pdos-proj-elem(%s)-lm-%s-spin-%d.data' using ($1-%f):%d title '%s' w l, \\\n" % (element, self.magnetic_status, i+1, self.efermi, col, element+"("+lm+")"))
                            col += 1
                os.system("gnuplot pdos-proj-elem-lm-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1))


    def plot_proj_elem_l(self, plotrange=[0, 1], engine="matplotlib"):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])

        # construct data
        """
        in case of non-soc-ispin-1
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "p": [], 'd': [] .......},
            },
            .....
        }
        in case of non-soc-ispin-2
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "p": [], 'd': [], .......},
                "spin_2": {"energy": [], "s": [], "p": [], 'd': [].......},                    
            },
            .....
        }
        in case of soc-ispin-1 or soc-ispin-2
        data:{
            'element-label': {
                "ion": "element-label",
                "spin_1": {"energy": [], "s": [], "p": [], "d": [], .......},
                "spin_2": {"energy": [], "s": [], "p": [], "d": [], .......},                    
                "spin_3": {"energy": [], "s": [], "p": [], "d": [], .......},    
                "spin_4": {"energy": [], "s": [], "p": [], "d": [], .......},                                            
            },
            .....
        }
        """        
        if self.magnetic_status == "non-soc-ispin-1":
            spin_n = 1
        elif self.magnetic_status == "non-soc-ispin-2":
            spin_n = 2
        elif self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            spin_n = 4
        data = {}
        for ion in self.data:
            if ion["ion"] not in data:
                data[ion["ion"]] = ion
            else:
                for i in range(spin_n):
                    for key in data[ion["ion"]]["spin_%d" % (i+1)].keys():
                        if key != 'energy':
                            for j in range(len(ion["spin_%d" % (i+1)][key])):
                                data[ion["ion"]]["spin_%d" % (i+1)][key][j] += ion["spin_%d" % (i+1)][key][j] 
                        else:
                            continue
        # now l are decomposed to m, we should merge m for every l
        for element in data:
            n_orb = len(data[element]["spin_1"]) - 1
            break
        if n_orb == 1:
            # only s exits, no need to adjust
            pass
        elif n_orb == 4:
            # s and p exist, need to  merge px py pz to p
            for element in data:
                for i in range(spin_n):
                    data[element]["spin_%d" % (i+1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                    for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                        data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"][j]
                    data[element]["spin_%d" % (i+1)].pop("px")
                    data[element]["spin_%d" % (i+1)].pop("py")
                    data[element]["spin_%d" % (i+1)].pop("pz")                                
        elif n_orb == 9:
            # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
            for element in data:
                for i in range(spin_n):
                    data[element]["spin_%d" % (i+1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                    data[element]["spin_%d" % (i+1)]["d"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["dxy"])                    
                    for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                        data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"][j]
                        data[element]["spin_%d" % (i+1)]["d"][j] += data[element]["spin_%d" % (i+1)]["dyz"][j] + data[element]["spin_%d" % (i+1)]["dz2"][j] + data[element]["spin_%d" % (i+1)]["dxz"][j] + data[element]["spin_%d" % (i+1)]["x2-y2"][j]
                    data[element]["spin_%d" % (i+1)].pop("px")
                    data[element]["spin_%d" % (i+1)].pop("py")
                    data[element]["spin_%d" % (i+1)].pop("pz")                    
                    data[element]["spin_%d" % (i+1)].pop("dxy")     
                    data[element]["spin_%d" % (i+1)].pop("dyz")     
                    data[element]["spin_%d" % (i+1)].pop("dz2")     
                    data[element]["spin_%d" % (i+1)].pop("dxz")     
                    data[element]["spin_%d" % (i+1)].pop("x2-y2")                                                                                                         
        elif n_orb == 16:
            # s, p and d and f exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d and fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3 to f
            for element in data:
                for i in range(spin_n):
                    data[element]["spin_%d" % (i+1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                    data[element]["spin_%d" % (i+1)]["d"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["dxy"])
                    data[element]["spin_%d" % (i+1)]["f"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["fy3x2"])
                    for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                        data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"][j]
                        data[element]["spin_%d" % (i+1)]["d"][j] += data[element]["spin_%d" % (i+1)]["dyz"][j] + data[element]["spin_%d" % (i+1)]["dz2"][j] + data[element]["spin_%d" % (i+1)]["dxz"][j] + data[element]["spin_%d" % (i+1)]["x2-y2"][j]
                        data[element]["spin_%d" % (i+1)]["f"][j] += data[element]["spin_%d" % (i+1)]["fxyz"][j] + data[element]["spin_%d" % (i+1)]["fyz2"][j] + data[element]["spin_%d" % (i+1)]["fz3"][j] + data[element]["spin_%d" % (i+1)]["fxz2"][j] + data[element]["spin_%d" % (i+1)]["fzx2"][j] + data[element]["spin_%d" % (i+1)]["fx3"][j]
                    data[element]["spin_%d" % (i+1)].pop("px")
                    data[element]["spin_%d" % (i+1)].pop("py")
                    data[element]["spin_%d" % (i+1)].pop("pz")                    
                    data[element]["spin_%d" % (i+1)].pop("dxy")     
                    data[element]["spin_%d" % (i+1)].pop("dyz")     
                    data[element]["spin_%d" % (i+1)].pop("dz2")     
                    data[element]["spin_%d" % (i+1)].pop("dxz")     
                    data[element]["spin_%d" % (i+1)].pop("x2-y2")     
                    data[element]["spin_%d" % (i+1)].pop("fy3x2")
                    data[element]["spin_%d" % (i+1)].pop("fxyz")                    
                    data[element]["spin_%d" % (i+1)].pop("fyz2")     
                    data[element]["spin_%d" % (i+1)].pop("fz3")     
                    data[element]["spin_%d" % (i+1)].pop("fxz2")     
                    data[element]["spin_%d" % (i+1)].pop("fzx2")     
                    data[element]["spin_%d" % (i+1)].pop("fx3")                       

        # export data in gnuplot format
        for i in range(spin_n):
            for element in data:
                with open("pdos-proj-elem(%s)-l-%s-spin-%d.data" % (element, self.magnetic_status, i+1), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_%d" % (i+1)]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                        for item in data[element]["spin_%d" % (i+1)]:
                            fout.write("%f " % data[element]["spin_%d" % (i+1)][item][j])
                        fout.write("\n")

        if engine == "matplotlib":
            # make the plot
            for i in range(spin_n):
                for element in data:
                    for l in data[element]["spin_%d" % (i+1)]:
                        if l != "energy":
                            plt.plot([en - self.efermi for en in data[element]["spin_%d" % (i+1)]["energy"][begin:end]], data[element]["spin_%d" % (i+1)][l][begin:end], label=element+"(%s)" % l)
                plt.title("Partial Density of States")
                plt.xlabel("Energy(ev)")
                plt.ylabel("PDOS")
                plt.legend()
                plt.savefig("pdos-proj-elem-l-%s-spin-%d.png" % (self.magnetic_status, i+1))
                #plt.show()
                plt.close()
        elif engine == "gnuplot":
            # make the plot
            for i in range(spin_n):
                with open("pdos-proj-elem-l-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1), 'w') as fout:
                    fout.write("set terminal gif\n")
                    fout.write("set output 'pdos-proj-elem-l-%s-spin-%d.gif'\n" % (self.magnetic_status, i+1))
                    #fout.write("unset key\n")
                    fout.write("set parametric\n")
                    fout.write("set title 'Partial Density of States'\n")
                    fout.write("set xlabel 'Energy(eV)'\n")
                    fout.write("set ylabel 'PDOS'\n")
                    fout.write("set grid xtics ytics\n")
                    fout.write("set autoscale\n")
                    fout.write("plot ")
                    for element in data:
                        col = 2
                        for l in data[element]["spin_%d" % (i+1)]:
                            if l == "energy":
                                continue
                            fout.write("'pdos-proj-elem(%s)-l-%s-spin-%d.data' using ($1-%f):%d title '%s' w l, \\\n" % (element, self.magnetic_status, i+1, self.efermi, col, element+"("+l+")"))
                            col += 1
                os.system("gnuplot pdos-proj-elem-l-%s-spin-%d.gnuplot" % (self.magnetic_status, i+1))                        


    def export(self, directory="tmp-vasp-static", plotrange=[0, 1], engine="matplotlib"):
        """
        :param engine:
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
        self.plot_proj_elem_l_m(plotrange=plotrange, engine=engine)
        self.plot_proj_elem_l(plotrange=plotrange, engine=engine)
        os.chdir("../../")

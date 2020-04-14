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
    def matplotlib_proj_elem_l_m(self, plotrange=[0, 1]):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])
        if self.magnetic_status == "non-soc-ispin-1":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                        else:
                            continue
            # make the plot
            for element in data:
                for lm in data[element]["spin_1"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_1"]["energy"][begin:end]], data[element]["spin_1"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s.png" % self.magnetic_status)
            #plt.show()
            plt.close()

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")

        if self.magnetic_status == "non-soc-ispin-2":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "py": [], .......},
                    "spin_2": {"energy": [], "s": [], "py": [], .......},                    
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]                                                             
                        else:
                            continue
            # make the plot
            for element in data:
                for lm in data[element]["spin_1"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_1"]["energy"][begin:end]], data[element]["spin_1"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 1)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%sspin-1.png" % self.magnetic_status)
            #plt.show()
            plt.close()
            for element in data:
                for lm in data[element]["spin_2"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_2"]["energy"][begin:end]], data[element]["spin_2"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 2)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s-spin-2.png" % self.magnetic_status)
            #plt.show()
            plt.close()

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")

        if self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            """
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
            data = {}
            for ion in self.data:
                if ion["ion"] not in data:
                    data[ion["ion"]] = ion
                else:
                    for key in data[ion["ion"]]["spin_1"].keys():
                        if key != 'energy':
                            data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]                            
                            data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]                            
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]     
                            #data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]
                            for i in range(len(ion["spin_3"][key])):
                                data[ion["ion"]]["spin_3"][key][i] += ion["spin_3"][key][i] 
                            #data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]
                            for i in range(len(ion["spin_4"][key])):
                                data[ion["ion"]]["spin_4"][key][i] += ion["spin_4"][key][i]                                                                 
                        else:
                            continue
            # make the plot
            for element in data:
                for lm in data[element]["spin_1"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_1"]["energy"][begin:end]], data[element]["spin_1"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 1)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s-spin-1.png" % self.magnetic_status)
            #plt.show()
            plt.close()
            for element in data:
                for lm in data[element]["spin_2"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_2"]["energy"][begin:end]], data[element]["spin_2"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 2)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s-spin-2.png" % self.magnetic_status)
            #plt.show()
            plt.close()
            for element in data:
                for lm in data[element]["spin_3"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_3"]["energy"][begin:end]], data[element]["spin_3"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 3)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s-spin-3.png" % self.magnetic_status)
            #plt.show()
            plt.close()
            for element in data:
                for lm in data[element]["spin_4"]:
                    if lm != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_4"]["energy"][begin:end]], data[element]["spin_4"][lm][begin:end], label=element+"(%s)" % lm)
            plt.title("Partial Density of States (Spin 4)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-lm-%s-spin-4.png" % self.magnetic_status)
            #plt.show()
            plt.close()            

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-3.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_3"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_3"]["energy"])):
                        for item in data[element]["spin_3"]:
                            fout.write("%f " % data[element]["spin_3"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-4.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_4"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_4"]["energy"])):
                        for item in data[element]["spin_4"]:
                            fout.write("%f " % data[element]["spin_4"][item][i])
                        fout.write("\n")


    def gnuplot_proj_elem_l_m(self, plotrange=[0, 1]):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])
        if self.magnetic_status == "non-soc-ispin-1":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                        else:
                            continue

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
            # make the plot
            for i in range(1):
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

        if self.magnetic_status == "non-soc-ispin-2":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "py": [], .......},
                    "spin_2": {"energy": [], "s": [], "py": [], .......},                    
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]                                                             
                        else:
                            continue

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")
            # make the plot
            for i in range(2):
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


        if self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            """
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
            data = {}
            for ion in self.data:
                if ion["ion"] not in data:
                    data[ion["ion"]] = ion
                else:
                    for key in data[ion["ion"]]["spin_1"].keys():
                        if key != 'energy':
                            data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]                            
                            data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]                            
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]     
                            #data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]
                            for i in range(len(ion["spin_3"][key])):
                                data[ion["ion"]]["spin_3"][key][i] += ion["spin_3"][key][i] 
                            #data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]
                            for i in range(len(ion["spin_4"][key])):
                                data[ion["ion"]]["spin_4"][key][i] += ion["spin_4"][key][i]                                                                 
                        else:
                            continue

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-lm-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)                    
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-3.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)                    
                    fout.write("#")
                    for item in data[element]["spin_3"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_3"]["energy"])):
                        for item in data[element]["spin_3"]:
                            fout.write("%f " % data[element]["spin_3"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-lm-%s-spin-4.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)                    
                    fout.write("#")
                    for item in data[element]["spin_4"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_4"]["energy"])):
                        for item in data[element]["spin_4"]:
                            fout.write("%f " % data[element]["spin_4"][item][i])
                        fout.write("\n")
            
            # make the plot
            for i in range(4):
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


    def matplotlib_proj_elem_l(self, plotrange=[0, 1]):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])
        if self.magnetic_status == "non-soc-ispin-1":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "p": [], 'd': [] .......},
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
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
                    for i in range(1):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(1):
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
                    for i in range(1):
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

            # make the plot
            for element in data:
                for l in data[element]["spin_1"]:
                    if l != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_1"]["energy"][begin:end]], data[element]["spin_1"][l][begin:end], label=element+"(%s)" % l)
            plt.title("Partial Density of States")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-l-%s.png" % self.magnetic_status)
            #plt.show()
            plt.close()

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-l-%s.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")

        if self.magnetic_status == "non-soc-ispin-2":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "py": [], .......},
                    "spin_2": {"energy": [], "s": [], "py": [], .......},                    
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]                                                             
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
                    for i in range(2):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(2):
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
                    for i in range(2):
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

            # make the plot
            for element in data:
                for l in data[element]["spin_1"]:
                    if l != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_1"]["energy"][begin:end]], data[element]["spin_1"][l][begin:end], label=element+"(%s)" % l)
            plt.title("Partial Density of States (Spin 1)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-l-%sspin-1.png" % self.magnetic_status)
            #plt.show()
            plt.close()
            for element in data:
                for l in data[element]["spin_2"]:
                    if l != "energy":
                        plt.plot([en - self.efermi for en in data[element]["spin_2"]["energy"][begin:end]], data[element]["spin_2"][l][begin:end], label=element+"(%s)" % l)
            plt.title("Partial Density of States (Spin 2)")
            plt.xlabel("Energy(ev)")
            plt.ylabel("PDOS")
            plt.legend()
            plt.savefig("pdos-proj-elem-l-%s-spin-2.png" % self.magnetic_status)
            #plt.show()
            plt.close()

            # export data in gnuplot format
            for element in data:
                with open("pdos-proj-elem(%s)-l-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-l-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")

        if self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            """
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
            data = {}
            for ion in self.data:
                if ion["ion"] not in data:
                    data[ion["ion"]] = ion
                else:
                    for key in data[ion["ion"]]["spin_1"].keys():
                        if key != 'energy':
                            data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]                            
                            data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]                            
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]     
                            #data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]
                            for i in range(len(ion["spin_3"][key])):
                                data[ion["ion"]]["spin_3"][key][i] += ion["spin_3"][key][i] 
                            #data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]
                            for i in range(len(ion["spin_4"][key])):
                                data[ion["ion"]]["spin_4"][key][i] += ion["spin_4"][key][i]                                                                 
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
                    for i in range(4):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(4):
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
                    for i in range(4):
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
            
            # make the plot
            for i in range(4):
                for element in data:
                    for l in data[element]["spin_%d" % (i+1)]:
                        if l != "energy":
                            plt.plot([en - self.efermi for en in data[element]["spin_%d" % (i+1)]["energy"][begin:end]], data[element]["spin_%d" % (i+1)][l][begin:end], label=element+"(%s)" % l)
                plt.title("Partial Density of States (Spin %d)" % (i+1))
                plt.xlabel("Energy(ev)")
                plt.ylabel("PDOS")
                plt.legend()
                plt.savefig("pdos-proj-elem-l-%s-spin-%d.png" % (self.magnetic_status, i+1))
                #plt.show()
                plt.close()
    
            # export data in gnuplot format
            for i in range(4):
                for element in data:
                    with open("pdos-proj-elem(%s)-l-%s-spin-%d.data" % (element, self.magnetic_status, i+1), 'w') as fout:
                        fout.write("# Efermi: %f\n" % self.efermi)
                        fout.write("#")
                        for item in data[element]["spin_%d" % (i+1)]:
                            fout.write(" %s" % item)
                        fout.write("\n")
                        for i in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            for item in data[element]["spin_%d" % (i+1)]:
                                fout.write("%f " % data[element]["spin_%d" % (i+1)][item][i])
                            fout.write("\n")

    def gnuplot_proj_elem_l(self, plotrange=[0, 1]):
        begin = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[0])
        end = int(len(self.data[0]["spin_1"]["energy"]) * plotrange[1])
        if self.magnetic_status == "non-soc-ispin-1":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "p": [], 'd': [] .......},
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
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
                    for i in range(1):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(1):
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
                    for i in range(1):
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
            for element in data:
                with open("pdos-proj-elem(%s)-l-%s.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
            
            # make the plot
            for i in range(1):
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

        if self.magnetic_status == "non-soc-ispin-2":
            """
            data:{
                'element-label': {
                    "ion": "element-label",
                    "spin_1": {"energy": [], "s": [], "py": [], .......},
                    "spin_2": {"energy": [], "s": [], "py": [], .......},                    
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
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]                                                             
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
                    for i in range(2):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(2):
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
                    for i in range(2):
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
            for element in data:
                with open("pdos-proj-elem(%s)-l-%s-spin-1.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_1"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_1"]["energy"])):
                        for item in data[element]["spin_1"]:
                            fout.write("%f " % data[element]["spin_1"][item][i])
                        fout.write("\n")
                with open("pdos-proj-elem(%s)-l-%s-spin-2.data" % (element, self.magnetic_status), 'w') as fout:
                    fout.write("# Efermi: %f\n" % self.efermi)
                    fout.write("#")
                    for item in data[element]["spin_2"]:
                        fout.write(" %s" % item)
                    fout.write("\n")
                    for i in range(len(data[element]["spin_2"]["energy"])):
                        for item in data[element]["spin_2"]:
                            fout.write("%f " % data[element]["spin_2"][item][i])
                        fout.write("\n")

            # make the plot
            for i in range(2):
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

        if self.magnetic_status == "soc-ispin-1" or self.magnetic_status == "soc-ispin-2":
            """
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
            data = {}
            for ion in self.data:
                if ion["ion"] not in data:
                    data[ion["ion"]] = ion
                else:
                    for key in data[ion["ion"]]["spin_1"].keys():
                        if key != 'energy':
                            data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]                            
                            data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]                            
                            #data[ion["ion"]]["spin_1"][key] = [data[ion["ion"]]["spin_1"][key][i] + ion["spin_1"][key][i] for i in range(len(ion["spin_1"][key])) ]
                            for i in range(len(ion["spin_1"][key])):
                                data[ion["ion"]]["spin_1"][key][i] += ion["spin_1"][key][i] 
                            #data[ion["ion"]]["spin_2"][key] = [data[ion["ion"]]["spin_2"][key][i] + ion["spin_2"][key][i] for i in range(len(ion["spin_2"][key])) ]
                            for i in range(len(ion["spin_2"][key])):
                                data[ion["ion"]]["spin_2"][key][i] += ion["spin_2"][key][i]     
                            #data[ion["ion"]]["spin_3"][key] = [data[ion["ion"]]["spin_3"][key][i] + ion["spin_3"][key][i] for i in range(len(ion["spin_3"][key])) ]
                            for i in range(len(ion["spin_3"][key])):
                                data[ion["ion"]]["spin_3"][key][i] += ion["spin_3"][key][i] 
                            #data[ion["ion"]]["spin_4"][key] = [data[ion["ion"]]["spin_4"][key][i] + ion["spin_4"][key][i] for i in range(len(ion["spin_4"][key])) ]
                            for i in range(len(ion["spin_4"][key])):
                                data[ion["ion"]]["spin_4"][key][i] += ion["spin_4"][key][i]                                                                 
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
                    for i in range(4):
                        data[element]["spin_%d" % (1)]["p"] = copy.deepcopy(data[element]["spin_%d" % (i+1)]["px"])
                        for j in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            data[element]["spin_%d" % (i+1)]["p"][j] += data[element]["spin_%d" % (i+1)]["py"][j] + data[element]["spin_%d" % (i+1)]["pz"]
                        data[element]["spin_%d" % (i+1)].pop("px")
                        data[element]["spin_%d" % (i+1)].pop("py")
                        data[element]["spin_%d" % (i+1)].pop("pz")                                
            elif n_orb == 9:
                # s, p and d exist, need to merge px py pz to p and dxy dyz dz2 dxz x2-y2 to d
                for element in data:
                    for i in range(4):
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
                    for i in range(4):
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
            for i in range(4):
                for element in data:
                    with open("pdos-proj-elem(%s)-l-%s-spin-%d.data" % (element, self.magnetic_status, i+1), 'w') as fout:
                        fout.write("# Efermi: %f\n" % self.efermi)
                        fout.write("#")
                        for item in data[element]["spin_%d" % (i+1)]:
                            fout.write(" %s" % item)
                        fout.write("\n")
                        for i in range(len(data[element]["spin_%d" % (i+1)]["energy"])):
                            for item in data[element]["spin_%d" % (i+1)]:
                                fout.write("%f " % data[element]["spin_%d" % (i+1)][item][i])
                            fout.write("\n")

            # make the plot
            for i in range(4):
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
        if option == "matplotlib":
            #self.matplotlib_proj_elem_l_m(plotrange=plotrange)
            self.matplotlib_proj_elem_l(plotrange=plotrange)            
        elif option == "gnuplot":
            #self.gnuplot_proj_elem_l_m(plotrange=plotrange)
            self.gnuplot_proj_elem_l(plotrange=plotrange)            
        os.chdir("../../")

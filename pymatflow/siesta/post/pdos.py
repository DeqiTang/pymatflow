#!/usr/bin/env python
# _*_ coding: utf-8

import os
import xml.etree.ElementTree
import numpy as np
import matplotlib.pyplot as plt

class orbital:
    """
    """
    def __init__(self, orbital_xml):
        self.index = orbital_xml.attrib["index"]
        self.atom_index = orbital_xml.attrib["atom_index"]
        self.species = orbital_xml.attrib["species"]
        self.Z = orbital_xml.attrib["Z"]
        self.position = orbital_xml.attrib["position"]
        self.n = orbital_xml.attrib["n"]
        self.l = orbital_xml.attrib["l"]
        self.m = orbital_xml.attrib["m"]
        self.z = orbital_xml.attrib["z"]
        self.P = orbital_xml.attrib["P"]
        self.data = np.array([float(orbital_xml[0].text.split()[i]) for i in range(len(orbital_xml[0].text.split()))])

class pdos:
    """
    """
    def __init__(self):
        self.file = None
        
    def get_info(self, pdos_f):
        self.file = os.path.abspath(pdos_f)
        tree = xml.etree.ElementTree.parse(self.file)
        root = tree.getroot()
        
        self.nspin = int(root[0].text)
        self.fermi_energy = float(root[1].text)
        self.energy_values = [float(root[2].text.split()[i]) for i in range(len(root[2].text.split()))]
        self.norbitals = len(root) - 3 # number of all orbitals in the unit cell

        self.orbitals = []
        for i in range(3, len(root)):
            self.orbitals.append(orbital(root[i]))

    def projected_to_element(self):
        tree = xml.etree.ElementTree.parse(self.file)
        root = tree.getroot()
        elements = {}
        for i in root[3:]:
            elements[i.attrib["species"]] = np.zeros(len(self.energy_values))

        for i in self.orbitals:
            elements[i.species] += i.data
        
        for element in elements:
            plt.plot(self.energy_values, elements[element], label=element)

        plt.vlines(self.fermi_energy, 0, 1, linestyle="dashed", label="Fermi Energy")
        plt.legend()
        #plt.show()
        plt.savefig("pdos.png")

    def export(self, directory):
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        self.projected_to_element()
        os.chdir("../../")

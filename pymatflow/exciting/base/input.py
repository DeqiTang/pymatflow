#from xml.etree.ElementTree import parse
import xml.etree.ElementTree as ET


input = ET.Element("input")


title = ET.Element("title", )
title.text = "test input"
structure = ET.Element("structure")
groundstate = ET.Element("groundstate")

groundstate.set("do", "from_scratch")
groundstate.set("ngridk", "4 4 4")


input.append(title)
input.append(structure)
input.append(groundstate)

input_tree = ET.ElementTree(input)
#input_tree.indent(space=" ")
input_tree.write("input.xml", method="xml", short_empty_elements=False)

class input:
    """
    """
    def __init__(self):
        pass

    def to_input(self, fout):
        pass

    def to_string(self):
        pass


    def set_params(self, params):
        """
        """
        pass
#!/usr/bin/env pyton
# _*_ coding: utf-8 _*_

import qmpy_rester as qr


def crystal_to_cartesian(cell, coord):
    """
    cell: [[x, x, x], [x, x, x], [x, x, x]]
    coord: [x, x, x] in crystal

    return [x, y, z] in cartesian
    """
    x = cell[0][0] * coord[0] + cell[1][0] * coord[1] + cell[2][0] * coord[2]
    y = cell[0][1] * coord[0] + cell[1][1] * coord[1] + cell[2][1] * coord[2]
    z = cell[0][2] * coord[0] + cell[1][2] * coord[1] + cell[2][2] * coord[2]
    return [x, y, z]


def oqmd_data_to_xyz(data, xyzfile="oqmd.xyz"):
    cell = data["unit_cell"]
    with open(xyzfile, 'w') as fout:
        fout.write("%d\n" % data['natoms'])
        fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0][0], cell[0][1], cell[0][2], cell[1][0], cell[1][1], cell[1][2], cell[2][0], cell[2][1], cell[2][2]))
        for atom in data['sites']:
            cart = crystal_to_cartesian(data["unit_cell"], [float(atom.split()[2]), float(atom.split()[3]), float(atom.split()[4])])
            fout.write(atom.split()[0])
            fout.write("\t%.9f\t%.9f\t%.9f\n" % (cart[0], cart[1], cart[2]))

def run():
    ## Return list of data
    with qr.QMPYRester() as q:
        kwargs = {
                'element_set': '(Li),H',  # composition include (Fe OR Mn) AND O
                                'stability': '0',            # hull distance smaller than -0.1 eV
                                        'natom': '<10',                  # number of atoms less than 10
                }
        list_of_data = q.get_oqmd_phases(verbose=False, **kwargs)
    #print(list_of_data)
    #print(list_of_data["data"][0])
    oqmd_data_to_xyz(list_of_data["data"][0])

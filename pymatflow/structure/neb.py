"""
implement inter image construction
"""

import copy
import numpy as np

def interpolate(initial, final, nimage, moving_atom):
    """
    :param initial -> instance of pymatflow.structure.crystal.crystal()
    :param final -> instance of pymatflow.structure.crystal.crystal()
    :param nimage -> number of intermediate images(not including initial and final image)
    :param moving_atom -> the index of moving atoms (index starts from 0)
    :return images -> a list of instance of pymatflow.structure.crystal.crystal() 
            as inter images
    Feature:
        only interpolate(linearly) the specified moving atoms. other atoms may have
        different positions between initial and final image, however we use the coordinate
        of the initial image as the cooresponding coordinate for the inter image.
    """
    # use fractional coordinates to interpolate the images
    initial.natom = len(initial.atoms)
    final.natom = len(final.atoms)
    initial_frac = initial.get_fractional()
    final_frac = final.get_fractional()

    #images = []
    images_frac = []
    for i in range(nimage):
        #images.append(copy.deepcopy(initial))
        images_frac.append(copy.deepcopy(initial_frac))

    for i in moving_atom:        
        dx = final_frac[i][1] - initial_frac[i][1]
        dy = final_frac[i][2] - initial_frac[i][2]
        dz = final_frac[i][3] - initial_frac[i][3]

        for j in range(nimage):
            x = initial_frac[i][1] + j * dx / (nimage+1)
            y = initial_frac[i][2] + j * dy / (nimage+1)            
            z = initial_frac[i][3] + j * dz / (nimage+1)
            #print("fractional x, y, z: %f %f %f\n" % (x, y, z))
            images_frac[j][i][1] = x
            images_frac[j][i][2] = y
            images_frac[j][i][3] = z
    # transform from images_frac to images
    # convert frac to cartesian again
    images = []
    latcell = np.array(initial.cell)
    convmat = latcell.T
    from pymatflow.base.atom import Atom
    from pymatflow.structure.crystal import Crystal
    for i in range(nimage):
        img = Crystal()
        img.atoms = []
        img.cell = initial.cell
        for atom in images_frac[i]:
            cartesian = list(convmat.dot(np.array([atom[1], atom[2], atom[3]])))
            img.atoms.append(Atom(name=atom[0], x=cartesian[0], y=cartesian[1], z=cartesian[2]))
        images.append(img)
        #
    return images
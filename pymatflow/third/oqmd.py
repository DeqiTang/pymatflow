"""
wrapper for OQMD RESTFUL API
"""
import os
import json
import numpy as np

import configparser

def _get_from_url(url):
    """
    Note: url is in format like this(following OQMD  RESTful API) and the result format should be set to json:
        http://oqmd.org/oqmdapi/formationenergy?fields=name,entry_id,delta_e&filter=stability=0&format=json

        Namely url should be in the form supported by OQMD RESTful API
    """
    os.system("mkdir -p /tmp/pymatflow/third")
    #os.system("wget \"%s\" -O /tmp/pymatflow/third/oqmd_restful_api_results.json" % (url))
    #os.system("curl \"%s\" -Lo /tmp/pymatflow/third/oqmd_restful_api_results.json" % (url))
    # silent output of curl
    os.system("curl \"%s\" -s -Lo /tmp/pymatflow/third/oqmd_restful_api_results.json" % (url))
    with open("/tmp/pymatflow/third/oqmd_restful_api_results.json", 'r') as fin:
        out = json.loads(fin.read())
    return out


def get_structure_from_url(url):
    from pymatflow.structure.crystal import crystal
    from pymatflow.base.atom import Atom

    result = _get_from_url(url)
    out = []
    for item in result["data"]:
        a = crystal()
        a.cell = item["unit_cell"]
        latcell = np.array(a.cell)
        convmat_frac_to_cartesian = latcell.T
        a.atoms = []
        for atm in item["sites"]:
            name = atm.split()[0]
            cartesian = list(convmat_frac_to_cartesian.dot(np.array([
                float(atm.split()[2]),
                float(atm.split()[3]),
                float(atm.split()[4])
                ])))
            a.atoms.append(Atom(name=name, x=cartesian[0], y=cartesian[1], z=cartesian[2]))
        out.append(a)
    return out

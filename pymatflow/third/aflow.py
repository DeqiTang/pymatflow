"""
wrapper for AFLOW RESTFUL API
AFLOW preovide REST-API and Search-API at the same time

References:
    * http://aflowlib.duke.edu/aflowwiki/
"""
import os
import json
import numpy as np

from urllib.request import urlopen
from pymatflow.structure.crystal import Crystal
from pymatflow.base.atom import Atom

def search_from_url(url):
    """
    Note: get search result from url directly through AFLOW Search API.
        URL example:
        http://aflowlib.duke.edu/search/API/?species((Na:K),Cl),nspecies(2),Egap(2*,*5),energy_cell,$paging(0)
    """
    response = json.loads(urlopen(url).read().decode("utf-8"))
    return response

def search(matchbook, directives="$paging(0)"):
    """
    :directives params: default "$paging(0)" means return all entries satisfying query
    :matchbook:
        example:
        "species((Na:K),Cl),nspecies(2),Egap(2*,*5),energy_cell"
    Note:
        search by accessing AFLUX search-API
    """
    MATCHBOOK = matchbook.replace(" ", "") # no space is allowed in the string for url
    #MATCHBOOK = matchbook
    DIRECTIVES = directives
    SERVER = "http://aflowlib.duke.edu"
    API = "/search/API/?"
    return search_from_url(SERVER + API + MATCHBOOK + "," + DIRECTIVES)


def get_val_from_aurl(aurl, key):
    """
    key :param ->  keyword
    Note: get result from aurl through AFLOW RESTful API.
        aflowlib.duke.edu:AFLOWDATA/LIB3_RAW/Bi_dRh_pvTi_sv/T0003.ABC:LDAU2
        aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/CUB/Cl1Na1_ICSD_622368

        run get_val_from_aurl(aurl="xxx", key="keywords")
        to see all the keywords supported to be passed to key
    """
    return urlopen("http://" + aurl.replace(":AFLOWDATA", "/AFLOWDATA") + "/?%s" % key).read().decode("utf-8")


def get_final_structure_from_aurl(aurl):
    """
    Note: get final structure(relaxed using vasp) from aurl, aurl is in format like this:
        aflowlib.duke.edu:AFLOWDATA/LIB3_RAW/Bi_dRh_pvTi_sv/T0003.ABC:LDAU2
        aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/CUB/Cl1Na1_ICSD_622368
    """
    out = crystal()

    
    geometry = urlopen("http://" + aurl.replace(":AFLOWDATA", "/AFLOWDATA") + "/?geometry").read().decode("utf-8")
    a, b, c, alpha, beta, gamma = [float(item) for item in geometry.replace("\n", "").split(",")]
    alpha = alpha / 180 * np.pi
    beta = beta / 180 * np.pi
    gamma = gamma / 180 * np.pi
    
    out.cell = []
    out.cell.append([a, 0, 0])
    out.cell.append([np.cos(gamma) * b, np.sin(gamma) * b, 0])
    new_c1 = np.cos(beta)
    new_c2 = (np.cos(alpha - np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
    new_c3_square = 1.0 - new_c1 * new_c1 - new_c2 * new_c2
    new_c3 = np.sqrt(new_c3_square)
    out.cell.append([new_c1 * c, new_c2 * c, new_c3 * c])

    elements_list = []
    composition = urlopen("http://" + aurl.replace(":AFLOWDATA", "/AFLOWDATA") + "/?composition").read().decode("utf-8")
    species = urlopen("http://" + aurl.replace(":AFLOWDATA", "/AFLOWDATA") + "/?species").read().decode("utf-8")
    for i, number in enumerate(composition.replace("\n", "").split(",")):
        for j in range(int(number)):
            elements_list.append(species.replace("\n", "").split(",")[i])
    
    #
    convmat_frac_to_cartesian = np.array(out.cell).T
    out.atoms = []
    positions_fractional = urlopen("http://" + aurl.replace(":AFLOWDATA", "/AFLOWDATA") + "/?positions_fractional").read().decode("utf-8")
    for i, item in enumerate(positions_fractional.replace("\n", "").split(";")):
        cartesian = convmat_frac_to_cartesian.dot([float(k) for k in item.split(",")])
        out.atoms.append(Atom(name=elements_list[i], x=cartesian[0], y=cartesian[1], z=cartesian[2]))
    #
    return out

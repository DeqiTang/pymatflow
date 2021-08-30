"""
wrapper for Materials Project RESTFUL API
"""
import os
import json

import requests
import configparser

def _get_from_uri(uri):
    """
    Note: uri is in format like this(following Materials Project RESTful API):
        https://www.materialsproject.org/rest/v2/materials/MoS2/vasp?API_KEY=XXXXXXXXXXXX
        or https://www.materialsproject.org/rest/v2/materials/MoS2/vasp/energy?API_KEY=XXXXXXX
        Namely uri should be in the form supported by Materials Project RESTful API, and the API_KEY
        should be contained.
    """
    os.system("mkdir -p /tmp/pymatflow/third")
    os.system("curl %s -s -Lo /tmp/pymatflow/third/mp_restful_api_results.json" % (uri))
    with open("/tmp/pymatflow/third/mp_restful_api_results.json", 'r') as fin:
        out = json.loads(fin.read())
    return out

def get_from_api(request_type, identifier, parameters=None):
    """
    :param request_type: can be materials, tasks, battery
    :param identifier: a materials composition, like Fe2O3
    """
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.expanduser("~"), ".pymatflow/mp.conf"))
    api_key = config["MaterialsProject"]["API_KEY"]
    
    if parameters == None:
        uri = "https://www.materialsproject.org/rest/v2/%s/%s/vasp?API_KEY=%s" % (request_type, identifier, api_key)
    else:
        uri = "https://www.materialsproject.org/rest/v2/%s/%s/vasp/%s?API_KEY=%s" % (request_type, identifier, parameters, api_key)
    
    return _get_from_uri(uri)

def get_from_post_query(data):
    """
    Note:
        data is in format like this:
        data = {
            'criteria': {
                    'elements': {'$in': ['Li', 'Na', 'K'], '$all': ['O']},
                            'nelements': 2,
            },
                'properties': [
                                'formula',
                                        'formation_energy_per_atom',                                   
                        ]
        }
    """
    import json
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.expanduser("~"), ".pymatflow/mp.conf"))
    api_key = config["MaterialsProject"]["API_KEY"]
    
    r = requests.post("https://materialsproject.org/rest/v2/query",
            headers={"X-API-KEY": api_key},
            data={k: json.dumps(v) for k, v in data.items()})
    return r.json()


def get_cif_from_api(request_type, identifier, directory="./", parameters=None):
    result = get_from_api(request_type=request_type, identifier=identifier, parameters=parameters)
    try:
        for item in result["response"]:
            with open(os.path.join(directory, "%s.cif" % item["material_id"]), "w") as fout:
                fout.write(item["cif"])
        successful = True
    except KeyError:
        successful = False
        print("KeyError in get 'response' key from json output of API provided by Materials Project")

    return successful

    #


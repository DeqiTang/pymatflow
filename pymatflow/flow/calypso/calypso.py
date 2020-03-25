"""
providing interface to calypso software
"""

class calypso:
    """
    """
    def __init__(self):
        self._initialize()

    def _initialize(self):
        self.params = {
            "SystemName": None, "NumberOfSpecies": None, "NameOfAtoms": None, "NumberOfFormula": None,
            "Volume": None, "Ialgo": None, "ICode": None, "NumberOfLocalOptim": None, "PsoRatio": None,
            "PopSize": None, "Kgrid": None, "Command": None, "MaxStep": None, "PickUp": None, "PickStep": None,
            "MaxTime": None, "LMC": None, "2D": None, "Areal": None, "DeltaZ": None, "MultiLayer": None,
            "LayerGap": None, "VacuumGap": None, "LAtom_Dis": None, "Cluster": None, "Vacancy": None,
            "MOL": None, "NumberOfTypeMolecule": None, "NumberOfMolecule": None, "DistOfMol": None,
            "SpeSpaceGroup": None, "FixCell": None, "FixAtom": None, "VSC": None, "MaxNumAtom": None,
            "LSurface": None, "CalSubstrate": None, "SurfaceThickness": None, "SubstrateThickness": None,
            "Substrate": None, "ECR": None, "CifFilePath": None, "MillerIndex": None, "SlabVacuumThick": None,
            "SlabDepth": None, "SlabNumLayers": None, "NumRelaxedLayers": None, "CapBondsWithH": None,
            "Hardness": None, "Band_edge": None, "TarBandGap": None, "Adsorption": None, "AdsorptionStyle": None,
            "NumberOfTypeAtom": None, "SuperCell": None, "RangeOfBondLenth": None, "AdsorptionSymmetry": None,
            "TypeOfSubstrate": None, "BothSide": None,
        }

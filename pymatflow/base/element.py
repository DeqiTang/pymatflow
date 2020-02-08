# -----------------------------------------------------------
# Information of the atomic_mass from:
# https://en.wikipedia.org/wiki/Periodic_table
# -----------------------------------------------------------

class Element:
    def __init__(self, number=None, mass=None, symbol=None):
        self.number = int(number) if number != None else None
        self.mass = float(mass) if mass != None else None
        self.symbol = str(symbol) if symbol != None else None


element = {
        "H":  Element(1, 1.008, "H"),
        "He": Element(2, 4.002602, "He"),
        "Li": Element(3, 6.94, "Li"),
        "Be": Element(4, 9.012183, "Be"),
        "B":  Element(5, 10.81, "B"),
        "C":  Element(6, 12.011, "C"),
        "N":  Element(7, 14.007, "N"),
        "O":  Element(8, 15.999, "O"),
        "F":  Element(9, 18.99840316, "F"),
        "Ne": Element(10, 20.1797, "Ne"),
        "Na": Element(11, 22.98976928, "Na"),
        "Mg": Element(12, 24.305, "Mg"),
        "Al": Element(13, 26.9815384, "Al"),
        "Si": Element(14, 28.085, "Si"),
        "P":  Element(15, 30.973761998, "P"),
        "S":  Element(16, 32.06, "S"),
        "Cl": Element(17, 35.45, "Cl"),
        "Ar": Element(18, 39.95, "Ar"),
        "K":  Element(19, 39.0983, "K"),
        "Ca": Element(20, 40.078, "Ca"),
        "Sc": Element(21, 44.955908, "Sc"),
        "Ti": Element(22, 47.867, "Ti"),
        "V":  Element(23, 50.9415, "V"),
        "Cr": Element(24, 51.9961, "Cr"),
        "Mn": Element(25, 54.938043, "Mn"),
        "Fe": Element(26, 55.845, "Fe"),
        "Co": Element(27, 58.933154, "Co"),
        "Ni": Element(28, 58.6934, "Ni"),
        "Cu": Element(29, 63.546, "Cu"),
        "Zn": Element(30, 65.38, "Zn"),
        "Ga": Element(31, 69.723, "Ga"),
        "Ge": Element(32, 72.630, "Ge"),
        "As": Element(33, 74.921595, "As"),
        "Se": Element(34, 78.971, "Se"),
        "Br": Element(35, 79.904, "Br"),
        "Kr": Element(36, 83.798, "Kr"),
        "Rb": Element(37, 85.4678, "Rb"),
        "Sr": Element(38, 87.62, "Sr"),
        "Y":  Element(39, 88.90584, "Y"),
        "Zr": Element(40, 91.224, "Zr"),
        "Nb": Element(41, 92.90637, "Nb"),
        "Mo": Element(42, 95.95, "Mo"),
        "Tc": Element(43, 97.0, "Tc"),           # mass number
        "Ru": Element(44, 101.07, "Ru"),
        "Rh": Element(45, 102.90549, "Rh"),
        "Pd": Element(46, 106.42, "Pd"),
        "Ag": Element(47, 107.8682, "Ag"),
        "Cd": Element(48, 112.414, "Cd"),
        "In": Element(49, 114.818, "In"),
        "Sn": Element(50, 118.710, "Sn"),
        "Sb": Element(51, 121.760, "Sb"),
        "Te": Element(52, 127.60, "Te"),
        "I":  Element(53, 126.90447, "I"),
        "Xe": Element(54, 131.293, "Xe"),
        "Cs": Element(55, 132.90545196, "Cs"),
        "Ba": Element(56, 137.327, "Ba"),
        "La": Element(57, 138.90547, "La"),
        "Ce": Element(58, 140.116, "Ce"),
        "Pr": Element(59, 140.90766, "Pr"),
        "Nd": Element(60, 144.242, "Nd"),
        "Pm": Element(61, 145.0, "Pm"),          # mass number
        "Sm": Element(62, 150.36, "Sm"),
        "Eu": Element(63, 151.964, "Eu"),
        "Gd": Element(64, 157.25, "Gd"),
        "Tb": Element(65, 158.925354, "Tb"),
        "Dy": Element(66, 162.500, "Dy"),
        "Ho": Element(67, 164.930328, "Ho"),
        "Er": Element(68, 167.259, "Er"),
        "Tm": Element(69, 168.934218, "Tm"),
        "Yb": Element(70, 173.045, "Yb"),
        "Lu": Element(71, 174.9668, "Lu"),
        "Hf": Element(72, 178.49, "Hf"),
        "Ta": Element(73, 180.94788, "Ta"),
        "W":  Element(74, 183.84, "W"),
        "Re": Element(75, 186.207, "Re"),
        "Os": Element(76, 190.23, "Os"),
        "Ir": Element(77, 192.217, "Os"),
        "Pt": Element(78, 195.084, "Pt"),
        "Au": Element(79, 196.966570, "Au"),
        "Hg": Element(80, 200.592, "Hg"),
        "Tl": Element(81, 204.38, "Tl"),
        "Pb": Element(82, 207.2, "Pb"),
        "Bi": Element(83, 208.98040, "Bi"),
        "Po": Element(84, 209, "Po"),           # mass number
        "At": Element(85, 210, "At"),           # mass number
        "Rn": Element(86, 222, "Rn"),           # mass number
        "Fr": Element(87, 223, "Fr"),           # mass number
        "Ra": Element(88, 226, "Ra"),           # mass number
        "Ac": Element(89, 227, "Ac"),           # mass number
        "Th": Element(90, 232.0377, "Th"),
        "Pa": Element(91, 231.03588, "Pa"),
        "U":  Element(92, 238.02891, "U"),
        "Np": Element(93, 237, "Np"),
        "Pu": Element(94, 244, "Pu"),           # mass number
        "Am": Element(95, 243, "Am"),           # mass number
        "Cm": Element(96, 247, "Cm"),           # mass number
        "Bk": Element(97, 247, "Bk"),           # mass number
        "Cf": Element(98, 251, "Cf"),           # mass number
        "Es": Element(99, 252, "Es"),           # mass number
        "Fm": Element(100, 257, "Fm"),          # mass number
        "Md": Element(101, 258, "Md"),          # mass number
        "No": Element(102, 259, "No"),          # mass number
        "Lr": Element(103, 266, "Lr"),          # mass number
        "Rf": Element(104, 267, "Rf"),          # mass number
        "Db": Element(105, 268, "Db"),          # mass number
        "Sg": Element(106, 269, "Sg"),          # mass number
        "Bh": Element(107, 270, "Bh"),          # mass number, unconfirmed: 278
        "Hs": Element(108, 269, "Hs"),          # mass number
        "Mt": Element(109, 278, "Mt"),          # mass number, unconfirmed: 282
        "Ds": Element(110, 281, "Ds"),          # mass number
        "Rg": Element(111, 282, "Rg"),          # mass number, unconfirmed: 286
        "Cn": Element(112, 285, "Cn"),          # mass number
        "Nh": Element(113, 286, "Nh"),          # mass number
        "Fl": Element(114, 289, "Fl"),          # mass number, unconfirmed: 290
        "Mc": Element(115, 290, "Mc"),          # mass number
        "Lv": Element(116, 293, "Lv"),          # mass number
        "Ts": Element(117, 294, "Ts"),          # mass number
        "Og": Element(118, 294, "Og"),          # mass number, unconfirmed: 295
        }

/*
-----------------------------------------------------------
Information of the atomic_mass from:
https://en.wikipedia.org/wiki/Periodic_table
-----------------------------------------------------------
*/

#ifndef ATOMSCIKIT_INCLUDE_ASKIT_BASE_ELEMENT_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_BASE_ELEMENT_H_

#include <map>
#include <string>

namespace askit {

class Element {
  public:
    Element(int number, double mass, std::string symbol): number(number), mass(mass), symbol(symbol) {}
    Element();
    ~Element();
    int number;
    double mass;
    std::string symbol;
};





inline Element build_element(int number, double mass, std::string symbol) {
  Element element(number, mass, symbol);
  return element;
}

inline std::map<std::string, Element> get_element_number_map() {
  std::map<std::string, Element> element_number_map;
  element_number_map["H"] =  build_element(1, 1.008, "H");
  element_number_map["He"] = build_element(2, 4.002602, "He");
  element_number_map["Li"] = build_element(3, 6.94, "Li");
  element_number_map["Be"] = build_element(4, 9.012183, "Be");
  element_number_map["B"] =  build_element(5, 10.81, "B");
  element_number_map["C"] =  build_element(6, 12.011, "C");
  element_number_map["N"] =  build_element(7, 14.007, "N");
  element_number_map["O"] =  build_element(8, 15.999, "O");
  element_number_map["F"] =  build_element(9, 18.99840316, "F");
  element_number_map["Ne"] = build_element(10, 20.1797, "Ne");
  element_number_map["Na"] = build_element(11, 22.98976928, "Na");
  element_number_map["Mg"] = build_element(12, 24.305, "Mg");
  element_number_map["Al"] = build_element(13, 26.9815384, "Al");
  element_number_map["Si"] = build_element(14, 28.085, "Si");
  element_number_map["P"] =  build_element(15, 30.973761998, "P");
  element_number_map["S"] =  build_element(16, 32.06, "S");
  element_number_map["Cl"] = build_element(17, 35.45, "Cl");
  element_number_map["Ar"] = build_element(18, 39.95, "Ar");
  element_number_map["K"] =  build_element(19, 39.0983, "K");
  element_number_map["Ca"] = build_element(20, 40.078, "Ca");
  element_number_map["Sc"] = build_element(21, 44.955908, "Sc");
  element_number_map["Ti"] = build_element(22, 47.867, "Ti");
  element_number_map["V"] =  build_element(23, 50.9415, "V");
  element_number_map["Cr"] = build_element(24, 51.9961, "Cr");
  element_number_map["Mn"] = build_element(25, 54.938043, "Mn");
  element_number_map["Fe"] = build_element(26, 55.845, "Fe");
  element_number_map["Co"] = build_element(27, 58.933154, "Co");
  element_number_map["Ni"] = build_element(28, 58.6934, "Ni");
  element_number_map["Cu"] = build_element(29, 63.546, "Cu");
  element_number_map["Zn"] = build_element(30, 65.38, "Zn");
  element_number_map["Ga"] = build_element(31, 69.723, "Ga");
  element_number_map["Ge"] = build_element(32, 72.630, "Ge");
  element_number_map["As"] = build_element(33, 74.921595, "As");
  element_number_map["Se"] = build_element(34, 78.971, "Se");
  element_number_map["Br"] = build_element(35, 79.904, "Br");
  element_number_map["Kr"] = build_element(36, 83.798, "Kr");
  element_number_map["Rb"] = build_element(37, 85.4678, "Rb");
  element_number_map["Sr"] = build_element(38, 87.62, "Sr");
  element_number_map["Y"] =  build_element(39, 88.90584, "Y");
  element_number_map["Zr"] = build_element(40, 91.224, "Zr");
  element_number_map["Nb"] = build_element(41, 92.90637, "Nb");
  element_number_map["Mo"] = build_element(42, 95.95, "Mo");
  element_number_map["Tc"] = build_element(43, 97.0, "Tc");          // mass number
  element_number_map["Ru"] = build_element(44, 101.07, "Ru");
  element_number_map["Rh"] = build_element(45, 102.90549, "Rh");
  element_number_map["Pd"] = build_element(46, 106.42, "Pd");
  element_number_map["Ag"] = build_element(47, 107.8682, "Ag");
  element_number_map["Cd"] = build_element(48, 112.414, "Cd");
  element_number_map["In"] = build_element(49, 114.818, "In");
  element_number_map["Sn"] = build_element(50, 118.710, "Sn");
  element_number_map["Sb"] = build_element(51, 121.760, "Sb");
  element_number_map["Te"] = build_element(52, 127.60, "Te");
  element_number_map["I"] =  build_element(53, 126.90447, "I");
  element_number_map["Xe"] = build_element(54, 131.293, "Xe");
  element_number_map["Cs"] = build_element(55, 132.90545196, "Cs");
  element_number_map["Ba"] = build_element(56, 137.327, "Ba");
  element_number_map["La"] = build_element(57, 138.90547, "La");
  element_number_map["Ce"] = build_element(58, 140.116, "Ce");
  element_number_map["Pr"] = build_element(59, 140.90766, "Pr");
  element_number_map["Nd"] = build_element(60, 144.242, "Nd");
  element_number_map["Pm"] = build_element(61, 145.0, "Pm");          // mass number
  element_number_map["Sm"] = build_element(62, 150.36, "Sm");
  element_number_map["Eu"] = build_element(63, 151.964, "Eu");
  element_number_map["Gd"] = build_element(64, 157.25, "Gd");
  element_number_map["Tb"] = build_element(65, 158.925354, "Tb");
  element_number_map["Dy"] = build_element(66, 162.500, "Dy");
  element_number_map["Ho"] = build_element(67, 164.930328, "Ho");
  element_number_map["Er"] = build_element(68, 167.259, "Er");
  element_number_map["Tm"] = build_element(69, 168.934218, "Tm");
  element_number_map["Yb"] = build_element(70, 173.045, "Yb");
  element_number_map["Lu"] = build_element(71, 174.9668, "Lu");
  element_number_map["Hf"] = build_element(72, 178.49, "Hf");
  element_number_map["Ta"] = build_element(73, 180.94788, "Ta");
  element_number_map["W"] =  build_element(74, 183.84, "W");
  element_number_map["Re"] = build_element(75, 186.207, "Re");
  element_number_map["Os"] = build_element(76, 190.23, "Os");
  element_number_map["Ir"] = build_element(77, 192.217, "Os");
  element_number_map["Pt"] = build_element(78, 195.084, "Pt");
  element_number_map["Au"] = build_element(79, 196.966570, "Au");
  element_number_map["Hg"] = build_element(80, 200.592, "Hg");
  element_number_map["Tl"] = build_element(81, 204.38, "Tl");
  element_number_map["Pb"] = build_element(82, 207.2, "Pb");
  element_number_map["Bi"] = build_element(83, 208.98040, "Bi");
  element_number_map["Po"] = build_element(84, 209, "Po");           // mass number
  element_number_map["At"] = build_element(85, 210, "At");           // mass number
  element_number_map["Rn"] = build_element(86, 222, "Rn");           // mass number
  element_number_map["Fr"] = build_element(87, 223, "Fr");           // mass number
  element_number_map["Ra"] = build_element(88, 226, "Ra");           // mass number
  element_number_map["Ac"] = build_element(89, 227, "Ac");           // mass number
  element_number_map["Th"] = build_element(90, 232.0377, "Th");
  element_number_map["Pa"] = build_element(91, 231.03588, "Pa");
  element_number_map["U"] =  build_element(92, 238.02891, "U");
  element_number_map["Np"] = build_element(93, 237, "Np");
  element_number_map["Pu"] = build_element(94, 244, "Pu");           // mass number
  element_number_map["Am"] = build_element(95, 243, "Am");           // mass number
  element_number_map["Cm"] = build_element(96, 247, "Cm");           // mass number
  element_number_map["Bk"] = build_element(97, 247, "Bk");           // mass number
  element_number_map["Cf"] = build_element(98, 251, "Cf");           // mass number
  element_number_map["Es"] = build_element(99, 252, "Es");           // mass number
  element_number_map["Fm"] = build_element(100, 257, "Fm");          // mass number
  element_number_map["Md"] = build_element(101, 258, "Md");          // mass number
  element_number_map["No"] = build_element(102, 259, "No");          // mass number
  element_number_map["Lr"] = build_element(103, 266, "Lr");          // mass number
  element_number_map["Rf"] = build_element(104, 267, "Rf");          // mass number
  element_number_map["Db"] = build_element(105, 268, "Db");          // mass number
  element_number_map["Sg"] = build_element(106, 269, "Sg");          // mass number
  element_number_map["Bh"] = build_element(107, 270, "Bh");          // mass number, unconfirmed] = 278
  element_number_map["Hs"] = build_element(108, 269, "Hs");          // mass number
  element_number_map["Mt"] = build_element(109, 278, "Mt");          // mass number, unconfirmed] = 282
  element_number_map["Ds"] = build_element(110, 281, "Ds");          // mass number
  element_number_map["Rg"] = build_element(111, 282, "Rg");          // mass number, unconfirmed] = 286
  element_number_map["Cn"] = build_element(112, 285, "Cn");          // mass number
  element_number_map["Nh"] = build_element(113, 286, "Nh");          // mass number
  element_number_map["Fl"] = build_element(114, 289, "Fl");          // mass number, unconfirmed] = 290
  element_number_map["Mc"] = build_element(115, 290, "Mc");          // mass number
  element_number_map["Lv"] = build_element(116, 293, "Lv");          // mass number
  element_number_map["Ts"] = build_element(117, 294, "Ts");          // mass number
  element_number_map["Og"] = build_element(118, 294, "Og");          // mass number, unconfirmed] = 295

  return element_number_map;
}

} // end namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_H_
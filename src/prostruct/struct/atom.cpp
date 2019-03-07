/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "prostruct/struct/atom.h"
#include "atom.h"

#include <algorithm>
#include <memory>

static const std::map<std::string, std::vector<std::string>> elementName
{
	//       name
	{"H" , {"Hydrogen",    }},
	{"He", {"Helium",      }},
	{"Li", {"Lithium",     }},
	{"Be", {"Beryllium",   }},
	{"B" , {"Boron",       }},
	{"C" , {"Carbon",      }},
	{"N" , {"Nitrogen",    }},
	{"O" , {"Oxygen",      }},
	{"F" , {"Fluorine",    }},
	{"Ne", {"Neon",        }},
	{"Na", {"Sodium",      }},
	{"Mg", {"Magnesium",   }},
	{"Al", {"Aluminium",   }},
	{"Si", {"Silicon",     }},
	{"P" , {"Phosphorus",  }},
	{"S" , {"Sulphur",     }},
	{"Cl", {"Chlorine",    }},
	{"Ar", {"Argon",       }},
	{"K" , {"Potassium",   }},
	{"Ca", {"Calcium",     }},
	{"Sc", {"Scandium",    }},
	{"Ti", {"Titanium",    }},
	{"V" , {"Vanadium",    }},
	{"Cr", {"Chromium",    }},
	{"Mn", {"Manganese",   }},
	{"Fe", {"Iron",        }},
	{"Co", {"Cobalt",      }},
	{"Ni", {"Nickel",      }},
	{"Cu", {"Copper",      }},
	{"Zn", {"Zinc",        }},
	{"Ga", {"Gallium",     }},
	{"Ge", {"Germanium",   }},
	{"As", {"Arsenic",     }},
	{"Se", {"Selenium",    }},
	{"Br", {"Bromine",     }},
	{"Kr", {"Krypton",     }},
	{"Rb", {"Rubidium",    }},
	{"Sr", {"Strontium",   }},
	{"Y" , { "Yttrium",    }},
	{"Zr", {"Zirconium",   }},
	{"Nb", {"Niobium",     }},
	{"Mo", {"Molybdenum",  }},
	{"Tc", {"Technetium",  }},
	{"Ru", {"Ruthenium",   }},
	{"Rh", {"Rhodium",     }},
	{"Pd", {"Palladium",   }},
	{"Ag", {"Silver",      }},
	{"Cd", {"Cadmium",     }},
	{"In", {"Indium",      }},
	{"Sn", {"Tin",         }},
	{"Sb", {"Antimony",    }},
	{"Te", {"Tellurium",   }},
	{"I" , { "Iodine",     }},
	{"Xe", {"Xenon",       }},
	{"Cs", {"Cesium",      }},
	{"Ba", {"Barium",      }},
	{"La", {"Lanthanum",   }},
	{"Ce", {"Cerium",      }},
	{"Pr", {"Praseodymium",}},
	{"Nd", {"Neodymium",   }},
	{"Pm", {"Promethium",  }},
	{"Sm", {"Samarium",    }},
	{"Eu", {"Europium",    }},
	{"Gd", {"Gadolinium",  }},
	{"Tb", {"Terbium",     }},
	{"Dy", {"Dysprosium",  }},
	{"Ho", {"Holmium",     }},
	{"Er", {"Erbium",      }},
	{"Tm", {"Thulium",     }},
	{"Yb", {"Ytterbium",   }},
	{"Lu", {"Lutetium",    }},
	{"Hf", {"Hafnium",     }},
	{"Ta", {"Tantalum",    }},
	{"W" , {"Tungsten",    }},
	{"Re", {"Rhenium",     }},
	{"Os", {"Osmium",      }},
	{"Ir", {"Iridium",     }},
	{"Pt", {"Platinum",    }},
	{"Au", {"Gold",        }},
	{"Hg", {"Mercury",     }},
	{"Tl", {"Thallium",    }},
	{"Pb", {"Lead",        }},
	{"Bi", {"Bismuth",     }},
	{"Po", {"Polonium",    }},
	{"At", {"Astatine",    }},
	{"Rn", {"Radon",       }},
	{"Fr", {"Francium",    }},
	{"Ra", {"Radium",      }},
	{"Ac", {"Actinium",    }},
	{"Th", {"Thorium",     }},
	{"Pa", {"Protactinium",}},
	{"U" , {"Uranium",     }},
	{"Np", {"Neptunium",   }},
	{"Pu", {"Plutonium",   }},
	{"Am", {"Americium",   }},
	{"Cm", {"Curium",      }},
	{"Bk", {"Berkelium",   }},
	{"Cf", {"Californium", }},
	{"Es", {"Einsteinium", }},
	{"Fm", {"Fermium",     }},
	{"Md", {"Mendelevium", }},
	{"No", {"Nobelium",    }},
	{"Lr", {"Lawrencium",  }},
	{"Rf", {"Rutherfordium"}},
	{"Db", {"Dubnium",     }},
	{"Sg", {"Seaborgium",  }},
	{"Bh", {"Bohrium",     }},
	{"Hs", {"Hassium",     }},
	{"Mt", {"Meitnerium",  }},
	{"Ds", {"Darmstadtium",}},
	{"Rg", {"Roentgenium", }},
	{"Cn", {"Copernicium", }},
	{"Ut", {"Ununtrium",   }},
	{"Fl", {"Flerovium",   }},
	{"Up", {"Ununpentium", }},
	{"Lv", {"Livermorium", }},
	{"Us", {"Ununseptium", }},
	{"Uo", {"Ununoctium",  }}
};

static const std::map<std::string, std::vector<double>> elementDescription {
	// atomic number, atomic weight
	{ "H", { 1, 1.0079 } },
	{ "He", { 2, 4.0026 } },
	{ "Li", { 3, 6.941 } },
	{ "Be", { 4, 9.0122 } },
	{ "B", { 5, 10.811 } },
	{ "C", { 6, 12.0107 } },
	{ "N", { 7, 14.0067 } },
	{ "O", { 8, 15.9994 } },
	{ "F", { 9, 18.9984 } },
	{ "Ne", { 10, 20.1797 } },
	{ "Na", { 11, 22.9897 } },
	{ "Mg", { 12, 24.305 } },
	{ "Al", { 13, 26.9815 } },
	{ "Si", { 14, 28.0855 } },
	{ "P", { 15, 30.9738 } },
	{ "S", { 16, 32.065 } },
	{ "Cl", { 17, 35.453 } },
	{ "Ar", { 18, 39.948 } },
	{ "K", { 19, 39.0983 } },
	{ "Ca", { 20, 40.078 } },
	{ "Sc", { 21, 44.9559 } },
	{ "Ti", { 22, 47.867 } },
	{ "V", { 23, 50.9415 } },
	{ "Cr", { 24, 51.9961 } },
	{ "Mn", { 25, 54.938 } },
	{ "Fe", { 26, 55.845 } },
	{ "Co", { 27, 58.9332 } },
	{ "Ni", { 28, 58.6934 } },
	{ "Cu", { 29, 63.546 } },
	{ "Zn", { 30, 65.39 } },
	{ "Ga", { 31, 69.723 } },
	{ "Ge", { 32, 72.64 } },
	{ "As", { 33, 74.9216 } },
	{ "Se", { 34, 78.96 } },
	{ "Br", { 35, 79.904 } },
	{ "Kr", { 36, 83.8 } },
	{ "Rb", { 37, 85.4678 } },
	{ "Sr", { 38, 87.62 } },
	{ "Y", { 39, 88.9059 } },
	{ "Zr", { 40, 91.224 } },
	{ "Nb", { 41, 92.9064 } },
	{ "Mo", { 42, 95.94 } },
	{ "Tc", { 43, 98 } },
	{ "Ru", { 44, 101.07 } },
	{ "Rh", { 45, 102.9055 } },
	{ "Pd", { 46, 106.42 } },
	{ "Ag", { 47, 107.8682 } },
	{ "Cd", { 48, 112.411 } },
	{ "In", { 49, 114.818 } },
	{ "Sn", { 50, 118.71 } },
	{ "Sb", { 51, 121.76 } },
	{ "Te", { 52, 127.6 } },
	{ "I", { 53, 126.9045 } },
	{ "Xe", { 54, 131.293 } },
	{ "Cs", { 55, 132.9055 } },
	{ "Ba", { 56, 137.327 } },
	{ "La", { 57, 138.9055 } },
	{ "Ce", { 58, 140.116 } },
	{ "Pr", { 59, 140.9077 } },
	{ "Nd", { 60, 144.24 } },
	{ "Pm", { 61, 145 } },
	{ "Sm", { 62, 150.36 } },
	{ "Eu", { 63, 151.964 } },
	{ "Gd", { 64, 157.25 } },
	{ "Tb", { 65, 158.9253 } },
	{ "Dy", { 66, 162.5 } },
	{ "Ho", { 67, 164.9303 } },
	{ "Er", { 68, 167.259 } },
	{ "Tm", { 69, 168.9342 } },
	{ "Yb", { 70, 173.04 } },
	{ "Lu", { 71, 174.967 } },
	{ "Hf", { 72, 178.49 } },
	{ "Ta", { 73, 180.9479 } },
	{ "W", { 74, 183.84 } },
	{ "Re", { 75, 186.207 } },
	{ "Os", { 76, 190.23 } },
	{ "Ir", { 77, 192.217 } },
	{ "Pt", { 78, 195.078 } },
	{ "Au", { 79, 196.9665 } },
	{ "Hg", { 80, 200.59 } },
	{ "Tl", { 81, 204.3833 } },
	{ "Pb", { 82, 207.2 } },
	{ "Bi", { 83, 208.9804 } },
	{ "Po", { 84, 209 } },
	{ "At", { 85, 210 } },
	{ "Rn", { 86, 222 } },
	{ "Fr", { 87, 223 } },
	{ "Ra", { 88, 226 } },
	{ "Ac", { 89, 227 } },
	{ "Th", { 90, 232.0381 } },
	{ "Pa", { 91, 231.0359 } },
	{ "U", { 92, 238.0289 } },
	{ "Np", { 93, 237 } },
	{ "Pu", { 94, 244 } },
	{ "Am", { 95, 243 } },
	{ "Cm", { 96, 247 } },
	{ "Bk", { 97, 247 } },
	{ "Cf", { 98, 251 } },
	{ "Es", { 99, 252 } },
	{ "Fm", { 100, 257 } },
	{ "Md", { 101, 258 } },
	{ "No", { 102, 259 } },
	{ "Lr", { 103, 262 } },
	{ "Rf", { 104, 261 } },
	{ "Db", { 105, 262 } },
	{ "Sg", { 106, 263 } },
	{ "Bh", { 107, 262 } },
	{ "Hs", { 108, 265 } },
	{ "Mt", { 109, 266 } },
	{ "Ds", { 110, 281 } },
	{ "Rg", { 111, 272 } },
	{ "Cn", { 112, 285 } },
	{ "Ut", { 113, 284 } },
	{ "Fl", { 114, 289 } },
	{ "Up", { 115, 288 } },
	{ "Lv", { 116, 292 } },
	{ "Us", { 117, 292 } },
	{ "Uo", { 118, 294 } }
};

template <typename T>
void Atom<T>::load_atom(const std::string& element_)
{

	/// Private method of Atom to load all the information
	/// of the given atom if it exists
	/// @param [std::string] name Name of the atom.

	element = element_;

	// check if it is a valid element
	if (elementDescription.find(element) == elementDescription.end())
		throw "Unknown element: " + element;

	atomicNumber = static_cast<int>(elementDescription.at(element)[0]);
	atomicWeight = elementDescription.at(element)[1];
	name = elementName.at(element)[0];
}

template <typename T>
void Atom<T>::load_atom(const std::string& element_, const std::string& name_)
{

	/// Private method of Atom to load all the information
	/// of the given atom if it exists
	/// @param [std::string] name Name of the atom.

	element = element_;

	// check if it is a valid element
	if (elementDescription.find(element) == elementDescription.end())
		throw "Unknown element: " + element;

	atomicNumber = static_cast<int>(elementDescription.at(element)[0]);
	atomicWeight = elementDescription.at(element)[1];
	name = name_;
}

template <typename T>
void Atom<T>::load_atom(const std::string& element_, const std::string& name_, T x_, T y_, T z_)
{

	/// Private method of Atom to load all the information
	/// of the given atom if it exists
	/// @param [std::string] name Name of the atom.

	element = element_;

	// check if it is a valid element
	if (elementDescription.find(element) == elementDescription.end())
		throw "Unknown element: " + element;

	atomicNumber = static_cast<int>(elementDescription.at(element)[0]);
	atomicWeight = elementDescription.at(element)[1];
	name = name_;
	x = x_;
	y = y_;
	z = z_;
}

template <typename T>
void Atom<T>::addBond(std::shared_ptr<Atom<T>> atom, int atomType)
{

	// assumes that an atom can at most form 4 bonds
	if (bonds.size() < 4) {
		auto newBond = std::make_shared<Bond<T>>(getAtom(), atom, atomType);
		bonds.emplace_back(newBond);
		// adds Bond to second Atom
		atom->addBond(newBond);
	} else
		throw "Atom can form at most 4 bonds";
}

template <typename T>
void Atom<T>::addBond(std::shared_ptr<Bond<T>> bond)
{

	// assumes that an atom can at most form 4 bonds
	if (bonds.size() < 4) {
		bonds.emplace_back(bond);
	} else
		throw "Atom can form at most 4 bonds";
}

template <typename T>
void Atom<T>::destroyBond(int bondIndex)
{
	// destroy reference to this bond from pairing atom
	bonds[bondIndex]->getAtom1()->destroyBond(bonds[bondIndex]);
	// destroy reference to this bond from this atom
	bonds.erase(bonds.begin() + bondIndex);
}

template <typename T>
void Atom<T>::destroyBond(std::shared_ptr<Bond<T>> bondP)
{

	auto iter = std::find(bonds.begin(), bonds.end(), bondP);

	if (iter == bonds.end()) {
		throw "The given bond pointer was not found in this atom";
	}

	bonds.erase(iter);
}

template <typename T>
bool Atom<T>::hasBond(const std::shared_ptr<Atom<T>>& atom2)
{

	// Checks if there is a bond between this and atom2
	// Note that the bond has no direction, therefore need
	// to check atom2 in Bond class member Atom1 and Atom2
	for (auto const& bond : bonds) {
		if (bond->getAtom1() == atom2 or bond->getAtom2() == atom2)
			return true;
	}
	return false;
}

template class Atom<float>;
template class Atom<double>;
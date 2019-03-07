//
// Created by gil on 28/03/18.
//

#include "gtest/gtest.h"

#include "prostruct/struct/residue.h"

using namespace prostruct;

template <typename T>
class ResidueTestClass: public Residue<T>
{
public:
	template <class... Args>
	ResidueTestClass(Args... m): Residue<T>(m...) {};

	auto get_atom_indices(const std::string& name) {
		return Residue<T>::get_atom_indices(name);
	}
};

TEST(ResidueTest, Alanine) {

	auto N = std::make_shared<Atom<double>>("N", "N", 3.653, 52.210, 0.963);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 32.975, 51.134, 0.036);
	auto C = std::make_shared<Atom<double>>("C", "C", 33.489, 49.933, 0.818);
	auto O = std::make_shared<Atom<double>>("O", "O", 34.699, 49.738, 0.953);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 34.026, 51.601, -0.964);

	auto ala = ResidueTestClass<double>(atomVector<double>({N, CA, C, O, CB}), "ALA", "ALA1");

	ASSERT_EQ(CA->getNumberOfBonds(), 3);
	ASSERT_EQ(N->getNumberOfBonds(), 1);
	CA->destroyBond(2);

	ASSERT_EQ(CA->getNumberOfBonds(), 2);
	ASSERT_EQ(CB->getNumberOfBonds(), 0);

	ala.createBonds();
	ASSERT_EQ(CA->getNumberOfBonds(), 3);
	ASSERT_EQ(CB->getNumberOfBonds(), 1);

	ASSERT_EQ(ala.get_atom_indices("N")(0), 0);
	ASSERT_EQ(ala.get_atom_indices("N").size(), 1);
	ASSERT_EQ(ala.get_atom_indices("C")(0), 2);
	ASSERT_EQ(ala.get_atom_indices("C").size(), 1);
}

TEST(ResidueTest, Argenine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 32.964, 52.298, 5.433);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 31.521, 52.533, 5.366);
	auto C = std::make_shared<Atom<double>>("C", "C", 31.011, 52.774, 3.947);
	auto O = std::make_shared<Atom<double>>("O", "O", 30.021, 52.173, 3.525);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 31.143, 53.743, 6.231);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 31.514, 53.618, 7.696);
	auto CD = std::make_shared<Atom<double>>("C", "CD", 31.171, 54.883, 8.482);
	auto NE = std::make_shared<Atom<double>>("N", "NE", 29.739, 55.187, 8.480);
	auto CZ = std::make_shared<Atom<double>>("C", "CZ", 29.141, 56.009, 7.622);
	auto NH1 = std::make_shared<Atom<double>>("N", "NH1", 29.847, 56.622, 6.682);
	auto NH2 = std::make_shared<Atom<double>>("N", "NH2", 27.832, 56.219, 7.704);

	auto arg = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD, CD, NE, CZ, NH1, NH2 }), "ARG", "ARG1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(CD->getNumberOfBonds(), 2);
	ASSERT_EQ(NE->getNumberOfBonds(), 2);
	ASSERT_EQ(CZ->getNumberOfBonds(), 3);
	ASSERT_EQ(NH1->getNumberOfBonds(), 1);
	ASSERT_EQ(NH2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Asparagine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 3.649, 51.692, -4.274);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 3.079, 51.607, -5.619);
	auto C = std::make_shared<Atom<double>>("C", "C", 1.960, 52.621, -5.867);
	auto O = std::make_shared<Atom<double>>("O", "O", 1.326, 52.628, -6.928);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 2.586, 50.183, -5.914);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 1.469, 49.736, -4.987);
	auto OD1 = std::make_shared<Atom<double>>("O", "OD1", 1.048, 50.471, -4.091);
	auto ND2 = std::make_shared<Atom<double>>("N", "ND2", 0.980, 48.518, -5.205);

	auto asn = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, OD1, ND2 }), "ASN", "ASN1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(OD1->getNumberOfBonds(), 1);
	ASSERT_EQ(ND2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Aspartate)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -6.741, 72.216, -0.006);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -6.999, 72.347, 1.416);
	auto C = std::make_shared<Atom<double>>("C", "C", -7.347, 73.805, 1.685);
	auto O = std::make_shared<Atom<double>>("O", "O", -7.709, 74.541, 0.764);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -5.769, 71.916, 2.225);
	auto CG = std::make_shared<Atom<double>>("C", "CG", -6.134, 71.383, 3.599);
	auto OD1 = std::make_shared<Atom<double>>("O", "OD1", -6.578, 72.177, 4.455);
	auto OD2 = std::make_shared<Atom<double>>("O", "OD2", -5.989, 70.163, 3.821);

	auto asp = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, OD1, OD2 }), "ASP", "ASP1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(OD1->getNumberOfBonds(), 1);
	ASSERT_EQ(OD2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Cysteine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -21.992, 79.723, 2.270);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -23.340, 79.793, 2.814);
	auto C = std::make_shared<Atom<double>>("C", "C", -23.848, 81.236, 2.905);
	auto O = std::make_shared<Atom<double>>("O", "O", -24.993, 81.514, 2.556);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -23.385, 79.131, 4.192);
	auto SG = std::make_shared<Atom<double>>("S", "SG", -23.535, 77.308, 4.200);

	auto cys = Residue<double>(atomVector<double>({ N, CA, C, O, CB, SG }), "CYS", "CYS1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(SG->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Glutamine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -20.534, 60.952, -3.420);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -19.881, 59.679, -3.725);
	auto C = std::make_shared<Atom<double>>("C", "C", -20.225, 59.151, -5.117);
	auto O = std::make_shared<Atom<double>>("O", "O", -19.341, 58.715, -5.857);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -20.227, 58.615, -2.682);
	auto CG = std::make_shared<Atom<double>>("C", "CG", -19.623, 58.862, -1.314);
	auto CD = std::make_shared<Atom<double>>("C", "CD", -19.831, 57.689, -0.378);
	auto OE1 = std::make_shared<Atom<double>>("O", "OE1", -19.110, 56.693, -0.443);
	auto NE2 = std::make_shared<Atom<double>>("N", "NE2", -20.832, 57.795, 0.490);

	auto gln = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD, OE1, NE2 }), "GLN", "GLN1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(CD->getNumberOfBonds(), 3);
	ASSERT_EQ(OE1->getNumberOfBonds(), 1);
	ASSERT_EQ(NE2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Glutamate)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -22.195, 56.100, 5.1976);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -22.549, 55.061, 4.2346);
	auto C = std::make_shared<Atom<double>>("C", "C", -21.673, 55.006, 2.9745);
	auto O = std::make_shared<Atom<double>>("O", "O", -22.031, 55.549, 1.9244);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -24.020, 55.234, 3.8415);
	auto CG = std::make_shared<Atom<double>>("C", "CG", -24.549, 54.201, 2.8649);
	auto CD = std::make_shared<Atom<double>>("C", "CD", -26.053, 54.305, 2.6729);
	auto OE1 = std::make_shared<Atom<double>>("O", "OE1", -26.556, 55.433, 2.4655);
	auto OE2 = std::make_shared<Atom<double>>("O", "OE2", -26.732, 53.259, 2.7241);

	auto glu = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD, OE1, OE2 }), "GLU", "GLU1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(CD->getNumberOfBonds(), 3);
	ASSERT_EQ(OE1->getNumberOfBonds(), 1);
	ASSERT_EQ(OE2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Glycine)
{

	auto CA = std::make_shared<Atom<double>>("C", "CA", 24.891, 31.124, -1.374);
	auto C = std::make_shared<Atom<double>>("C", "C", 25.888, 31.452, -2.472);
	auto O = std::make_shared<Atom<double>>("O", "O", 25.658, 32.337, -3.298);
	auto N = std::make_shared<Atom<double>>("N", "N", 24.005, 32.219, -1.031);

	auto gly = Residue<double>(atomVector<double>({ CA, C, O, N }), "GLY", "GLY1");

	ASSERT_EQ(gly.getBackbone()[0]->getX(), 24.005);
	ASSERT_EQ(gly.getBackbone()[0]->getY(), 32.219);
	ASSERT_EQ(gly.getBackbone()[0]->getZ(), -1.031);
	ASSERT_EQ(CA->getNumberOfBonds(), 2);
	ASSERT_EQ(C->getNumberOfBonds(), 2);
	ASSERT_EQ(gly.getBackbone()[0].get(), N.get());

	CA->destroyBond(1);
	ASSERT_EQ(CA->getNumberOfBonds(), 1);
	ASSERT_EQ(C->getNumberOfBonds(), 1);

	gly.createBonds();
	ASSERT_EQ(CA->getNumberOfBonds(), 2);
	ASSERT_EQ(C->getNumberOfBonds(), 2);

	ASSERT_EQ(C->getRadius(), 1.76);
	ASSERT_EQ(CA->getRadius(), 1.87);
	ASSERT_EQ(N->getRadius(), 1.65);
	ASSERT_EQ(O->getRadius(), 1.40);
}

TEST(ResidueTest, Histidine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 8.641, 40.445, 10.288);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 8.547, 41.617, 11.155);
	auto C = std::make_shared<Atom<double>>("C", "C", 9.195, 41.422, 12.521);
	auto O = std::make_shared<Atom<double>>("O", "O", 8.833, 40.518, 13.278);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 7.083, 42.025, 11.325);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 6.904, 43.335, 12.023);
	auto ND1 = std::make_shared<Atom<double>>("N", "ND1", 6.919, 43.456, 13.396);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD2", 6.730, 44.587, 11.536);
	auto CE1 = std::make_shared<Atom<double>>("C", "CE1", 6.759, 44.726, 13.726);
	auto NE2 = std::make_shared<Atom<double>>("N", "NE2", 6.644, 45.433, 12.616);

	auto his = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2 }), "HIS", "HIS1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(ND1->getNumberOfBonds(), 2);
	ASSERT_EQ(CD2->getNumberOfBonds(), 2);
	ASSERT_EQ(CE1->getNumberOfBonds(), 2);
	ASSERT_EQ(NE2->getNumberOfBonds(), 2);
}

TEST(ResidueTest, Isoleucine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 20.710, 62.644, 1.319);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 19.376, 62.572, 1.899);
	auto C = std::make_shared<Atom<double>>("C", "C", 19.455, 62.782, 3.403);
	auto O = std::make_shared<Atom<double>>("O", "O", 19.954, 61.929, 4.130);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 18.712, 61.210, 1.629);
	auto CG1 = std::make_shared<Atom<double>>("C", "CG1", 18.583, 60.982, 0.118);
	auto CG2 = std::make_shared<Atom<double>>("C", "CG2", 17.339, 61.159, 2.301);
	auto CD1 = std::make_shared<Atom<double>>("C", "CD1", 17.929, 59.660, -0.258);

	auto ile = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG1, CG2, CD1 }), "ILE", "ILE1");

	ASSERT_EQ(CB->getNumberOfBonds(), 3);
	ASSERT_EQ(CG1->getNumberOfBonds(), 2);
	ASSERT_EQ(CG2->getNumberOfBonds(), 1);
	ASSERT_EQ(CD1->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Leucine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 20.621, 54.814, -13.141);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 20.942, 54.764, -11.719);
	auto C = std::make_shared<Atom<double>>("C", "C", 22.406, 54.421, -11.450);
	auto O = std::make_shared<Atom<double>>("O", "O", 23.004, 54.929, -10.502);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 20.079, 53.714, -11.010);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 18.602, 54.001, -10.736);
	auto CD1 = std::make_shared<Atom<double>>("C", "CD1", 17.987, 52.800, -10.026);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD2", 18.467, 55.250, -9.879);

	auto ile = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD1, CD2 }), "LEU", "LEU1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(CD1->getNumberOfBonds(), 1);
	ASSERT_EQ(CD2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Lysine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 5.800, 40.909, -14.865);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 5.579, 42.240, -15.418);
	auto C = std::make_shared<Atom<double>>("C", "C", 6.695, 43.231, -15.104);
	auto O = std::make_shared<Atom<double>>("O", "O", 7.513, 43.002, -14.215);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 4.238, 42.791, -14.925);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 3.038, 41.989, -15.414);
	auto CD = std::make_shared<Atom<double>>("C", "CD", 1.744, 42.750, -15.212);
	auto CE = std::make_shared<Atom<double>>("C", "CE", 0.562, 41.997, -15.804);
	auto NZ = std::make_shared<Atom<double>>("N", "NZ", -0.711, 42.759, -15.643);

	auto lys = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD, CE, NZ }), "LYS", "LYS1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(CD->getNumberOfBonds(), 2);
	ASSERT_EQ(CE->getNumberOfBonds(), 2);
	ASSERT_EQ(NZ->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Methionine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 1.706, 59.788, 0.807);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 0.967, 60.974, 0.418);
	auto C = std::make_shared<Atom<double>>("C", "C", 1.840, 62.205, 0.579);
	auto O = std::make_shared<Atom<double>>("O", "O", 2.647, 62.287, 1.507);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -0.276, 61.149, 1.297);
	auto CG = std::make_shared<Atom<double>>("C", "CG", -1.147, 62.346, 0.906);
	auto SD = std::make_shared<Atom<double>>("S", "SD", -2.207, 62.955, 2.255);
	auto CE = std::make_shared<Atom<double>>("C", "CE", -1.060, 64.135, 3.031);

	auto met = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, SD, CE }), "MET", "MET1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(SD->getNumberOfBonds(), 2);
	ASSERT_EQ(CE->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Phenylalanine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 3.325, 74.582, -1.218);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 3.355, 74.056, -2.568);
	auto C = std::make_shared<Atom<double>>("C", "C", 4.743, 73.569, -2.957);
	auto O = std::make_shared<Atom<double>>("O", "O", 5.750, 74.161, -2.575);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 2.904, 75.136, -3.553);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 2.939, 74.698, -4.986);
	auto CD1 = std::make_shared<Atom<double>>("C", "CD1", 2.057, 73.726, -5.452);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD2", 3.857, 75.255, -5.874);
	auto CE1 = std::make_shared<Atom<double>>("C", "CE1", 2.086, 73.315, -6.781);
	auto CE2 = std::make_shared<Atom<double>>("C", "CE2", 3.896, 74.854, -7.206);
	auto CZ = std::make_shared<Atom<double>>("C", "CZ", 3.008, 73.880, -7.663);

	auto phe = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ }), "PHE", "PHE1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(CD1->getNumberOfBonds(), 2);
	ASSERT_EQ(CD2->getNumberOfBonds(), 2);
	ASSERT_EQ(CE1->getNumberOfBonds(), 2);
	ASSERT_EQ(CE2->getNumberOfBonds(), 2);
	ASSERT_EQ(CZ->getNumberOfBonds(), 2);
}

TEST(ResidueTest, Proline)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 29.627, 36.149, -5.045);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 28.814, 37.262, -5.534);
	auto C = std::make_shared<Atom<double>>("C", "C", 29.085, 37.599, -6.991);
	auto O = std::make_shared<Atom<double>>("O", "O", 29.997, 37.049, -7.614);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 29.236, 38.403, -4.625);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 30.713, 38.151, -4.500);
	auto CD = std::make_shared<Atom<double>>("C", "CD", 30.764, 36.646, -4.248);

	auto pro = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD }), "PRO", "PRO1");

	ASSERT_EQ(N->getNumberOfBonds(), 2);
	ASSERT_EQ(CA->getNumberOfBonds(), 3);
	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 2);
	ASSERT_EQ(CD->getNumberOfBonds(), 2);
}

TEST(ResidueTest, Serine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -5.690, 69.889, -2.546);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -6.410, 71.066, -2.086);
	auto C = std::make_shared<Atom<double>>("C", "C", -6.672, 71.030, -0.596);
	auto O = std::make_shared<Atom<double>>("O", "O", -6.806, 69.963, 0.005);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -7.741, 71.208, -2.830);
	auto OG = std::make_shared<Atom<double>>("O", "OG", -8.654, 70.206, -2.429);

	auto ser = Residue<double>(atomVector<double>({ N, CA, C, O, CB, OG }), "SER", "SER1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(OG->getNumberOfBonds(), 1);
}

TEST(ResidueTest, TerminalThreonine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 25.016, 52.737, 10.944);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 26.274, 52.269, 10.375);
	auto C = std::make_shared<Atom<double>>("C", "C", 27.057, 53.476, 9.871);
	auto O = std::make_shared<Atom<double>>("O", "O", 26.424, 54.533, 9.656);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 26.025, 51.306, 9.184);
	auto OG1 = std::make_shared<Atom<double>>("O", "OG1", 25.343, 50.133, 9.648);
	auto CG2 = std::make_shared<Atom<double>>("O", "CG2", 27.343, 50.897, 8.535);
	auto OXT = std::make_shared<Atom<double>>("O", "OXT", 28.285, 53.345, 9.687);

	auto thr = Residue<double>(atomVector<double>({ N, CA, C, O, CB, OG1, CG2, OXT }), "THR", "THR1");

	ASSERT_EQ(CB->getNumberOfBonds(), 3);
	ASSERT_EQ(OG1->getNumberOfBonds(), 1);
	ASSERT_EQ(CG2->getNumberOfBonds(), 1);
}

TEST(ResidueTest, Tryptophan)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 23.018, 53.170, 13.888);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 23.537, 53.665, 12.614);
	auto C = std::make_shared<Atom<double>>("C", "C", 24.927, 53.203, 12.188);
	auto O = std::make_shared<Atom<double>>("O", "O", 25.887, 53.273, 12.953);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 23.482, 55.194, 12.602);
	auto CG = std::make_shared<Atom<double>>("C", "CG", 22.093, 55.717, 12.821);
	auto CD1 = std::make_shared<Atom<double>>("C", "CD1", 21.623, 56.359, 13.929);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD2", 20.982, 55.592, 11.927);
	auto NE1 = std::make_shared<Atom<double>>("N", "NE1", 20.286, 56.642, 13.783);
	auto CE2 = std::make_shared<Atom<double>>("C", "CE2", 19.866, 56.182, 12.563);
	auto CE3 = std::make_shared<Atom<double>>("C", "CE3", 20.819, 55.037, 10.649);
	auto CZ2 = std::make_shared<Atom<double>>("C", "CZ2", 18.602, 56.233, 11.964);
	auto CZ3 = std::make_shared<Atom<double>>("C", "CZ3", 19.561, 55.088, 10.054);
	auto CH2 = std::make_shared<Atom<double>>("C", "CH2", 18.469, 55.682, 10.713);

	auto trp = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3,
								   CH2 }),
		"TRP", "TRP1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(CD1->getNumberOfBonds(), 2);
	ASSERT_EQ(CD2->getNumberOfBonds(), 3);
	ASSERT_EQ(NE1->getNumberOfBonds(), 2);
	ASSERT_EQ(CE2->getNumberOfBonds(), 3);
	ASSERT_EQ(CE3->getNumberOfBonds(), 2);
	ASSERT_EQ(CZ2->getNumberOfBonds(), 2);
	ASSERT_EQ(CZ3->getNumberOfBonds(), 2);
	ASSERT_EQ(CH2->getNumberOfBonds(), 2);
}

TEST(ResidueTest, Tyrosine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", -6.507, 82.172, -8.076);
	auto CA = std::make_shared<Atom<double>>("C", "CA", -6.693, 82.879, -6.809);
	auto C = std::make_shared<Atom<double>>("C", "C", -6.471, 81.928, -5.627);
	auto O = std::make_shared<Atom<double>>("O", "O", -5.663, 82.201, -4.741);
	auto CB = std::make_shared<Atom<double>>("C", "CB", -8.104, 83.477, -6.723);
	auto CG = std::make_shared<Atom<double>>("C", "CG", -8.302, 84.388, -5.528);
	auto CD1 = std::make_shared<Atom<double>>("C", "CD1", -7.742, 85.668, -5.500);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD2", -9.028, 83.962, -4.415);
	auto CE1 = std::make_shared<Atom<double>>("C", "CE1", -7.902, 86.501, -4.387);
	auto CE2 = std::make_shared<Atom<double>>("C", "CE2", -9.193, 84.785, -3.302);
	auto CZ = std::make_shared<Atom<double>>("C", "CZ", -8.629, 86.048, -3.293);
	auto OH = std::make_shared<Atom<double>>("O", "OH", -8.783, 86.853, -2.187);

	auto tyr = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH }),
		"TYR", "TYR1");

	ASSERT_EQ(CB->getNumberOfBonds(), 2);
	ASSERT_EQ(CG->getNumberOfBonds(), 3);
	ASSERT_EQ(CD1->getNumberOfBonds(), 2);
	ASSERT_EQ(CD2->getNumberOfBonds(), 2);
	ASSERT_EQ(CE1->getNumberOfBonds(), 2);
	ASSERT_EQ(CE2->getNumberOfBonds(), 2);
	ASSERT_EQ(CZ->getNumberOfBonds(), 2);
	ASSERT_EQ(OH->getNumberOfBonds(), 2);
}

TEST(ResidueTest, Valine)
{

	auto N = std::make_shared<Atom<double>>("N", "N", 0.691, 78.404, 3.398);
	auto CA = std::make_shared<Atom<double>>("C", "CA", 0.935, 78.627, 4.815);
	auto C = std::make_shared<Atom<double>>("C", "C", 2.047, 79.660, 4.958);
	auto O = std::make_shared<Atom<double>>("O", "O", 3.073, 79.572, 4.281);
	auto CB = std::make_shared<Atom<double>>("C", "CB", 1.361, 77.315, 5.521);
	auto CG1 = std::make_shared<Atom<double>>("C", "CG1", 1.989, 77.618, 6.885);
	auto CG2 = std::make_shared<Atom<double>>("C", "CG2", 0.147, 76.407, 5.687);

	auto val = Residue<double>(atomVector<double>({ N, CA, C, O, CB, CG1, CG2 }), "VAL", "VAL1");

	ASSERT_EQ(CB->getNumberOfBonds(), 3);
	ASSERT_EQ(CG1->getNumberOfBonds(), 1);
	ASSERT_EQ(CG2->getNumberOfBonds(), 1);
}
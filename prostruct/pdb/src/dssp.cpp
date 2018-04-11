//
// Created by gil on 06/04/18.
//
//
//The DSSP algorithm is the standard method for assigning secondary structure to the amino acids of a protein, given the
// atomic-resolution coordinates of the protein. The abbreviation is only mentioned once in the 1983 paper describing
// this algorithm,[1] where it is the name of the Pascal program that implements the algorithm Define Secondary Structure of Proteins.
//
//DSSP begins by identifying the intra-backbone hydrogen bonds of the protein using a purely electrostatic definition,
// assuming partial charges of -0.42 e and +0.20 e to the carbonyl oxygen and amide hydrogen respectively, their opposites
// assigned to the carbonyl carbon and amide nitrogen. A hydrogen bond is identified if E in the following equation
// is less than -0.5 kcal/mol:
//
//        E = 0.084 { 1 / rON + 1 / rCH − 1 / rOH − 1 / rCN } ⋅ 332 kcal/mol
//
//where the r A B {\displaystyle r_{AB}} r_{{AB}} terms indicate the distance between atoms A and B, taken from the
// carbon (C) and oxygen (O) atoms of the C=O group and the nitrogen (N) and hydrogen (H) atoms of the N-H group.
//

#include "../include/geometry.h"

#define ARMA_NO_DEBUG

enum SS_Types {
    Helix_310,
    Helix_alpha,
    Helix_pi,
    Bridge_beta,
    Bridge_beta_buldge,
    Loop_turn,
    Loop_high_curvature,
    Blank
};

static void predict_H_coords(arma::mat &H_coords, const arma::mat &C_coords, const arma::mat &O_coords,
                      const arma::mat &N_coords) {

//    for (arma::uword j = 1; j < H_coords.n_cols; ++j) {
//        H_coords.at(0, j - 1) = j;
//    }
//
//    H_coords.each_col([&](arma::vec &col) {col = (C_coords.col(static_cast<const arma::uword>(col[0])-1) -
//                                                  O_coords.col(static_cast<const arma::uword>(col[0])-1)) /
//                                                  arma::norm(C_coords.col(static_cast<const arma::uword>(col[0])-1) -
//                                                             O_coords.col(static_cast<const arma::uword>(col[0])-1), 2) +
//                                                  N_coords.col(static_cast<const arma::uword>(col[0]));});

    arma::vec co(3);

    for (arma::uword i = 1; i < H_coords.n_cols; ++i) {

        co.at(0) = C_coords.at(0, i - 1) - O_coords.at(0, i - 1);
        co.at(1) = C_coords.at(1, i - 1) - O_coords.at(1, i - 1);
        co.at(2) = C_coords.at(2, i - 1) - O_coords.at(2, i - 1);

        // CO norm -> distance between C and O in backbone of residue i - 1
        double co_norm = arma::norm(co, 2);

        co.at(0) /= co_norm;
        co.at(1) /= co_norm;
        co.at(2) /= co_norm;

        // calculate coordinates of H bound to N in backbone of residue i
        H_coords.at(0, i) = co.at(0) + N_coords.at(0, i);
        H_coords.at(1, i) = co.at(1) + N_coords.at(1, i);
        H_coords.at(2, i) = co.at(2) + N_coords.at(2, i);

    }
}


void kabsch_sander(const arma::mat &C_coords, const arma::mat &O_coords, const arma::mat &N_coords, const arma::mat &CA_coords,
                   std::vector<bool> &hasHbond, arma::mat &E, const arma::uword n_residues) {

    double ca_dist_squared = 81;
    static arma::vec zerosVec = std::vector<double>({0,0,0});

    arma::mat H_coords(3, n_residues);
//    H_coords.insert_cols(0, zerosVec);

    predict_H_coords(H_coords, C_coords, O_coords, N_coords);
#pragma omp parallel for collapse(2)
    for (arma::uword acceptor = 0; acceptor < n_residues; ++acceptor) {
        for (arma::uword donor = 0; donor < n_residues; ++donor) {

            if (std::abs(static_cast<int>(acceptor - donor)) != 1 && acceptor != donor) {

                if (arma::dot(CA_coords.col(donor) - CA_coords.col(acceptor),
                              CA_coords.col(donor) - CA_coords.col(acceptor)) < ca_dist_squared) {
                    // E = 0.084 { 1 / rON + 1 / rCH − 1 / rOH − 1 / rCN } ⋅ 332 kcal/mol
                    // where r is the distance between A and B sqrt(dot(A-B, A-B)
                    // and we do this for each possible combination -> gives a matrix residue x residue
                    E.at(acceptor, donor) = (1 / arma::norm(N_coords.col(donor) - O_coords.col(acceptor), 2) +
                                             1 / arma::norm(H_coords.col(donor) - C_coords.col(acceptor), 2) -
                                             1 / arma::norm(H_coords.col(donor) - O_coords.col(acceptor), 2) -
                                             1 / arma::norm(N_coords.col(donor) - C_coords.col(acceptor), 2)) * 27.88;
                }

                if (E.at(acceptor, donor) < -0.5) {
                    hasHbond[acceptor] = true;
                    hasHbond[donor] = true;
                }
            }
        }
    }
}


static void predict_alpha_helix() {

    //Based on this, eight types of secondary structure are assigned. The 310 helix, α helix and π helix have symbols G,
    // H and I and are recognized by having a repetitive sequence of hydrogen bonds in which the residues are three, four,
    // or five residues apart respectively. Two types of beta sheet structures exist; a beta bridge has symbol B while
    // longer sets of hydrogen bonds and beta bulges have symbol E. T is used for turns, featuring hydrogen bonds typical
    // of helices, S is used for regions of high curvature (where the angle between
    // Ciα C(i+2)α and C(i−2)α Ciα is at least 70°), and a blank (or space) is used if no other rule applies,
    // referring to loops.[2] These eight types are usually grouped into three larger classes: helix (G, H and I),
    // strand (E and B) and loop (S, T, and C, where C sometimes is represented also as blank space).




}

static void predict_beta_sheet() {



}


void dssp(const arma::mat &C_coords, const arma::mat &O_coords, const arma::mat &N_coords, const arma::mat &CA_coords) {

    arma::uword n_residues = C_coords.n_cols;
    std::vector<bool> has_Hbond(n_residues, false);

    // All the code is in column major -> each cartesian point is stored in a column (rather than a row)
    arma::mat E(n_residues, n_residues);

    kabsch_sander(C_coords, O_coords, N_coords, CA_coords, has_Hbond, E, n_residues);

    std::vector<SS_Types> secondaryStructure(n_residues);

    predict_alpha_helix();

    predict_beta_sheet();

}
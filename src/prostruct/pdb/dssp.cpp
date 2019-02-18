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

#include "prostruct/pdb/geometry.h"

namespace prostruct {
    namespace geometry {

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

        template<typename T>
        void predict_H_coords(arma::Mat<T> &H_coords, const arma::Mat<T> &C_coords, const arma::Mat<T> &O_coords,
                              const arma::Mat<T> &N_coords) {

            auto co = (C_coords(arma::span::all, arma::span(0, H_coords.n_cols - 2)) -
                       O_coords(arma::span::all, arma::span(0, H_coords.n_cols - 2))).eval();
            auto co_norm = arma::sqrt(arma::sum(arma::square(co), 0));
            co.each_row() /= co_norm;

            H_coords(arma::span::all, arma::span(1, H_coords.n_cols - 1)) =
                    co + N_coords(arma::span::all, arma::span(1, H_coords.n_cols - 1));
        }

        template<typename T>
        void kabsch_sander(const arma::Mat<T> &C_coords, const arma::Mat<T> &O_coords, const arma::Mat<T> &N_coords,
                           const arma::Mat<T> &CA_coords, arma::Mat<T> &E, arma::uword n_residues) {

            T ca_dist_squared = 81;
            static arma::Col<T> zerosVec = std::vector<T>({0, 0, 0});

            arma::Mat<T> H_coords(3, n_residues, arma::fill::zeros);
            //    H_coords.insert_cols(0, zerosVec);

            predict_H_coords(H_coords, C_coords, O_coords, N_coords);
#pragma omp parallel for collapse(2)
            for (arma::uword acceptor = 0; acceptor < n_residues; ++acceptor) {
                for (arma::uword donor = 0; donor < n_residues; ++donor) {

                    if (std::abs(static_cast<int>(acceptor - donor)) != 1 && acceptor != donor) {

                        if (arma::dot(CA_coords.col(donor) - CA_coords.col(acceptor),
                                      CA_coords.col(donor) - CA_coords.col(acceptor))
                            < ca_dist_squared) {
                            // E = 0.084 { 1 / rON + 1 / rCH − 1 / rOH − 1 / rCN } ⋅ 332 kcal/mol
                            // where r is the distance between A and B sqrt(dot(A-B, A-B)
                            // and we do this for each possible combination -> gives a matrix residue x residue
                            T rev_rON = 1 / std::sqrt(arma::dot(N_coords.col(donor) - O_coords.col(acceptor), N_coords.col(donor) - O_coords.col(acceptor)));
                            T rev_rCH = 1 / std::sqrt(arma::dot(H_coords.col(donor) - C_coords.col(acceptor), H_coords.col(donor) - C_coords.col(acceptor)));
                            T rev_rOH = 1 / std::sqrt(arma::dot(H_coords.col(donor) - O_coords.col(acceptor), H_coords.col(donor) - O_coords.col(acceptor)));
                            T rev_rCN = 1 / std::sqrt(arma::dot(N_coords.col(donor) - C_coords.col(acceptor), N_coords.col(donor) - C_coords.col(acceptor)));
                            E.at(acceptor, donor) = (rev_rON + rev_rCH - rev_rOH - rev_rCN) * 27.88;
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

//template <typename T>
        void predict_beta_sheet() {
        }

        template<typename T>
        void dssp(const arma::Mat<T> &C_coords, const arma::Mat<T> &O_coords, const arma::Mat<T> &N_coords,
                  const arma::Mat<T> &CA_coords) {

            arma::uword n_residues = C_coords.n_cols;
            std::vector<bool> has_Hbond(n_residues, false);

            // All the code is in column major -> each cartesian point is stored in a column (rather than a row)
            arma::Mat<T> E(n_residues, n_residues);

            kabsch_sander(C_coords, O_coords, N_coords, CA_coords, E, n_residues);

            std::vector<SS_Types> secondaryStructure(n_residues);

            predict_alpha_helix();

            predict_beta_sheet();
        }

        template void kabsch_sander(const arma::Mat<float> &, const arma::Mat<float> &, const arma::Mat<float> &,
                                    const arma::Mat<float> &,
                                    arma::Mat<float> &, const arma::uword);

        template void kabsch_sander(const arma::Mat<double> &, const arma::Mat<double> &, const arma::Mat<double> &,
                                    const arma::Mat<double> &,
                                    arma::Mat<double> &, const arma::uword);

        template void
        dssp(const arma::Mat<float> &, const arma::Mat<float> &, const arma::Mat<float> &, const arma::Mat<float> &);

        template void dssp(const arma::Mat<double> &, const arma::Mat<double> &, const arma::Mat<double> &,
                           const arma::Mat<double> &);
    }
}
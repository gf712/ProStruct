//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#include "chain.h"
#include "PDBparser.h"
#include "geometry.h"

class PDB {

public:

    PDB(std::string filename);

    ~PDB() = default;

    static PDB fetch(std::string);

    std::string getFilename() { return filename;}

    int n_chains() { return numberOfChains; }

    std::vector<std::string> getChainIDs() {
        return chainOrder;
    }

    std::shared_ptr<Chain> getChain(std::string name_) { return chainMap[name_]; }

    arma::mat calculate_KabschSander();

    arma::mat predict_backboneHbonds();

    void calculate_dssp();

    arma::mat getXYZ() { return xyz; }

    int n_residues() { return nResidues; }

    int n_atoms() { return nAtoms; }

    arma::vec getRadii() { return radii; }

    arma::vec calculate_ASA(double probe);

    double calculate_RMSD(PDB& other);


private:

    arma::mat xyz;
    std::string filename;
    int numberOfChains;
    std::map<std::string, std::shared_ptr<Chain>> chainMap;
    int nAtoms;
    std::vector<std::string> chainOrder;
    arma::uword nResidues;
    arma::vec radii;
    void getBackboneAtoms(arma::mat&, arma::mat&, arma::mat&, arma::mat&);
    void internalKS(arma::mat&);
};


#endif //PROSTRUCT_PDB_H

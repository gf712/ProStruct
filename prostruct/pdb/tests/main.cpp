//
// Created by gil on 29/03/18.
//

#include "PDB.h"
#include "geometry.h"
#include <armadillo>

int main() {

    std::string file = "test.pdb";


    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 100; ++i) {

        auto pdb = PDB(file);
    }

//    for (int i = 0; i < 1000; ++i) {
//        pdb.calculate_dssp();
//    }

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();

    printf ("Elapsed time is %li microseconds.", duration / 100);

}



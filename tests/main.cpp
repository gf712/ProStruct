//
// Created by gil on 29/03/18.
//

#include "prostruct/prostruct.h"
#include <armadillo>

using namespace prostruct;

int main() {

    std::string file = "test.pdb";

    auto pdb = PDB<double>(file);

//    auto start = std::chrono::high_resolution_clock::now();

//    for (int i = 0; i < 100; ++i) {
//
//        auto pdb = PDB(file);
//    }

//    for (int i = 0; i < 10000; ++i) {
        pdb.compute_kabsch_sander();
//    }

//    arma::Mat<T> result(3, 1000);
//    generate_sphere(1000, result);

//    std::cout << std::accumulate(pdb.predict_backboneHbonds().begin(), pdb.predict_backboneHbonds().end(), 0) << std::endl;

//    arma::Col<T> test;

//    for (int i = 0; i < 100; ++i) {

//        test = pdb.calculate_ASA(1.4);

//    }

//    pdb.getXYZ().col(1).print();
//    arma::Mat<T> xyz = pdb.getXYZ();
//    for (int j = 0; j < 10000; ++j) {
//
//        double sum = 0.0;

//    arma::Mat<T>::col_iterator col_it     = xyz.begin_col(0);  // start of column 1
//    arma::Mat<T>::col_iterator col_it_end = xyz.end_col(xyz.n_cols);    //   end of column 3


//    }


//    std::cout << sum;

//    std::cout << arma::norm(pdb.getXYZ(), 2);

//    test.print();

//    calculate_SASA(pdb.getXYZ(), pdb.n_atoms(), 2);

//    auto end = std::chrono::high_resolution_clock::now();
//
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
//
//    printf ("Elapsed time is %lli microseconds.\n", duration/1000);

}



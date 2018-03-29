//
// Created by gil on 29/03/18.
//

#include "residue.h"
#include <chrono>
#include <iostream>

void main_helper() {

    auto N = std::make_shared<Atom>("N", "N", 23.018, 53.170, 13.888);
    auto CA = std::make_shared<Atom>("C", "CA", 23.537, 53.665, 12.614);
    auto C = std::make_shared<Atom>("C", "C", 24.927, 53.203, 12.188);
    auto O = std::make_shared<Atom>("O", "O", 25.887, 53.273, 12.953);
    auto CB = std::make_shared<Atom>("C", "CB", 23.482, 55.194, 12.602);
    auto CG = std::make_shared<Atom>("C", "CG", 22.093, 55.717, 12.821);
    auto CD1 = std::make_shared<Atom>("C", "CD1", 21.623, 56.359, 13.929);
    auto CD2 = std::make_shared<Atom>("C", "CD2", 20.982, 55.592, 11.927);
    auto NE1 = std::make_shared<Atom>("N", "NE1", 20.286, 56.642, 13.783);
    auto CE2 = std::make_shared<Atom>("C", "CE2", 19.866, 56.182, 12.563);
    auto CE3 = std::make_shared<Atom>("C", "CE3", 20.819, 55.037, 10.649);
    auto CZ2 = std::make_shared<Atom>("C", "CZ2", 18.602, 56.233, 11.964);
    auto CZ3 = std::make_shared<Atom>("C", "CZ3", 19.561, 55.088, 10.054);
    auto CH2 = std::make_shared<Atom>("C", "CH2", 18.469, 55.682, 10.713);

    auto trp = Residue(std::vector<std::shared_ptr<Atom>>({N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3,
                                                           CH2}),
                       "TRP", "TRP1");
}

int main() {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int i = 0;
    while (i < 1000) {
        main_helper();
        i++;
    }

    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<std::endl;

    return 0;
}
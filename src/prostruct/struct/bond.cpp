//
// Created by gil on 27/03/18.
//

#include <utility>
#include <cmath>
#include <iostream>

#include "prostruct/struct/bond.h"
#include "prostruct/struct/atom.h"

template <typename T>
Bond<T>::Bond(std::shared_ptr<Atom<T>> Atom1, std::shared_ptr<Atom<T>> Atom2, int bondType_) {
    
    atom1 = std::weak_ptr<Atom<T>>(Atom1);
    atom2 = std::weak_ptr<Atom<T>>(Atom2);

    initialiseBond(bondType_);
}

template <typename T>
Bond<T>::Bond(T x, T y, T z, int bondType) {

    auto atom1_ = std::make_shared<Atom<T>>("H", "Atom1", 0., 0., 0.);
    auto atom2_ = std::make_shared<Atom<T>>("H", "Atom2", x, y, z);

    atom1 = atom1_;
    atom2 = atom2_;

    initialiseBond(bondType);
}

template <typename T>
Bond<T>::Bond(T x1, T y1, T z1, T x2, T y2, T z2, int bondType) {

    atom1 = std::make_shared<Atom<T>>("H", "Atom1", x1, y1, z1);
    atom2 = std::make_shared<Atom<T>>("H", "Atom2", x2, y2, z2);

    initialiseBond(bondType);
}

template <typename T>
void Bond<T>::initialiseBond(int bondType_) {

    // x, y, z is the vector corresponding to the bond
    x = atom2.lock()->getX() - atom1.lock()->getX();
    y = atom2.lock()->getY() - atom1.lock()->getY();
    z = atom2.lock()->getZ() - atom1.lock()->getZ();

    bondVector = {x, y, z};

    bondType = bondType_;

    // calculate magnitude/length
    length = std::sqrt(x*x + y*y + z*z);
}

template class Bond<float>;
template class Bond<double>;
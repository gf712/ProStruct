//
// Created by gil on 27/03/18.
//

#include <utility>
#include <iostream>

#include "../include/bond.h"
#include "../include/atom.h"

Bond::Bond(std::shared_ptr<Atom> Atom1, std::shared_ptr<Atom> Atom2, int bondType_) {

    std::cout << "Forming bond between: " << Atom1->getName() << " and " << Atom2->getName() << std::endl;

    atom1 = std::weak_ptr<Atom>(Atom1);
    atom2 = std::weak_ptr<Atom>(Atom2);

    initialiseBond(bondType_);
}

//Bond::Bond(double x, double y, double z, int bondType) {
//
//    atom1 = std::weak_ptr<Atom>("H", "Atom1", 0, 0, 0);
//    atom2 = std::weak_ptr<Atom>("H", "Atom2", x, y, z);
//
//    initialiseBond(bondType);
//}
//
//Bond::Bond(double x1, double y1, double z1, double x2, double y2, double z2, int bondType) {
//
//    atom1 = std::make_shared<Atom>("H", "Atom1", x1, y1, z1);
//    atom2 = std::make_shared<Atom>("H", "Atom2", x2, y2, z2);
//
//    initialiseBond(bondType);
//}

void Bond::initialiseBond(int bondType_) {

    x = atom2.lock()->getX() - atom1.lock()->getX();
    y = atom2.lock()->getY() - atom1.lock()->getY();
    z = atom2.lock()->getZ() - atom1.lock()->getZ();

    bondVector = {x, y, z};

    bondType = bondType_;

}
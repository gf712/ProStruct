//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_BOND_H
#define PROSTRUCT_BOND_H

#include <vector>
#include <memory>

class Atom;

class Bond {

public:

    Bond(double x, double y, double z, int bondType);
    Bond(double x1, double y1, double z1, double x2, double y2, double z2, int bondType);
    Bond(std::shared_ptr<Atom> Atom1, std::shared_ptr<Atom> Atom2, int bondType);

    void initialiseBond(int);
    std::vector<double> getBondVector() { return bondVector; }
    double getX() { return x; }
    double getY() { return y; }
    double getZ() { return z; }
    int getBondType() { return bondType; }

    std::shared_ptr<Atom> getAtom1() { return atom1.lock(); }
    std::shared_ptr<Atom> getAtom2() { return atom2.lock(); }

private:
    // keep a weak reference to std::shared<Atom>
    // otherwise Bond owns Atom, which doesn't make sense
    // and can cause memory issues
    std::weak_ptr<Atom> atom1, atom2;
    double x, y, z;
    std::vector<double> bondVector;
    int bondType;
    double length;
};


#endif //PROSTRUCT_BOND_H

//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_ATOM_H
#define PROSTRUCT_ATOM_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include "bond.h"
#include <armadillo>


class Atom: public std::enable_shared_from_this< Atom > {

public:
    Atom(std::string element) {
        load_atom(element);
    }
    Atom(std::string element, std::string name) {
        load_atom(element, name);
    }
    Atom(std::string element, std::string name, double x, double y, double z) {
        load_atom(element, name, x, y, z);
    };

//
//    Atom(std::string element, std::string name, double x, double y, double z) {
//        load_atom(element, name, x, y, z);
//    };

    void addBond(std::shared_ptr<Atom> atom, int bondType);
    void addBond(std::shared_ptr<Bond> bond);
    void destroyBond(int);
    void destroyBond(std::shared_ptr<Bond>);

    void setRadius(double radius_) {radius=radius_;}

    std::vector<std::shared_ptr<Bond>> getBonds() { return bonds; };
    int getNumberOfBonds() { return bonds.size(); }

    bool hasBond(const std::shared_ptr<Atom>&);

    double getAtomicWeight() { return atomicWeight; }
    int getAtomicNumber() { return atomicNumber; }

    double getX() { return x; }
    double getY() { return y; }
    double getZ() { return z; }
    double getRadius() { return radius; }
    std::string getElement() { return name; }
    arma::vec getXYZ() {return arma::vec(std::vector<double>({x, y, z}));}

    std::string getName() { return name; }

private:

    double x, y, z;
    double radius;
    void load_atom(std::string element);
    void load_atom(std::string element, std::string name);
    void load_atom(std::string element, std::string name, double x, double y, double z);
    double atomicWeight;
    int atomicNumber;
    std::string element;
    std::string name;
    std::vector<std::shared_ptr<Bond>> bonds;

};


#endif //PROSTRUCT_ATOM_H

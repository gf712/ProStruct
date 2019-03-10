/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_ATOM_H
#define PROSTRUCT_ATOM_H

#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include <prostruct/struct/bond.h>
#include <armadillo>

template <typename T>
using atomVector = std::vector<std::shared_ptr<Atom<T>>>;

template <typename T>
class Atom : public std::enable_shared_from_this<Atom<T>>
{

public:
	Atom(const std::string& element) { load_atom(element); }
	Atom(const std::string&, const std::string&) { load_atom(element, name); }
	Atom(const std::string& element, const std::string& name, T x, T y, T z)
	{
		load_atom(element, name, x, y, z);
	};

	std::shared_ptr<Atom<T>> getAtom() { return this->shared_from_this(); };

	void addBond(std::shared_ptr<Atom<T>> atom, int bondType);
	void addBond(std::shared_ptr<Bond<T>> bond);
	void destroyBond(int);
	void destroyBond(std::shared_ptr<Bond<T>>);

	void setRadius(double radius_) { radius = radius_; }

	std::vector<std::shared_ptr<Bond<T>>> getBonds() { return bonds; };
	int getNumberOfBonds() { return bonds.size(); }

	bool hasBond(const std::shared_ptr<Atom>&);

	double getAtomicWeight() { return atomicWeight; }
	int getAtomicNumber() { return atomicNumber; }

	T getX() { return x; }
	T getY() { return y; }
	T getZ() { return z; }
	T getRadius() { return radius; }
	std::string getElement() { return name; }
	arma::Col<T> getXYZ() { return arma::Col<T>(std::vector<T>({ x, y, z })); }

	std::string get_name() const noexcept { return name; }

private:
	T x, y, z;
	T radius;
	void load_atom(const std::string& element);
	void load_atom(const std::string& element, const std::string& name);
	void load_atom(const std::string& element, const std::string& name, T x, T y, T z);
	T atomicWeight;
	int atomicNumber;
	std::string element;
	std::string name;
	std::vector<std::shared_ptr<Bond<T>>> bonds;
};

#endif // PROSTRUCT_ATOM_H

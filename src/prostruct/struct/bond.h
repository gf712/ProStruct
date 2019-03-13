/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_BOND_H
#define PROSTRUCT_BOND_H

#include <memory>
#include <vector>

namespace prostruct {
	template<typename>
	class Atom;

	template<typename T>
	class Bond {

	public:
		Bond(T x, T y, T z, int bondType);

		Bond(T x1, T y1, T z1, T x2, T y2, T z2, int bondType);

		Bond(std::shared_ptr<Atom<T>> Atom1, std::shared_ptr<Atom<T>> Atom2, int bondType);

		void initialiseBond(int);

		std::vector<T> getBondVector() { return bondVector; }

		T getX() { return x; }

		T getY() { return y; }

		T getZ() { return z; }

		int getBondType() { return bondType; }

		std::shared_ptr<Atom<T>> getAtom1() { return atom1.lock(); }

		std::shared_ptr<Atom<T>> getAtom2() { return atom2.lock(); }

	private:
		// keep a weak reference to std::shared<Atom>
		// otherwise Bond owns Atom, which doesn't make sense
		// and can cause memory issues
		std::weak_ptr<Atom<T>> atom1;
		std::weak_ptr<Atom<T>> atom2;
		T x, y, z;
		std::vector<T> bondVector;
		int bondType;
		T length;
	};
}

#endif // PROSTRUCT_BOND_H

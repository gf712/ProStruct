//
// Created by gil on 09/04/18.
//

#include "prostruct/pdb/geometry.h"

namespace prostruct {
	namespace geometry {
		template<typename T>
		inline T dihedrals(const arma::Col<T>& atom1, const arma::Col<T>& atom2, const arma::Col<T>& atom3, const arma::Col<T>& atom4, T coef) 
		{
			arma::Col<T> b1 = arma::normalise(atom1 - atom2);
			arma::Col<T> b2 = arma::normalise(atom2 - atom3);
			arma::Col<T> b3 = arma::normalise(atom3 - atom4);
			arma::Col<T> n1 = arma::cross(b1, b2);	
			arma::Col<T> n2 = arma::cross(b2, b3);
			return std::atan2(arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2)) * coef;
		}

		template float dihedrals(const arma::Col<float>&, const arma::Col<float>&, const arma::Col<float>&, const arma::Col<float>&, float);
		template double dihedrals(const arma::Col<double>&, const arma::Col<double>&, const arma::Col<double>&, const arma::Col<double>&, double);
	}
}
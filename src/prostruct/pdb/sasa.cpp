//
// Created by gil on 09/04/18.
//

#include "prostruct/pdb/geometry.h"
#include "prostruct/struct/atom.h"

constexpr double golden_angle = 2.399963229728653;

template <typename T>
void generate_sphere(int N, arma::Mat<T>& result)
{

	T offset = 2.0 / N;

	for (arma::uword i = 0; i < N; ++i) {

		T y = i * offset - 1.0 + (offset / 2.0);
		T r = std::sqrt(1.0 - y * y);

		T t = i * golden_angle;

		result.at(0, i) = r * std::cos(t); // x
		result.at(1, i) = y; // y
		result.at(2, i) = r * std::sin(t); // z
	}
}

template <typename T>
void calculate_atom_SASA(const arma::Mat<T>& xyz, const arma::Col<T>& radius, const arma::subview_col<T>& neighbours,
	const int current_atom_index, const T probe, const arma::Mat<T>& sphere_points,
	const T adjustment, arma::Col<T>& asa)
{

	arma::Col<T> atom_XYZ = xyz.col(current_atom_index);
	arma::uvec neighbourIndices = arma::find(neighbours, 0);

	double atomRadius = probe + radius.at(current_atom_index);
	int nNeighbours = neighbourIndices.size();
	int accessiblePoints = 0;

	int k = 0;
	for (int i = 0; i < 1000; ++i) {
		arma::Col<T> current_sphere_point(3);
		current_sphere_point.at(0) = sphere_points.at(0, i) * atomRadius + atom_XYZ.at(0);
		current_sphere_point.at(1) = sphere_points.at(1, i) * atomRadius + atom_XYZ.at(1);
		current_sphere_point.at(2) = sphere_points.at(2, i) * atomRadius + atom_XYZ.at(2);

		for (int j = k; j < nNeighbours + k; ++j) {

			int index = neighbourIndices.at(j % nNeighbours);
			double r = radius.at(index) + probe;
			double dist = arma::sum(arma::square(current_sphere_point - xyz.col(index)));

			if (dist < r * r) {
				k = j;
				goto BURRIED;
			}
		}

		accessiblePoints++;
	BURRIED:
		continue;
	}

	asa.at(current_atom_index) = adjustment * accessiblePoints * atomRadius * atomRadius;
}

template <typename T>
void get_neighbours(const arma::Mat<T>& xyz, arma::Mat<T>& neighbours, int n_atoms, const arma::Col<T>& radii)
{
#pragma omp parallel for
	for (arma::uword i = 0; i < n_atoms; ++i) {
		for (arma::uword j = 1; j < n_atoms; ++j) {
			arma::Col<T> dist(3);
			double cutoff = radii.at(i) + radii.at(j);
			double cutoff_2 = cutoff * cutoff;
			dist.at(0) = xyz.at(0, i) - xyz.at(0, j);
			dist.at(1) = xyz.at(1, i) - xyz.at(1, j);
			dist.at(2) = xyz.at(2, i) - xyz.at(2, j);
			if (arma::sum(arma::square(dist)) < cutoff_2) {
				neighbours.at(i, j) = 1;
				neighbours.at(j, i) = 1;
			}
		}
	}
}

template <typename T>
void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii, arma::Col<T>& asa, int n_atoms, T probe)
{

	arma::Mat<T> sphere_points(3, 1000);
	arma::Mat<T> neighbours(n_atoms, n_atoms, arma::fill::zeros);

	generate_sphere(1000, sphere_points);

	get_neighbours(xyz, neighbours, n_atoms, radii);

	constexpr T adjustment = 4.0 * M_PI / 1000;

#pragma omp parallel for
	for (int i = 0; i < n_atoms; ++i) {
		calculate_atom_SASA(xyz, radii, neighbours.col(i), i, probe, sphere_points, adjustment, asa);
	}
}

template void calculate_atom_SASA(const arma::Mat<float>&, const arma::Col<float>&, const arma::subview_col<float>&,
	const int, const float, const arma::Mat<float>&, const float, arma::Col<float>&);
template void calculate_atom_SASA(const arma::Mat<double>&, const arma::Col<double>&, const arma::subview_col<double>&,
	const int, const double, const arma::Mat<double>&, const double, arma::Col<double>&);

template void shrake_rupley(const arma::Mat<float>&, const arma::Col<float>&, arma::Col<float>&, int n_atoms, float probe);
template void shrake_rupley(const arma::Mat<double>&, const arma::Col<double>&, arma::Col<double>&, int n_atoms, double probe);

//
// Created by gil on 09/04/18.
//

#include "prostruct/pdb/geometry.h"
#include "prostruct/struct/atom.h"

// compile time sqrt from https://baptiste-wicht.com/posts/2014/07/compile-integer-square-roots-at-compile-time-in-cpp.html
static constexpr std::size_t ct_sqrt(double res, std::size_t l, std::size_t r){
    if(l == r){
        return r;
    } else {
        const auto mid = (r + l) / 2;

        if(mid * mid >= res){
            return ct_sqrt(res, l, mid);
        } else {
            return ct_sqrt(res, mid + 1, r);
        }
    }
}

static constexpr std::size_t ct_sqrt(double res){
    return ct_sqrt(res, 1, res);
}

constexpr double golden_angle = M_PI * (3 - ct_sqrt(5.0));

static void generate_sphere(int N, arma::mat& result) {

    double offset = 2.0 / N;

    for (arma::uword i = 0; i < N; ++i) {

        double y = i * offset - 1 + (offset / 2);
        double r = std::sqrt(1 - y*y);

        double t = i * golden_angle;

        result.at(0, i) = r * std::cos(t); // x
        result.at(1, i) = y; // y
        result.at(2, i) = r * std::sin(t); // z
    }
}


static void calculate_atom_SASA(const arma::mat& xyz, const arma::vec& radius, const arma::vec& neighbours,
                         const int current_atom_index, const double probe, const arma::mat &sphere_points,
                         const double adjustment, arma::vec& asa) {

    arma::vec atom_XYZ = xyz.col(current_atom_index);
    arma::uvec neighbourIndices = arma::find(neighbours, 0);

    double atomRadius = probe + radius.at(current_atom_index);
    int nNeighbours = neighbourIndices.size();
    int accessiblePoints=0;

    int k = 0;
    for (int i = 0; i < 1000; ++i) {
        arma::vec current_sphere_point(3);
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

static void get_neighbours(const arma::mat &xyz, arma::mat& neighbours, int n_atoms, const arma::vec &radii) {
#pragma omp parallel for
    for (arma::uword i = 0; i < n_atoms; ++i) {
        for (arma::uword j = 1; j < n_atoms; ++j) {
            arma::vec dist(3);
            double cutoff = radii.at(i) + radii.at(j);
            double cutoff_2 = cutoff * cutoff;
            dist.at(0) = xyz.at(0, i) - xyz.at(0, j);
            dist.at(1) = xyz.at(1, i) - xyz.at(1, j);
            dist.at(2) = xyz.at(2, i) - xyz.at(2, j);
            if (arma::sum(arma::square(dist)) < cutoff_2) {
                neighbours.at(i,j) = 1;
                neighbours.at(j,i) = 1;
            }
        }
    }
}

void shrake_rupley(const arma::mat &xyz, const arma::vec &radii, arma::vec &asa, int n_atoms, double probe) {

    arma::mat sphere_points(3, 1000);
    arma::mat neighbours(n_atoms, n_atoms, arma::fill::zeros);

    generate_sphere(1000, sphere_points);

    get_neighbours(xyz, neighbours, n_atoms, radii);

    constexpr double adjustment = 4.0 * M_PI / 1000;

#pragma omp parallel for
    for (int i = 0; i < n_atoms; ++i) {
        calculate_atom_SASA(xyz, radii, neighbours.col(i), i, probe, sphere_points, adjustment, asa);
    }

//    asa.t().print();
}
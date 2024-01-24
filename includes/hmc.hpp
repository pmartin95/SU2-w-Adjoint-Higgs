#include "global_decl.hpp"
#include <Eigen/Dense>
#include "utility.hpp"
#include "lattice_ops.hpp"
#include "action.hpp"
#include <iostream>
#include "observables.hpp"
#pragma once
// update lattice momenta of lattice 1
void updateLatticeMomenta(site **current_momenta, site **current_lattice, site **next_momenta, double step_size);
// update lattice fields of lattice 1
void updateLatticeFields(site **current_lattice, site **current_momenta, site **next_lattice, double step_size);

// perform an HMC with t time intervals
void HMC_warmup(int t);
matrix averageMomenta(site *plattice1);
double HMC(int t);

// Copy lattice 1 to lattice 2
void copyLattice(site *lattice1, site *lattice2);

// Link field force term i.e. \partial S / \partial A_{mu}
void linkFieldForce(site *lattice1, int site_index, int mu, matrix &force);
const matrix linkFieldForce(site *lattice1, int site_index, int mu);

// gauge field staple
void gaugeFieldStaple(site *lattice1, int site_index, int mu, matrix &staple);
const matrix gaugeFieldStaple(site *lattice1, int site_index, int mu);
void gaugeFieldUpperStaple(site *lattice1, int site_index, int mu, matrix &Ustaple);
void gaugeFieldLowerStaple(site *lattice1, int site_index, int mu, matrix &Lstaple);
// Cayley-Klein version of exponential function
void expCK(const matrix &m, matrix &expm);
const matrix expCK(const matrix &m);
double totalMomentum(site *plattice1);
// Hamiltonian of the system
double hamiltonian(site *lattice1, site *plattice1);

// generate random normally distributed Hermitian matrices
void randomHermitianMatrix(matrix &m);

// file a lattice with random normally distirbuted Hermitian matrices
void randomMomentumLattice(site *lattice1);
matrix antiHermitianTraceless(matrix &m);
matrix hermitianTraceless(matrix &m);
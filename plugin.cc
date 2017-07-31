// Author:Tarini Hardikar
// Summer 2017
/*
 * @BEGIN LICENSE
 *
 * myscf by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/* This is entirely built from the tutorial to develop a SCF sovler by the
 * DePrince Group using the Psi4 conda.  http://psicode.org/
 * https://www.chem.fsu.edu/~deprince/programming_projects/scf/index.php
 */



#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

namespace psi{ namespace myscf {

extern "C"
int read_options(std::string name, Options& options)
{
	if (name == "MYSCF"|| options.read_globals()) {
		/*- The amount of information printed to the output file -*/
		options.add_int("PRINT", 1);
	}

	return true;
}

extern "C"
SharedWavefunction myscf(SharedWavefunction ref_wfn, Options& options)
{
	int print = options.get_int("PRINT");
	outfile->Printf("Hello from my first Psi4 plugin! THIS IS TERRIFYING.\n");

	/*XXX STEP ONE: ONE ELECTRON INTEGRALS XXX*/

	//grab the one-electron integrals from MintsHelper:
	std::shared_ptr<MintsHelper> mints (new MintsHelper(ref_wfn));

	/*the prefix so is for "symmetry orbital." The SO basis is obtained
	  by transforming the AO basis so that the SOs are
	irreducible reps in the point group*/

	//one-electron kinetic energry integrals
	std::shared_ptr<Matrix> T = mints->so_kinetic();

	//one-electron potential energy integrals
	std::shared_ptr<Matrix> V = mints->so_potential();

	//overlap intergrals
	std::shared_ptr<Matrix> S = mints->so_overlap();

	//build the core Hamiltonian (kinetic + potential)
	std::shared_ptr<Matrix> h = (std::shared_ptr<Matrix>)(new Matrix(T));
	h->add(V);


	/*XXX STEP TWO: TWO ELECTRON INTEGRALS XXX*/

	//However, we will be using a three-index integral with density fitting
	//instead of the four-index electron repulsion

	//grab the molecule from the wavefunction that was passed
	std::shared_ptr<Molecule> mol = ref_wfn->molecule();

	//get the primary basis set
	std::shared_ptr<BasisSet> primary = ref_wfn->get_basisset("ORBITAL");

	//the number of basis functions
	int nso = primary->nbf();

	//get auxillary basis
	std::shared_ptr<BasisSet> auxillary = ref_wfn->get_basisset("DF_BASIS_SCF");

	//the number of auxillary functions
	int nQ = auxillary->nbf();

	//using the molecule to calculate the total number of electrons
	int charge = mol->molecular_charge();
	int nelectron = 0;
	for (int i = 0; i < mol->natom(); i++){
	  nelectron += (int)mol->Z(i);
	}
	nelectron -= charge;

	//this is only for closed shells
	if (nelectron % 2 != 0){
		throw PsiException ("plugin myscf only works for closed shell", __FILE__, __LINE__);
	}

	//the number of doubly occupied orbitas
	int na = nelectron/2;

	double e_convergence = options.get_double("E_CONVERGENCE");
	double d_convergence = options.get_double("D_CONVERGENCE");
	int maxiter = options.get_int("MAXITER");

	outfile->Printf("\n");
	outfile->Printf("    No. basis functions:            %5i\n",nso);
	outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
	outfile->Printf("    No. electrons:                  %5i\n",nelectron);
	outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
	outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
	outfile->Printf("    maxiter:                        %5i\n",maxiter);
	outfile->Printf("\n");
	outfile->Printf("\n");

	//Roothan coefficients matrix
	std::shared_ptr<Matrix> Ca = (std::shared_ptr<Matrix>)(new Matrix(nso, nso));
	//Fock matrix
	std::shared_ptr<Matrix> Fa = (std::shared_ptr<Matrix>)(new Matrix(nso, nso));
	//Density matrix
	std::shared_ptr<Matrix> Da = (std::shared_ptr<Matrix>)(new Matrix(nso, nso));

	//construct the three-index integral
	std::shared_ptr<DFTensor> DF (new DFTensor(primary, auxillary, Ca, na, nso-na, na, nso-na, options));
	//Qso contations the three-index integral now
	std::shared_ptr<Matrix> Qso = DF->Qso();

	/*XXX STEP THREE: LOWDIN SYMMETRIC ORTHOGONALIZATION XXX*/

	//http://booksite.elsevier.com/9780444594365/downloads/16755_10030.pdf
	//allocate memory for overlap matrix eigen stuff
	std::shared_ptr<Matrix> Sevec (new Matrix(nso, nso));
	std::shared_ptr<Vector> Seval (new Vector(nso));
	//compute Sevec and Seval
	//NOTE: the diagonalize command doesn't actually diagonalize
	S->diagonalize(Sevec, Seval);

	//build $S^{-1/2}$ symmetric orthogonalization Matrix (Shalf)
	//NOTE: Shalf is actually $S^{-1/2}$, not $sqrt(S)$
	std::shared_ptr<Matrix> Shalf = (std::shared_ptr<Matrix>) (new Matrix(nso, nso));
	for (int mu = 0; mu<nso; mu++){
	  Shalf->pointer()[mu][mu] = 1.0/(sqrt(Seval->pointer()[mu]));
	}
	//transform Seval back to the non orthogonal Matrix: $Sevec^{-1} Shalf Sevec
	Shalf->back_transform(Sevec);

	/* XXX STEP FOUR: GUESS INITIAL FOCK, ENERGY XXX*/

	//set the intial fock as the core Hamiltonian
	Fa->copy(h);
	//form F' = (S^{-1/2})^T F S^{-1/2}
	std::shared_ptr<Matrix> Fprime (new Matrix(Fa));
	//$F' = Shalf^T (Fprime = Fa) Shalf$
	Fprime->transform(Shalf);

	//allocate memory for F' eigen vectors and values
	//NOTE: Fevec and Feval are for F' NOT Fa
	std::shared_ptr<Matrix> Fevec (new Matrix(nso, nso));
	std::shared_ptr<Vector> Feval (new Vector(nso));
	/* This is the F' from the Roothan equations, and so Fevec are
		really the columns of C' (coefficients matrix) and Feval are really the
		the MO energies given by $\varepsilon$
	 */
	//diagonalize F'  ie, find the evec and evals really.
	Fprime->diagonalize(Fevec, Feval, ascending);

	//For one index, find Ca = Shalf * Cprime
	//Gemm is just a general matrix multiplier -- false so that no transpose
	//taken, Shalf is the left one, Fevec the right
	Ca->gemm(false, false, 1.0, Shalf, Fevec, 0.0);
	//Construct the density: $C^T * C + F', to understand this function:
	//http://www.netlib.org/lapack/
	C_DGEMM('n', 't', nso, nso, na, 1.0, &(Ca->pointer()[0][0]), nso,\
	&(Ca->pointer()[0][0]), nso, 0.0, &(Fprime->pointer()[0][0]), nso);
	Da->copy(Fprime);


	/* FIXME The problem here is that the coefficients matrix is not right
	   (3 columns have switched signs). But the density matrix is okay.
	 */

	/* XXX STEP FIVE: SCF XXX */

	//get nuclear energy
	double e_nuc = mol->nuclear_repulsion_energy();
	//Initialize energy with nuclear energy and the current hamiltonian and fock
	double e_current = e_nuc;
	e_current += Da->vector_dot(h);
	e_current += Da->vector_dot(Fa);

	//e_last is a tracker for energy between iterations
	double e_last = 0.0;
	//dele is the difference in energies between iterations
	double dele = 0.0;
	//deld is the difference in density (as calculated below) between iterations
	double deld = 0.0;
	//iteration counter
	int iter = 0;

	//construct an intermediate tensor for JK matrices:
	std::shared_ptr <Vector> JKI (new Vector(nso* nso*nQ));

	do {
	  e_last = e_current;

	  //construct the coloumb terms
	  C_DGEMV('n',nQ, nso*nso, 1.0, &(Qso->pointer()[0][0]), nso*nso, \
	  &(Da->pointer()[0][0]), 1, 0.0, &(JKI->pointer()[0]), 1);
	  C_DGEMV('t',nQ, nso*nso, 2.0, &(Qso->pointer()[0][0]), nso*nso, \
	  &(JKI->pointer()[0]), 1, 0.0, &(Fa->pointer()[0][0]), 1);

	  //construct the exchange terms
	  C_DGEMM('n','t', nQ*nso, nso, nso, 1.0, &(Qso->pointer()[0][0]), nso, \
	  &(Da->pointer()[0][0]), nso, 0.0, &(JKI->pointer()[0]), nso);
	  for(int Q = 0; Q< nQ; Q++){
		for (int nu = 0; nu < nso; nu ++){
		  C_DCOPY(nso, &(JKI->pointer()[Q*nso*nso + nu]), nso, &(S->pointer()[nu][0]),1);
		}
		C_DCOPY(nso*nso, &(S->pointer()[0][0]), 1, &(JKI->pointer()[Q*nso*nso]), 1);
	  }
	  C_DGEMM('t', 'n', nso, nso, nQ*nso, -1.0, &(Qso->pointer()[0][0]), nso, \
	  &(JKI->pointer()[0]), nso, 1.0, &(Fa->pointer()[0][0]), nso);
	  Fa->add(h);
	  Fprime->copy(Fa);
	  //doing the same things iteratively now
	  //$F' = Shalf^T (Fprime = Fa) Shalf$
	  Fprime->transform(Shalf);
	  //get the new eigvectors
	  Fprime->diagonalize(Fevec, Feval, ascending);
	  //get new coeffs C = Shalf * C'
	  Ca->gemm(false, false, 1.0, Shalf, Fevec, 0);

	  //calculate density from Ca
	  deld = 0.0;
	  C_DGEMM('n', 't', nso, nso, na, 1.0, &(Ca->pointer()[0][0]), nso,\
	  &(Ca->pointer()[0][0]), nso, 0.0, &(Fprime->pointer()[0][0]), nso);
	  //update del d
	  for (int mu  = 0; mu < nso; mu++){
		  for (int nu = 0; nu < nso; nu++){
			  double dum = Da->pointer()[mu][nu]-Fprime->pointer()[nu][mu];
			  deld += dum * dum;
		  }
	  }
	  Da->copy(Fprime);
	  deld = sqrt(deld);

	  //update energy
	  e_current = e_nuc;
	  e_current += Da->vector_dot(h);
	  e_current += Da->vector_dot(Fa);

	  //update del e
	  dele = fabs(e_last-e_current);
	  //print everything in a nice pretty table
	  outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,e_current,dele,deld);
	  iter++;
	  if (iter > maxiter){
		  break;
	  }
	} while(deld > d_convergence || dele > e_convergence);

	if (iter > maxiter){
		throw PsiException("Maximum allowed iterations exceeded", __FILE__, __LINE__);
	}

	outfile->Printf("\n");
	outfile->Printf("    SCF iterations converged!\n");
	outfile->Printf("\n");
	outfile->Printf("    * SCF total energy: %20.12lf\n",e_current);
	Process::environment.globals["SCF TOTAL ENERGY"] = e_current;

	// Typically you would build a new wavefunction and populate it with data
	return ref_wfn;
}

}} // End namespaces

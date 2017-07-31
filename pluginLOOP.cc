
/*
 * @BEGIN LICENSE
 *
 * myscf by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
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

    /* Your code goes here */

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    myscf                                            *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    A restricted Hartree-Fock plugin to Psi4         *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Eugene DePrince                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    // make sure we are running in c1 symmetry
    if ( ref_wfn->nirrep() > 1 ) {
        throw PsiException("plugin myscf only works with c1 symmetry. Set symmetry c1 in the molecule group in your input file.",__FILE__,__LINE__);
    }

    // grab the one-electron integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(ref_wfn));

    // one-electron kinetic energy integrals
    std::shared_ptr<Matrix> T = mints->so_kinetic();

    // one-electron potential energy integrals
    std::shared_ptr<Matrix> V = mints->so_potential();

    // overlap integrals
    std::shared_ptr<Matrix> S = mints->so_overlap();

    // build the core hamiltonian
    std::shared_ptr<Matrix> h = (std::shared_ptr<Matrix>)(new Matrix(T));
    h->add(V);

    // grab the molecule from the wavefunction that was passed into the plugin
    std::shared_ptr<Molecule> mol = ref_wfn->molecule();

    // get primary basis:
    std::shared_ptr<BasisSet> primary = ref_wfn->get_basisset("ORBITAL");

    // total number of basis functions
    int nso = primary->nbf();

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = ref_wfn->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    int nQ = auxiliary->nbf();

    // grab some input options
    double e_convergence = options.get_double("E_CONVERGENCE");
    double d_convergence = options.get_double("D_CONVERGENCE");
    int maxiter          = options.get_int("MAXITER");

    // use the molecule to determine the total number of electrons
    int charge     = mol->molecular_charge();
    int nelectron  = 0;
    for (int i = 0; i < mol->natom(); i++) {
        nelectron += (int)mol->Z(i);
    }
    nelectron -= charge;

    // this code only works for closed shells
    if ( nelectron % 2 != 0 ) {
        throw PsiException("plugin myscf only works for closed shells",__FILE__,__LINE__);
    }

    // the number of alpha electrons
    int na = nelectron / 2;

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
    outfile->Printf("    No. electrons:                  %5i\n",nelectron);
    outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
    outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // allocate memory for coefficients, density, fock matrix
    std::shared_ptr<Matrix> Ca = (std::shared_ptr<Matrix>)(new Matrix(nso,nso));
    std::shared_ptr<Matrix> Da = (std::shared_ptr<Matrix>)(new Matrix(nso,nso));
    std::shared_ptr<Matrix> Fa = (std::shared_ptr<Matrix>)(new Matrix(nso,nso));

    // construct the three-index integrals
    // since we want the SO-basis integrals, it is fine to pass empty Ca matrix
    // similarly, the number of active vs inactive orbitals isn't really important here.
    std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca,na,nso-na,na,nso-na,options));
    std::shared_ptr<Matrix> Qso = DF->Qso();
    double ** Qp = Qso->pointer();

    // allocate memory for eigenvectors and eigenvalues of the overlap matrix
    std::shared_ptr<Matrix> Sevec ( new Matrix(nso,nso) );
    std::shared_ptr<Vector> Seval ( new Vector(nso) );

    // build S^(-1/2) symmetric orthogonalization matrix
    S->diagonalize(Sevec,Seval);
    for (int mu = 0; mu < nso; mu++) {
        Seval->pointer()[mu] = 1.0 / sqrt(Seval->pointer()[mu]);
    }

    std::shared_ptr<Matrix> Shalf = (std::shared_ptr<Matrix>)( new Matrix(S) );
    double ** sp = Shalf->pointer();
    // transform Seval back to nonorthogonal basis
    for (int mu = 0; mu < nso; mu++) {
        for (int nu = 0; nu < nso; nu++) {
            double dum = 0.0;
            for (int i = 0; i < nso; i++) {
                dum += Seval->pointer()[i] * Sevec->pointer()[nu][i] * Sevec->pointer()[mu][i];
            }
            sp[mu][nu] = dum;
        }
    }

    // form F' = ST^(-1/2) F S^(-1/2), where F = h
    Fa->copy(h);
    std::shared_ptr<Matrix> Fprime ( new Matrix(Fa) );

    double ** fp  = Fa->pointer();
    double ** fpp = Fprime->pointer();

    for (int i = 0; i < nso; i++) {
        for (int j = 0; j < nso; j++) {
            double dum = 0.0;
            for (int mu = 0; mu < nso; mu++) {
                for (int nu = 0; nu < nso; nu++) {
                    dum += fp[mu][nu] * sp[mu][i] * sp[nu][j];
                }
            }
            fpp[i][j] = dum;
        }
    }


    // allocate memory for eigenvectors and eigenvalues of F'
    std::shared_ptr<Matrix> Fevec ( new Matrix(nso,nso) );
    std::shared_ptr<Vector> Feval ( new Vector(nso) );

    // diagonalize F' to obtain C'
    Fprime->diagonalize(Fevec,Feval,ascending);

    // Find C = S^(-1/2)C'
    double ** cp = Ca->pointer();
    for (int mu = 0; mu < nso; mu++) {
        for (int i = 0; i < nso; i++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso; nu++) {
                dum += sp[nu][mu] * Fevec->pointer()[nu][i];
            }
            cp[mu][i] = dum;
        }
    }
    //outfile->Printf("This is Da 1\n");
	//Da->print();
    // Construct density from C
    double ** dp = Da->pointer();
    for (int mu = 0; mu < nso; mu++) {
        for (int nu = 0; nu < nso; nu++) {
            double dum = 0.0;
            for (int i = 0; i < na; i++) {
                dum += cp[mu][i] * cp[nu][i];
            }
            dp[mu][nu] = dum;
        }
    }
    //outfile->Printf("This is Da 2\n");
	//Da->print();
    //Da->copy(Fprime);
    //outfile->Printf("This is Da 3\n");
	//Da->print();

    //Fprime->print();

    // initial energy, E = D(H+F) + Enuc
    double e_nuc = mol->nuclear_repulsion_energy();

    double e_current = e_nuc;
    e_current       += Da->vector_dot(h);
    e_current       += Da->vector_dot(Fa);

    //  SCF iterations

    double e_last    = 0.0;
    double dele      = 0.0;
    double deld      = 0.0;

    outfile->Printf("\n");
    outfile->Printf("    Guess energy:  %20.12lf\n",e_current);
    outfile->Printf("\n");
    outfile->Printf("    ==>  Begin SCF Iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("              energy ");
    outfile->Printf("                  dE ");
    outfile->Printf("                  dD ");
    outfile->Printf("\n");

    int iter = 0;

    // intermediate tensor for JK construction
    std::shared_ptr<Vector> JKI (new Vector(nso*nso*nQ) );
    double * Ip = JKI->pointer();

    do {

        e_last = e_current;

        // form J
        for (int Q = 0; Q < nQ; Q++) {
            double dum = 0.0;
            for (int lam = 0; lam < nso; lam++) {
                for (int sig = 0; sig < nso; sig++) {
                    dum += dp[lam][sig] * Qp[Q][lam*nso+sig];
                }
            }
            Ip[Q] = dum;
        }
        for (int mu = 0; mu < nso; mu++) {
            for (int nu = 0; nu < nso; nu++) {
                double dum = 0.0;
                for (int Q = 0; Q < nQ; Q++) {
                    dum += Ip[Q] * Qp[Q][mu*nso+nu];
                }
                fp[mu][nu] = 2.0 * dum;
            }
        }

        // form k
        for (int nu = 0; nu < nso; nu++) {
            for (int sig = 0; sig < nso; sig++) {
                for (int Q = 0; Q < nQ; Q++) {
                    double dum = 0.0;
                    for (int lam = 0; lam < nso; lam++) {
                        dum += dp[lam][sig] * Qp[Q][lam*nso+nu];
                    }
                    Ip[nu*nso*nQ+sig*nQ+Q] =dum;
                }
            }
        }
        for (int mu = 0; mu < nso; mu++) {
            for (int nu = 0; nu < nso; nu++) {
                double dum = 0.0;
                for (int sig = 0; sig < nso; sig++) {
                    for (int Q = 0; Q < nQ; Q++) {
                        dum += Ip[nu*nso*nQ+sig*nQ+Q] * Qp[Q][mu*nso+sig];
                    }
                }
                fp[mu][nu] -= dum;
            }
        }

        // form F = h + J - K
        Fa->add(h);

        // form F' = ST^(-1/2) F S^(-1/2)
        for (int i = 0; i < nso; i++) {
            for (int j = 0; j < nso; j++) {
                double dum = 0.0;
                for (int mu = 0; mu < nso; mu++) {
                    for (int nu = 0; nu < nso; nu++) {
                        dum += fp[mu][nu] * sp[mu][i] * sp[nu][j];
                    }
                }
                fpp[i][j] = dum;
            }
        }

        // diagonalize F' to obtain C'
        Fprime->diagonalize(Fevec,Feval,ascending);

        // Find C = S^(-1/2)C'
        for (int mu = 0; mu < nso; mu++) {
            for (int i = 0; i < nso; i++) {
                double dum = 0.0;
                for (int nu = 0; nu < nso; nu++) {
                    dum += sp[nu][mu] * Fevec->pointer()[nu][i];
                }
                cp[mu][i] = dum;
            }
        }

        // Construct density from C
        deld = 0.0;
        for (int mu = 0; mu < nso; mu++) {
            for (int nu = 0; nu < nso; nu++) {
                double dum = 0.0;
                for (int i = 0; i < na; i++) {
                    dum += cp[mu][i] * cp[nu][i];
                }
                double dum2 = dp[mu][nu] - dum;
                deld += dum2*dum2;
                dp[mu][nu] = dum;
            }
        }
        deld = sqrt(deld);

        // current energy, E = D(H+F) + Enuc
        e_current  = e_nuc;
        e_current += Da->vector_dot(h);
        e_current += Da->vector_dot(Fa);

        // dele
        dele = fabs(e_last - e_current);

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,e_current,dele,deld);

        iter++;
        if ( iter > maxiter ) break;

    }while(dele > e_convergence || deld > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }
    outfile->Printf("\n");
    outfile->Printf("    SCF iterations converged!\n");
    outfile->Printf("\n");
    outfile->Printf("    * SCF total energy: %20.12lf\n",e_current);

    Process::environment.globals["SCF TOTAL ENERGY"] = e_current;

    // Typically you would build a new wavefunction and populate it with data
    // (note this is just the original wavefunction passed into the module)
    return ref_wfn;
}

}} // End namespaces

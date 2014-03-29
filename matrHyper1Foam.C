/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    hyper1Foam

Description
    Solves a transport equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
 
#include "complex.H"
 
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "psi before = " << psi << nl << endl;
	Info<< "phi before = " << phi << nl << endl;
        
	fvScalarMatrix psiEqn
        (
	    fvm::ddt(psi)
	    +
	    fvm::div(phi,psi)
        );
	



// ===============================================================================

//                 arma::mat A = arma::zeros<arma::mat>(psi.size(),psi.size());
                List<List<scalar> > A;
// 		arma::mat b = arma::zeros<arma::mat>(psi.size(),psi.size());
		List<scalar> b;
// 		upSum = p.size();
		
		A.resize(psi.size());					// A and p have equal size
                b.resize(psi.size());
                forAll(A, i)
                {
            	  A[i].resize(psi.size());				// every element of A is array
            	  forAll(A[i],j)					// clearing A and b
            	  {
            	    A[i][j] = 0.0;
		  }
            	  b[i] = 0.0;
                }
		
		
                forAll(psi,i) // forAll(A,i)
                {
           	    A[i][i] = psiEqn.diag()[i];
//             	    A(i,i) = psiEqn.diag()[i];
// // 		    Info << "A("<<i<<","<< i << ") ===" << A(i,i) << nl << endl;
           	    b[i]    = psiEqn.source()[i];
// 		    b(i,0) = psiEqn.source()[i];
// // 		    Info << "b("<<i<<","<< 0 << ") ===" << b(i,0) << nl << endl;
		  
		}
                
                const lduAddressing& addr = psiEqn.lduAddr();
                const labelList& lowerAddr = addr.lowerAddr();
                const labelList& upperAddr = addr.upperAddr();                
                
                forAll(lowerAddr, i)
                {
// 		    A(lowerAddr[i],upperAddr[i]) = psiEqn.upper()[i];
// // 		    Info << "A(" << lowerAddr[i] << "," << upperAddr[i] << ")=" << psiEqn.upper()[i] << nl << endl;
// 		    A(upperAddr[i],lowerAddr[i]) = psiEqn.lower()[i];
// // 		    Info << "A(" << upperAddr[i] << "," << lowerAddr[i] << ")=" << psiEqn.lower()[i] << nl << endl;
           	    A[lowerAddr[i]][upperAddr[i]] = psiEqn.upper()[i];
           	    A[upperAddr[i]][lowerAddr[i]] = psiEqn.lower()[i];           	                	    
// // //             	    downSum += pEqn.upper()[i]* pEqn.upper()[i];
////////            	    downSum += pEqn.lower()[i]*pEqn.lower()[i];    
                }
                
                forAll(psi.boundaryField(),I) // what is it??
                {
            	    const fvPatch &ptch=psi.boundaryField()[I].patch();
		    Info << "psi.boundaryField()[" << I << "]" << psi.boundaryField()[I];
            	    forAll(ptch,J)
            	    {
           		int w=ptch.faceCells()[J];
// // 			Info << "ptch.faceCells() = " << ptch.faceCells() << nl << endl;
           		A[w][w]+=psiEqn.internalCoeffs()[I][J];
// 			A(w,w)+=psiEqn.internalCoeffs()[I][J];
// // 			Info << "A(" << w << "," << w << ") = " << A(w,w) << nl << endl;
           		b[w]   +=psiEqn.boundaryCoeffs()[I][J];
// 			b(w,0) +=psiEqn.boundaryCoeffs()[I][J];
// // 			Info << "b(" << w << "," << 0 << ") = " << b(w,0) << nl << endl;
            	    }
                
                }
                Info << "psi.boundaryField() = " << psi.boundaryField() << nl << endl;
		Info << "=== A(i,j) ===" << nl << endl;
		for (int i=0; i<psi.size(); i++) // forAll(A,i)
                {
		  for (int j=0; j<psi.size(); j++)
		  {
		    Info << A[i][j] << " ";
		  }
		  Info << nl << endl;
		  
		}
		
		
		Info << "=== b(i,j) ===" << nl << endl;
		for (int i=0; i<psi.size(); i++) // forAll(A,i)
                {
// 		  for (int j=0; j<psi.size(); j++)
// 		  {
		    Info << b[i] << " ";
// 		  }
		  Info << nl << endl;
		  
		}
// 		arma::inplace_trans(b);
// 		arma::mat Ai = arma::inv(A);
// 		arma::mat bt = trans(b);
// 		arma::mat X = inv(A)*b;
// 		Info << "=== Ai(i,j) ===" << nl << endl;
// 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
//                 {
// 		  for (int j=0; j<psi.size(); j++)
// 		  {
// 		    Info << Ai(i,j) << " ";
// 		  }
// 		  Info << nl << endl;
// 		  
// 		}
// 		Info << "=== X(i,j) ===" << nl << endl;
// 		for (int i=0; i<psi.size(); i++) // forAll(A,i)
//                 {
// 		  for (int j=0; j<psi.size(); j++)
// 		  {
// 		    Info << X(i,j) << " ";
// 		  }
// 		  Info << nl << endl;
// 		  
// 		}
		
// ===================================================================================
	
        psiEqn.solve();	
	
	Info<< "psiEqn = " << psiEqn << nl << endl;
	Info<< "psi after = " << psi << nl << endl;
	Info<< "phi after = " << phi << nl << endl;
	Info<< "fvm::ddt(psi) = " << fvm::ddt(psi) << nl << endl;
	Info<< "fvm::div(phi,psi) = " << fvm::div(phi,psi) << nl << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

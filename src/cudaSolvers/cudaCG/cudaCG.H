/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    cudaCG

Description
    Preconditioned Conjugate Gradient solver with run-time selectable
    preconditioning for use with GPU

Author
    Dominik Christ, Wikki Ltd.
    Based on Cufflink library by Daniel P. Combest

SourceFiles
    cudaCG.C

\*---------------------------------------------------------------------------*/

#ifndef cudaCG_H
#define cudaCG_H

#include "cudaSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
			     Class cudaCG Declaration
\*---------------------------------------------------------------------------*/

class cudaCG
:
    public cudaSolver
{
    // Private Member Functions

	//- Disallow default bitwise copy construct
	cudaCG(const cudaCG&);

	//- Disallow default bitwise assignment
	void operator=(const cudaCG&);


public:

    //- Runtime type information
    TypeName("cudaCG");


    // Constructors

	//- Construct from matrix components and solver data stream
	cudaCG
	(
	    const word& fieldName,
	    const lduMatrix& matrix,
	    const FieldField<Field, scalar>& coupleBouCoeffs,
	    const FieldField<Field, scalar>& coupleIntCoeffs,
	    const lduInterfaceFieldPtrsList& interfaces,
	    const dictionary& dict
	);


    // Destructor

	virtual ~cudaCG()
	{}


    // Member Functions

	//- Solve the matrix with this solver
	virtual lduSolverPerformance solve
	(
	    scalarField& x,
	    const scalarField& b,
	    const direction cmpt = 0
	) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

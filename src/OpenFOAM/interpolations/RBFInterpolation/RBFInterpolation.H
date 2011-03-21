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
    RBFInterpolation

Description
    Radial basis function interpolation class

Description
    Interpolation class which uses Radial Basis Functions to interpolate the
    fluid displacements for given boundary displacements.

	The coefficient vectors, alpha and beta are determined by solving
    the system:

	| db | = | Mbb Pb | | alpha |
	| 0  |   | Pb  0  | |  beta |

	where db are the given boundary displacements,
	Mbb the boundary RBF correlation matrix (NbxNb), containing RBF evaluations
    at the boundary nodes, and Pb some linear polynomial matrix (Nbx4).

	Those coefficients are calculated every timestep, with the current
    boundary displacements db, with the inverse of Mbb. Using those
    coefficients, the RBF is evaluated at all fluid points every
    timestep.

	The efficiency of this method is increased by:
		1) using control points which is a subset of the moving
           boundary points. Those control points are selected by
		   a coarsening function.
		2) The outer boundary points are neglected since a cutoff function
           is used toward the outer boundaries.

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

SourceFiles
    RBFInterpolation.C
    RBFInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include "dictionary.H"
#include "RBFFunction.H"
#include "point.H"
#include "Switch.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class RBFInterpolation Declaration
\*---------------------------------------------------------------------------*/

class RBFInterpolation
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Reference to control points
        const vectorField& controlPoints_;

        //- Rerefence to all points
        const vectorField& allPoints_;

        //- RBF function
        autoPtr<RBFFunction> RBF_;

        //- Interpolation matrix
        mutable scalarSquareMatrix* BPtr_;

        //- Focal point for cut-off radii
        point focalPoint_;

        //- Inner cut-off radius
        scalar innerRadius_;

        //- Outer cut-off radius
        scalar outerRadius_;

        //- Add polynomials to RBF matrix
        Switch polynomials_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RBFInterpolation&);


        //- Return interpolation matrix
        const scalarSquareMatrix& B() const;

        //- Calculate interpolation matrix
        void calcB() const;

        //- Clear out
        void clearOut();


public:

    // Constructors

        //- Construct from components
        RBFInterpolation
        (
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& allPoints
        );

        //- Construct as copy
        RBFInterpolation(const RBFInterpolation&);


    // Destructor

        ~RBFInterpolation();


    // Member Functions

        //- Interpolate
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& ctrlField) const;

        //- Move points
        void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RBFInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

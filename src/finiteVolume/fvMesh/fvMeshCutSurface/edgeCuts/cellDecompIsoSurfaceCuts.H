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
    Foam::cellDecompIsoSurfaceCuts

Description
    cellCuts resulting from isoSurface reconstruction on cellDecomp tet mesh

SourceFiles
    cellDecompIsoSurfaceCuts.C

\*---------------------------------------------------------------------------*/

#ifndef cellDecompIsoSurfaceCuts_H
#define cellDecompIsoSurfaceCuts_H

#include "cellDecompCuts.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class cellDecompIsoSurfaceCuts Declaration
\*---------------------------------------------------------------------------*/

class cellDecompIsoSurfaceCuts
:
    public cellDecompCuts
{
    // Private Member Functions

        //- Do all the hard work: determine edge crossings
        void constructEdgeCuts
        (
            const volScalarField&,
            const pointScalarField&,
            const scalar,
            const scalar
        );

public:

    // Constructors

        //- Construct from field on cells and field on vertices
        cellDecompIsoSurfaceCuts
        (
            const volScalarField&,
            const pointScalarField&,
            const scalar isoValue,
            const scalar tol
        );

        //- Construct from interpolator and field (does interpolation)
        cellDecompIsoSurfaceCuts
        (
            const volScalarField&,
            const volPointInterpolation&,
            const scalar isoValue,
            const scalar tol
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

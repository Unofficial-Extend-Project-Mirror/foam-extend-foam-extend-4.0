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
    Foam::scotchDecomp

Description
    Scotch domain decomposition

SourceFiles
    scotchDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef scotchDecomp_H
#define scotchDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scotchDecomp Declaration
\*---------------------------------------------------------------------------*/

class scotchDecomp
:
    public decompositionMethod
{
    // Private data

        //- Mesh reference
        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        scotchDecomp(const scotchDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const scotchDecomp&);


        //- Check and print error message
        static void check(const int, const char*);

        label decompose
        (
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& cWeights,
            List<int>& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("scotch");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        scotchDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~scotchDecomp()
        {}


    // Member Functions

        //- Return mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        virtual bool parallelAware() const
        {
            // Scotch does not know about proc boundaries
            return false;
        }

        //- Decompose cells with weights
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Decompose cell clusters with weights on clusters
        virtual labelList decompose
        (
            const labelList& fineToCoarse,
            const pointField& coarsePoints,
            const scalarField& coarseWeights
        );

        //- Decompose cells with weights with explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

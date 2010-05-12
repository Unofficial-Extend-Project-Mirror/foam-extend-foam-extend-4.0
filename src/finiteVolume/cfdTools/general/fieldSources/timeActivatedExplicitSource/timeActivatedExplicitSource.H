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
    Foam::timeActivatedExplicitSource

Description
    Creates a cell set source that is activated at timeStart_ for duration_

SourceFiles
    timeActivatedExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef timeActivatedExplicitSource_H
#define timeActivatedExplicitSource_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class timeActivatedExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedExplicitSource
{
    // Private data

        //- Properties dictionary
        IOdictionary dict_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to time database
        const Time& runTime_;

        //- Name of cell source
        word cellSource_;

        //- Time start [s]
        scalar timeStart_;

        //- Duration [s]
        scalar duration_;

        //- Value when "on"
        dimensionedScalar onValue_;

        //- Value when "off"
        dimensionedScalar offValue_;

        //- Current source value
        dimensionedScalar currentValue_;

        //- The method by which the cells will be selecetd
        autoPtr<topoSetSource> cellSelector_;

        //- The set of selected cells
        cellSet selectedCellSet_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeActivatedExplicitSource(const timeActivatedExplicitSource&);

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedExplicitSource&);


public:

    // Constructors

        //- Construct from explicit source name and mesh
        timeActivatedExplicitSource(const word&, const fvMesh&);


    // Member Functions

        // Access

            //- Return the start time
            scalar timeStart() const;

            //- Return the duration
            scalar duration() const;

            //- Return the current value of the source
            const dimensionedScalar& currentValue() const;

            //- Return a tmp field of the source
            tmp<DimensionedField<scalar, volMesh> > Su() const;


        //- Update
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

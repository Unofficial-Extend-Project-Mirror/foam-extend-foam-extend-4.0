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
    Foam::molConfig

Description

SourceFiles
    molConfigI.H
    molConfig.C
    molConfigIO.C

\*---------------------------------------------------------------------------*/

#ifndef molConfig_H
#define molConfig_H

#include "labelVector.H"
#include "scalar.H"
#include "vector.H"
#include "labelField.H"
#include "scalarField.H"
#include "vectorField.H"
#include "IOField.H"
#include "EulerCoordinateRotation.H"
#include "Random.H"

#include "Time.H"
#include "IOdictionary.H"
#include "IOstreams.H"
#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class molConfig Declaration
\*---------------------------------------------------------------------------*/

class molConfig
{
    // Private data

        const IOdictionary& molConfigDescription_;

        const polyMesh& mesh_;

        DynamicList<word> idList_;

        labelField id_;

        scalarField mass_;

        vectorField positions_;

        labelField cells_;

        vectorField U_;

        vectorField A_;

        labelField tethered_;

        vectorField tetherPositions_;

        label nMol_;


public:

    // Constructors

        //- Construct from IOdictionary and mesh
        molConfig(IOdictionary&, const polyMesh&);


    // Destructor

        ~molConfig();


    // Member Functions

        void createMolecules();


        // Access

            inline const List<word>& molIdList() const;

            inline const labelField& id() const;

            inline const scalarField& mass() const;

            inline const vectorField& positions() const;

            inline const labelField& cells() const;

            inline const vectorField& U() const;

            inline const vectorField& A() const;

            inline const labelField& tethered() const;

            inline const vectorField& tetherPositions() const;

            inline label nMol() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "molConfigI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

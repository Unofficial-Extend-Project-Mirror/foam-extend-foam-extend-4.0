/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    simpleEnginePiston

Description

SourceFiles
    simpleEnginePiston.C

\*---------------------------------------------------------------------------*/

#ifndef simpleEnginePiston_H
#define simpleEnginePiston_H

#include "polyPatchID.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                           Class simpleEnginePiston Declaration
\*---------------------------------------------------------------------------*/

class simpleEnginePiston
{
    // Private data

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Piston patch
        polyPatchID patchID_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Piston layering data

            //- Min layer thickness
            const scalar minLayer_;

            //- Max layer thickness
            const scalar maxLayer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleEnginePiston(const simpleEnginePiston&);

        //- Disallow default bitwise assignment
        void operator=(const simpleEnginePiston&);


public:

    // Static data members


    // Constructors

        //- Construct from components
        simpleEnginePiston
        (
            const polyMesh& mesh,
            const word& pistonPatchName,
            const autoPtr<coordinateSystem>& pistonCS,
            const scalar minLayer,
            const scalar maxLayer

        );

        //- Construct from dictionary
        simpleEnginePiston
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return ID of piston patch
        const polyPatchID& patchID() const
        {
            return patchID_;
        }

        // Piston layering thickness

            scalar minLayer() const
            {
                return minLayer_;
            }

            scalar maxLayer() const
            {
                return maxLayer_;
            }


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

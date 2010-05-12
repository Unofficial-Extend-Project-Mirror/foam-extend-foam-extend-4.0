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
    Foam::enginePiston

Description
    Foam::enginePiston

SourceFiles
    enginePiston.C

\*---------------------------------------------------------------------------*/

#ifndef enginePiston_H
#define enginePiston_H

#include "polyPatchID.H"
#include "coordinateSystem.H"
#include "simpleEnginePiston.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                           Class enginePiston Declaration
\*---------------------------------------------------------------------------*/

class enginePiston
:
    public simpleEnginePiston
{
    // Private data

        //- Name of piston point set
        word pistonPointSetName_;

        //- Name of piston face set
        word pistonFaceSetName_;

        //- Name of piston cell set
        word pistonCellSetName_;


        //- Bowl in piston patch
        polyPatchID bowlInPistonPatchID_;

        //- Bowl in cylinder patch
        polyPatchID bowlInCylinderPatchID_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        enginePiston(const enginePiston&);

        //- Disallow default bitwise assignment
        void operator=(const enginePiston&);


public:

    // Static data members


    // Constructors

        //- Construct from components
        enginePiston
        (
            const polyMesh& mesh,
            const word& pistonPatchName,
            const autoPtr<coordinateSystem>& pistonCS,
            const scalar minLayer,
            const scalar maxLayer,
            const word& pistonPointSetName,
            const word& pistonFaceSetName,
            const word& pistonCellSetName,
            const word& bowlInPistonPatchName,
            const word& bowlInCylinderPatchName
        );

        //- Construct from dictionary
        enginePiston
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        const word& pistonPointSetName() const
        {
            return pistonPointSetName_;
        }

        const word& pistonFaceSetName() const
        {
            return pistonFaceSetName_;
        }

        const word& pistonCellSetName() const
        {
            return pistonCellSetName_;
        }


        const polyPatchID& bowlInPistonPatchID() const
        {
                return bowlInPistonPatchID_;
        }

        const polyPatchID& bowlInCylinderPatchID() const
        {
            return bowlInCylinderPatchID_;
        }


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

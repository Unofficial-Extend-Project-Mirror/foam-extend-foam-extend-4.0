/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::faceSet

Description
    A list of face labels.

SourceFiles
    faceSet.C

\*---------------------------------------------------------------------------*/

#ifndef faceSet_H
#define faceSet_H

#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceSet Declaration
\*---------------------------------------------------------------------------*/

class faceSet
:
    public topoSet
{
public:

    //- Runtime type information
    TypeName("faceSet");


    // Constructors


        //- Construct from IOobject
        faceSet(const IOobject& obj);

        //- Construct from objectRegistry and name
        faceSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r = MUST_READ,
            writeOption w = NO_WRITE
        );

        //- Construct from additional size of labelHashSet
        faceSet
        (
            const polyMesh& mesh,
            const word& name,
            const label,
            writeOption w = NO_WRITE
        );

        //- Construct from existing set
        faceSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet&,
            writeOption w = NO_WRITE
        );

        //- Construct from additional labelHashSet
        faceSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet&,
            writeOption w = NO_WRITE
        );


    // Destructor

        virtual ~faceSet();


    // Member functions

        //- Sync faceSet across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

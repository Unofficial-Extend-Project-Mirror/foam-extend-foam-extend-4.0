/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    componentReference

Description
    Class contains data for a component reference as used in stress
    analysis solvers.

SourceFiles
    componentReferenceI.H
    componentReference.C
    componentReferenceIO.C

\*---------------------------------------------------------------------------*/

#ifndef componentReference_H
#define componentReference_H

#include "polyPatchID.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class componentReference Declaration
\*---------------------------------------------------------------------------*/

class componentReference
{
    // Private data

        //- Patch ID
        polyPatchID patchID_;

        //- Face index
        label faceIndex_;

        //- Direction
        direction dir_;

        //- Value in direction
        scalar value_;

    // Private Member Functions

        //- Create direction given a name
        direction getDir(const dictionary& dict) const
        {
            word dirName(dict.lookup("direction"));

            if (dirName == "x" || dirName == "X")
            {
                return vector::X;
            }
            else if (dirName == "y" || dirName == "Y")
            {
                return vector::Y;
            }
            else if (dirName == "z" || dirName == "Z")
            {
                return vector::Z;
            }
            else
            {
                FatalIOErrorIn
                (
                    "vector::component componentReference::getComp("
                    "const word& dirName) const",
                    dict
                )   << "Direction " << dirName << " not recognised.  Please "
                    << "use x, y or z" << abort(FatalIOError);

                // Dummy return to keep compiler happy
                return vector::X;
            }
        }


        //- Check if patch face is in range
        void checkPatchFace(const fvMesh& mesh) const
        {
            if
            (
                !patchID_.active()
             || faceIndex_ >= mesh.boundaryMesh()[patchID_.index()].size()
            )
            {
                FatalErrorIn
                (
                    "void checkPatchFace(const componentReference::fvMesh&)"
                    "const"
                )   << "Non-existing patch or index out of range."
                    << abort(FatalError);
            }
        }


public:

    // Public classes

        //- Class used for the read-construction of
        //  PtrLists of componentReference
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<componentReference> operator()(Istream& is) const
            {
                dictionary crDict(is);

                autoPtr<componentReference> cr
                (
                    new componentReference(mesh_, crDict)
                );

                return cr;
            }
        };


    // Constructors

        //- Construct from components
        componentReference
        (
            const fvMesh& mesh,
            const word& patchName,
            const label faceIndex,
            const direction dir,
            const scalar value
        )
        :
            patchID_(patchName, mesh.boundaryMesh()),
            faceIndex_(faceIndex),
            dir_(dir),
            value_(value)
        {
            checkPatchFace(mesh);
        }
        

        //- Construct from dictionary
        componentReference
        (
            const fvMesh& mesh,
            const dictionary& dict
        )
        :
            patchID_(dict.lookup("patch"), mesh.boundaryMesh()),
            faceIndex_(readLabel(dict.lookup("face"))),
            dir_(getDir(dict)),
            value_(readScalar(dict.lookup("value")))
        {
            checkPatchFace(mesh);
        }

        //- Clone
        autoPtr<componentReference> clone() const
        {
            return autoPtr<componentReference>(new componentReference(*this));
        }


    // Destructor - default


    // Member Functions

        //- Return patch index
        label patchIndex() const
        {
            return patchID_.index();
        }

        //- Return face index
        label faceIndex() const
        {
            return faceIndex_;
        }

        //- Return direction
        direction dir() const
        {
            return dir_;
        }

        //- Return value
        scalar value() const
        {
            return value_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

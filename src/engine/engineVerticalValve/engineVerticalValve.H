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
    engineVerticalValve

Description
    class for engine meshes with vertical valves

SourceFiles
    engineVerticalValve.C

\*---------------------------------------------------------------------------*/

#ifndef engineVerticalValve_H
#define engineVerticalValve_H

#include "word.H"
#include "coordinateSystem.H"
#include "polyPatchID.H"
#include "graph.H"
#include "engineValve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                          Class engineVerticalValve Declaration
\*---------------------------------------------------------------------------*/

class engineVerticalValve
:
    public engineValve
{
    // Private data

        //- Valve head patch
        polyPatchID valveHeadPatch_;

        //- Offset for layer Addition/Removal
        scalar topLayerOffset_;

        //- Tolerance for adding layer face zone
        scalar topLayerTol_;

        //- Offset for layer Addition/Removal
        scalar bottomLayerOffset_;

        //- Tolerance for adding layer face zone
        scalar bottomLayerTol_;

        //- Distance for attach/detach boundary
        scalar detachDistance_;

        //- Tolerance for attach/detach boundary faces
        scalar detachTol_;

        //- Tolerance for deformation/layering
        scalar deformationLift_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        engineVerticalValve(const engineVerticalValve&);

        //- Disallow default bitwise assignment
        void operator=(const engineVerticalValve&);



public:

    // Static data members


    // Constructors

        //- Construct from components
        engineVerticalValve
        (
            const word& name,
            const polyMesh& mesh,
            const autoPtr<coordinateSystem>& valveCS,
            const word& bottomPatchName,
            const word& poppetPatchName,
            const word& stemPatchName,
            const word& curtainInPortPatchName,
            const word& curtainInCylinderPatchName,
            const word& detachInCylinderPatchName,
            const word& detachInPortPatchName,
            const labelList& detachFaces,
            const graph& liftProfile,
            const scalar minLift,
            const scalar minTopLayer,
            const scalar maxTopLayer,
            const scalar minBottomLayer,
            const scalar maxBottomLayer,
            const scalar diameter,
            const word& valveHeadPatchName,
            const scalar topLayerOffset,
            const scalar topLayerTol,
            const scalar bottomLayerOffset,
            const scalar bottomLayerTol,
            const scalar detachDistance,
            const scalar detachTol,
            const scalar deformationLift
        );

        //- Construct from dictionary
        engineVerticalValve
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions


        const scalar& topLayerOffset() const
        {
            return topLayerOffset_;
        }

        const scalar& topLayerTol() const
        {
            return topLayerTol_;
        }

        const scalar& bottomLayerOffset() const
        {
            return bottomLayerOffset_;
        }

        const scalar& bottomLayerTol() const
        {
            return bottomLayerTol_;
        }

        const scalar& detachDistance() const
        {
            return detachDistance_;
        }

        const scalar& detachTol() const
        {
            return detachTol_;
        }

        const scalar& deformationLift() const
        {
            return deformationLift_;
        }

   // Valve patches

       //- Return ID of bottom patch
       const polyPatchID& valveHeadID() const
       {
           return valveHeadPatch_;
       }


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

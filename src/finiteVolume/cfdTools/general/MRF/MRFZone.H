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
    Foam::MRFZone

Description
    MRF zone definition based on cell zone and optional face zone and
     parameters obtained from a control dictionary constructed from the
    given stream.

    The rotation of the MRF region is defined by an origin and axis of
    rotation and an angular speed.

    Coding for compressible flows and generalisation + efficiency improvements
    by Hrvoje Jasak, Wikki Ltd.  All rights reserved

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    MRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZone_H
#define MRFZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MRFZone Declaration
\*---------------------------------------------------------------------------*/

class MRFZone
{
    // Private data

        const fvMesh& mesh_;

        const word name_;

        const dictionary dict_;

        label cellZoneID_;

        //- label of face zone with faces on outside of cell zone.
        label faceZoneID_;

        //- Do patches move with frame (true) or are
        //  explicitly provided (false, old behaviour)
        Switch allPatchesMove_;

        const wordList patchNames_;
        labelList patchLabels_;

        //- Internal faces that are part of MRF
        labelList internalFaces_;

        //- Outside faces (per patch) that move with the MRF
        labelListList includedFaces_;

        //- Excluded faces (per patch) that do not move with the MRF
        labelListList excludedFaces_;

        const dimensionedVector origin_;
        dimensionedVector axis_;
        const dimensionedScalar omega_;
        dimensionedVector Omega_;


    // Private Member Functions

        //- Divide faces in frame according to patch
        void setMRFFaces
        (
            labelList& faceType,
            const labelList& excludedPatchIDs
        );

        //- Divide faces in frame according to patch
        void setMRFFaces();


        //- Disallow default bitwise copy construct
        MRFZone(const MRFZone&);

        //- Disallow default bitwise assignment
        void operator=(const MRFZone&);


public:

    // Constructors

        //- Construct from fvMesh and Istream
        MRFZone(const fvMesh& mesh, Istream& is);

        //- Return clone
        autoPtr<MRFZone> clone() const
        {
            notImplemented("autoPtr<MRFZone> clone() const");
            return autoPtr<MRFZone>(NULL);
        }

        //- Return a pointer to a new MRFZone created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<MRFZone> operator()(Istream& is) const
            {
                return autoPtr<MRFZone>(new MRFZone(mesh_, is));
            }
        };


    // Member Functions

        //- Update the mesh corresponding to given map
        void updateMesh(const mapPolyMesh& mpm)
        {
            // Only updates face addressing
            setMRFFaces();
        }


        // Incompressible MRF

            //- Add the Coriolis force contribution to the momentum equation
            void addCoriolis(fvVectorMatrix& UEqn) const;

            //- Make the given absolute flux relative within the MRF region
            //  Use only for volumetric flux.  HJ, 12/Dec/2009
            void relativeFlux(surfaceScalarField& phi) const;


        // Compressible MRF

            //- Add the Coriolis force contribution to the momentum equation
            //  Use for equations with prefactor fields, eg. density.
            //  This is typically visible in the ddt term
            //  HJ, 12/Dec/2009
            void addCoriolis
            (
                const volScalarField& rho,
                fvVectorMatrix& UEqn
            ) const;

            //- Make the given absolute flux relative within the MRF region
            //  Use for fluxes with prefactor fields, eg. density.
            //  For consistency, requires face-interpolated prefactor
            //  HJ. 12/Dec/2009
            void relativeFlux
            (
                const surfaceScalarField& rhof,
                surfaceScalarField& phi
            ) const;


        //- Correct the boundary velocity for the roation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const MRFZone&)
        {
            notImplemented("Ostream& operator<<(Ostream& os, const MRFZone&)");
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::fineMgFieldLevel

Description
    Fine multigrid field level. Allow for access to fvField at the fine level.

SourceFiles
    fineMgFieldLevel.C

\*---------------------------------------------------------------------------*/

#ifndef fineMgFieldLevel_H
#define fineMgFieldLevel_H

#include "mgFieldLevel.H"
#include "mgMeshLevel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fineMgFieldLevel Declaration
\*---------------------------------------------------------------------------*/

class fineMgFieldLevel
:
    public mgFieldLevel
{
    // Private data

        // reference to fine mesh level
        const mgMeshLevel& meshLevel_;

        // reference to fine mesh
        const fvMesh& mesh_;

        //- Number of cells
        label nCells_;

        //- Number of internal face
        label nInternalFaces_;


    // Primitive variables

        //- Fine p
        const volScalarField& p_;

        //- Fine U
        const volVectorField& U_;

        //- Fine T
        const volScalarField& T_;

        //- Fine Cv
        volScalarField const Cv_;

        //- Fine R
        volScalarField const R_;


    // Residuals

        //- fine rhoRes
        const volScalarField& rhoRes_;

        //- fine rhoURes
        const volVectorField& rhoURes_;

        //- fine rhoERes
        const volScalarField& rhoERes_;


    // Fluxes

        //- fine rhoFlux
        surfaceScalarField& rhoFlux_;

        //- fine rhoUFlux
        surfaceVectorField& rhoUFlux_;

        //- fine rhoEFlux
        surfaceScalarField& rhoEFlux_;


    // Conservative variables

        //- fine rho
        volScalarField& rho_;

        //- fine rhoU
        volVectorField& rhoU_;

        //- fine rhoE
        volScalarField& rhoE_;


    // Helper data

        label level_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineMgFieldLevel(const fineMgFieldLevel&);

        //- Disallow default bitwise assignment
        void operator=(const fineMgFieldLevel&);


public:

        //- Runtime type information
        TypeName("fineMgFieldLevel");


    // Constructors

        //- Construct from mesh
        fineMgFieldLevel
        (
            const mgMeshLevel& fineMeshLevel,
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& T,
            const volScalarField& Cv,
            const volScalarField& R,
            const volScalarField& rhoRes,
            const volVectorField& rhoURes,
            const volScalarField& rhoERes,
            surfaceScalarField& rhoFlux,
            surfaceVectorField& rhoUFlux,
            surfaceScalarField& rhoEFlux,
            volScalarField& rho,
            volVectorField& rhoU,
            volScalarField& rhoE
        )
        :
          meshLevel_(fineMeshLevel),
          mesh_(fineMeshLevel.mesh()),
          nCells_(fineMeshLevel.nCells()),
          nInternalFaces_(fineMeshLevel.nInternalFaces()),
          p_(p),
          U_(U),
          T_(T),
          Cv_(Cv),
          R_(R),
          rhoRes_(rhoRes),
          rhoURes_(rhoURes),
          rhoERes_(rhoERes),
          rhoFlux_(rhoFlux),
          rhoUFlux_(rhoUFlux),
          rhoEFlux_(rhoEFlux),
          rho_(rho),
          rhoU_(rhoU),
          rhoE_(rhoE),
          level_(0)
        {}


    //- Destructor

        virtual ~fineMgFieldLevel()
        {}


    // Member Functions

            //- Is this the finest level?
            virtual bool finest() const
            {
                return meshLevel_.finest();
            }

        // Addressing

            //- Return number of cells
            virtual label nCells() const
            {
                return mesh_.nCells();
            }

            //- Return number of internal faces
            virtual label nInternalFaces() const
            {
                return mesh_.nInternalFaces();
            }

            //- Return number of patches
            virtual label nPatches() const
            {
                return mesh_.boundary().size();
            }

            //- Access to child array
            virtual const labelList& child() const
            {
                return meshLevel_.child();
            }

            //- Access to fine level
            virtual const mgFieldLevel& fineLevel() const
            {
                return *this;
            }


        // Field data


            //- Return p
            virtual const scalarField& pVar() const
            {
                return p_.internalField();
            }

            virtual const volScalarField& p() const
            {
                return p_;
            }

            //- Return p value for patches
            virtual const scalarField& patchP
            (
                const label patchNo
            ) const
            {
                return p_.boundaryField()[patchNo];
            }

            //- Return U
            virtual const vectorField& UVar() const
            {
                return U_.internalField();
            }

            virtual const volVectorField& U() const
            {
                return U_;
            }

            //- Return U value for patches
            virtual const vectorField& patchU
            (
                const label patchNo
            ) const
            {
                return U_.boundaryField()[patchNo];
            }

            //- Return T
            virtual const scalarField& TVar() const
            {
                return T_.internalField();
            }

            virtual const volScalarField& T() const
            {
                return T_;
            }

            //- Return T value for patches
            virtual const scalarField& patchT
            (
                const label patchNo
            ) const
            {
                return T_.boundaryField()[patchNo];
            }

            //- Return Cv
            virtual const scalarField& CvVar() const
            {
                return Cv_.internalField();
            }

            virtual const volScalarField& Cv() const
            {
                return Cv_;
            }

            //- Return p value for patches
            virtual const scalarField& patchCv
            (
                const label patchNo
            ) const
            {
                return Cv_.boundaryField()[patchNo];
            }

            //- Return R
            virtual const scalarField& RVar() const
            {
                return R_.internalField();
            }

            virtual const volScalarField& R() const
            {
                return R_;
            }

            //- Return p value for patches
            virtual const scalarField& patchR
            (
                const label patchNo
            ) const
            {
                return R_.boundaryField()[patchNo];
            }


        // Fluxes

            //- Access to rhoFlux field
            virtual surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Access to rhoUFlux field
            virtual surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Access to rhoEFlux field
            virtual surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }


        // Residuals

            //- Access to rhoRes field
            virtual const volScalarField& rhoRes() const
            {
                return rhoRes_;
            }

            virtual const scalarField& rhoResVar() const
            {
                return rhoRes_.internalField();
            }

            virtual const scalarField& patchRhoRes
            (
                const label patchNo
            ) const
            {
                return rhoRes_.boundaryField()[patchNo];
            }

            //- Access to rhoURes field
            virtual const volVectorField& rhoURes() const
            {
                return rhoURes_;
            }

            virtual const vectorField& rhoUResVar() const
            {
                return rhoURes_.internalField();
            }

            virtual const vectorField& patchRhoURes
            (
                const label patchNo
            ) const
            {
                return rhoURes_.boundaryField()[patchNo];
            }

            //- Access to rhoERes field
            virtual const volScalarField& rhoERes() const
            {
                return rhoERes_;
            }

            virtual const scalarField& rhoEResVar() const
            {
                return rhoERes_.internalField();
            }

            virtual const scalarField& patchRhoERes
            (
                const label patchNo
            ) const
            {
                return rhoERes_.boundaryField()[patchNo];
            }


        // Conservative variables

            //- Access to rho field
            virtual const volScalarField& rho() const
            {
                return rho_;
            }

            virtual const scalarField& rhoVar() const
            {
                return rho_.internalField();
            }

            virtual const scalarField& patchRho
            (
                const label patchNo
            ) const
            {
                return rho_.boundaryField()[patchNo];
            }

            //- Access to rhoU field
            virtual const volVectorField& rhoU() const
            {
                return rhoU_;
            }

            virtual const vectorField& rhoUVar() const
            {
                return rhoU_.internalField();
            }

            virtual const vectorField& patchRhoU
            (
                const label patchNo
            ) const
            {
                return rhoU_.boundaryField()[patchNo];
            }

            //- Access to rhoE field
            virtual const volScalarField& rhoE() const
            {
                return rhoE_;
            }

            virtual const scalarField& rhoEVar() const
            {
                return rhoE_.internalField();
            }

            virtual const scalarField& patchRhoE
            (
                const label patchNo
            ) const
            {
                return rhoE_.boundaryField()[patchNo];
            }

            //- Access to helper data
            virtual label const& level() const
            {
                return level_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

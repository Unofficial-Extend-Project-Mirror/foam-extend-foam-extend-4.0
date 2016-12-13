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
    Foam::coarseMgFieldLevel

Description
    Coarse level of multigrid hierarchy, holding fields

SourceFiles
    coarseMgFieldLevel.C

\*---------------------------------------------------------------------------*/

#ifndef coarseMgFieldLevel_H
#define coarseMgFieldLevel_H

#include "mgFieldLevel.H"
#include "mgMeshLevel.H"
#include "coarseMgMeshLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coarseMgFieldLevel Declaration
\*---------------------------------------------------------------------------*/

class coarseMgFieldLevel
:
    public mgFieldLevel
{
    // Private data

        // Mesh levels

            //- Reference to the current mesh level
            //  corresponding to the current field level
            const mgMeshLevel& meshLevel_;

            //- Reference to fine data level
           const mgFieldLevel& fineFieldLevel_;


        // Geometry

            //- Number of cells
            label nCells_;

            //- Number of internal face
            label nInternalFaces_;


        // Fields

            //- Coarse p
            scalarField p_;

            //- Coarse U
            vectorField U_;

            //- Coarse T
            scalarField T_;


        // Thermodynamics fields

            //- Coarse Cv
            scalarField Cv_;

            //- Coarse R
            scalarField R_;


        // Residuals

            //- Coarse rhoRes
            scalarField rhoRes_;

            //- Coarse rhoURes
            vectorField rhoURes_;

            //- Coarse rhoERes
            scalarField rhoERes_;


        // Fluxes

            //- Fine rhoFlux
            surfaceScalarField& rhoFlux_;

            //- Fine rhoUFlux
            surfaceVectorField& rhoUFlux_;

            //- Fine rhoEFlux
            surfaceScalarField& rhoEFlux_;

            //- Coarse rho
            scalarField rho_;

            //- Coarse rhoU
            vectorField rhoU_;

            //- Coarse rhoE
            scalarField rhoE_;


        //- Helper data
        label level_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coarseMgFieldLevel(const coarseMgFieldLevel&);

        //- Disallow default bitwise assignment
        void operator=(const coarseMgFieldLevel&);

        // Restriction operators used for interpolation
       template <typename T>
       void interpolate
       (
           Field<T>& c,
           Field<T> f
       )
       {
           c = pTraits<T>::zero;
           f *= meshLevel_.fineLevel().cellVolumes();
           fineFieldLevel_.restrict(c, f);
           c /= meshLevel_.cellVolumes();
       }


public:

        //- Runtime type information
        TypeName("coarseMgFieldLevel");


    // Constructors

        //- Construct from fine mgFieldLevel
        coarseMgFieldLevel
        (
            const mgMeshLevel& coarseMeshLevel,
            const mgFieldLevel& fineFieldLevel
        )
        :
            meshLevel_(coarseMeshLevel),
            fineFieldLevel_(fineFieldLevel),
            nCells_(coarseMeshLevel.nCells()),
            nInternalFaces_(coarseMeshLevel.nInternalFaces()),
            p_(scalarField(coarseMeshLevel.nCells())),
            U_(vectorField(coarseMeshLevel.nCells())),
            T_(scalarField(coarseMeshLevel.nCells())),
            Cv_(scalarField(coarseMeshLevel.nCells())),
            R_(scalarField(coarseMeshLevel.nCells())),
            rhoRes_(scalarField(coarseMeshLevel.nCells())),
            rhoURes_(vectorField(coarseMeshLevel.nCells())),
            rhoERes_(scalarField(coarseMeshLevel.nCells())),
            rhoFlux_(fineFieldLevel.rhoFlux()),
            rhoUFlux_(fineFieldLevel.rhoUFlux()),
            rhoEFlux_(fineFieldLevel.rhoEFlux()),
            rho_(scalarField(coarseMeshLevel.nCells())),
            rhoU_(vectorField(coarseMeshLevel.nCells())),
            rhoE_(scalarField(coarseMeshLevel.nCells())),
            level_(fineFieldLevel.level())
        {
            ++level_;

            interpolate(p_, fineFieldLevel_.pVar());
            interpolate(U_, fineFieldLevel_.UVar());
            interpolate(T_, fineFieldLevel_.TVar());

            interpolate(Cv_, fineFieldLevel_.CvVar());
            interpolate(R_, fineFieldLevel_.RVar());

            interpolate(rho_, fineFieldLevel_.rhoVar());
            interpolate(rhoU_, fineFieldLevel_.rhoUVar());
            interpolate(rhoE_, fineFieldLevel_.rhoEVar());

            if (level_ == 1)
            {
                const scalarField& fineVol =
                    meshLevel_.fineLevel().cellVolumes();

                rhoRes_ = 0;
                scalarField fineRhoRes = fineFieldLevel_.rhoResVar();
                fineRhoRes *= fineVol;
                fineFieldLevel_.restrict(rhoRes_, fineRhoRes);

                rhoURes_ = vector::zero;
                vectorField fineRhoURes = fineFieldLevel_.rhoUResVar();
                fineRhoURes *= fineVol;
                fineFieldLevel_.restrict(rhoURes_, fineRhoURes);

                rhoERes_ = 0;
                scalarField fineRhoERes = fineFieldLevel_.rhoEResVar();
                fineRhoERes *= fineVol;
                fineFieldLevel_.restrict(rhoERes_, fineRhoERes);
            }
            else
            {
                rhoRes_ = 0;
                fineFieldLevel_.restrict
                (
                    rhoRes_,
                    fineFieldLevel_.rhoResVar()
                );

                rhoURes_ = vector::zero;
                fineFieldLevel_.restrict
                (
                     rhoURes_,
                     fineFieldLevel_.rhoUResVar()
                );

                rhoERes_ = 0;
                fineFieldLevel_.restrict
                (
                    rhoERes_,
                    fineFieldLevel_.rhoEResVar()
                );
            }
        }


    //- Destructor
    virtual ~coarseMgFieldLevel()
    {}


    // Member Functions

        // Access

            //- Is this the finest level?
            virtual bool finest() const
            {
                return meshLevel_.finest();
            }


        // Addressing

            //- Return number of cells
            virtual label nCells() const
            {
                return nCells_;
            }

            //- Return number of internal faces
            virtual label nInternalFaces() const
            {
                return nInternalFaces_;
            }

            //- Return number of patches
            virtual label nPatches() const
            {
                return fineFieldLevel_.nPatches();
            }

            //- Access to child array
            virtual const labelList& child() const
            {
                return meshLevel_.child();
            }

            //- Access to fine level
            virtual const mgFieldLevel& fineLevel() const
            {
                return fineFieldLevel_;
            }


        // Field data

            //- Return p
            virtual const scalarField& pVar() const
            {
                return p_;
            }

            //- Access to fine level pressure field
            virtual const volScalarField& p() const
            {
                return fineFieldLevel_.p();
            }

            //- Return p value for patches on fine level
            virtual const scalarField& patchP
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchP(patchNo);
            }

             //- Return U
            virtual const vectorField& UVar() const
            {
                return U_;
            }

            virtual const volVectorField& U() const
            {
                return fineFieldLevel_.U();
            }

            //- Return U value for patches
            virtual const vectorField& patchU
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchU(patchNo);
            }

            //- Return T
            virtual const scalarField& TVar() const
            {
                return T_;
            }

            virtual const volScalarField& T() const
            {
                return fineFieldLevel_.T();
            }

            //- Return T value for patches on fine level
            virtual const scalarField& patchT
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchT(patchNo);
            }

            //- Return Cv
            virtual const scalarField& CvVar() const
            {
                return Cv_;
            }

            virtual const volScalarField& Cv() const
            {
                return fineFieldLevel_.Cv();
            }

            //- Return Cv value for patches
            virtual const scalarField& patchCv
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchCv(patchNo);
            }

            //- Return R
            virtual const scalarField& RVar() const
            {
                return R_;
            }

            virtual const volScalarField& R() const
            {
                return fineFieldLevel_.R();
            }

            //- Return R value for patches on fine level
            virtual const scalarField& patchR
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchR(patchNo);
            }


        // Fluxes on fine level

            virtual surfaceScalarField& rhoFlux() const
            {
                return fineFieldLevel_.rhoFlux();
            }

            virtual surfaceVectorField& rhoUFlux() const
            {
                return fineFieldLevel_.rhoUFlux();
            }

            virtual surfaceScalarField& rhoEFlux() const
            {
                return fineFieldLevel_.rhoEFlux();
            }

        // Residuals

           //- Return rhoRes
            virtual const scalarField& rhoResVar() const
            {
                return rhoRes_;
            }

            virtual const volScalarField& rhoRes() const
            {
                return fineFieldLevel_.rhoRes();
            }

            //- Return rhoRes value for patches on fine level
            virtual const scalarField& patchRhoRes
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRhoRes(patchNo);
            }

            //- Return rhoURes
            virtual const vectorField& rhoUResVar() const
            {
                return rhoURes_;
            }

            virtual const volVectorField& rhoURes() const
            {
                return fineFieldLevel_.rhoURes();
            }

            //- Return rhoRes value for patches on fine level
            virtual const vectorField& patchRhoURes
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRhoURes(patchNo);
            }

            //- Return rhoERes
            virtual const scalarField& rhoEResVar() const
            {
                return rhoERes_;
            }

            virtual const volScalarField& rhoERes() const
            {
                return fineFieldLevel_.rhoERes();
            }

            //- Return rhoERes value for patches on fine level
            virtual const scalarField& patchRhoERes
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRhoERes(patchNo);
            }


        // Conservative variables

            //- Return rho
            virtual const scalarField& rhoVar() const
            {
                return rho_;
            }

            virtual const volScalarField& rho() const
            {
                return fineFieldLevel_.rho();
            }

            //- Return rho value for patches on fine level
            virtual const scalarField& patchRho
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRho(patchNo);
            }

            //- Return rhoU
            virtual const vectorField& rhoUVar() const
            {
                return rhoU_;
            }

            virtual const volVectorField& rhoU() const
            {
                return fineFieldLevel_.rhoU();
            }

            //- Return rhoU value for patches on fine level
            virtual const vectorField& patchRhoU
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRhoU(patchNo);
            }

            //- Return rhoE
            virtual const scalarField& rhoEVar() const
            {
                return rhoE_;
            }

            virtual const volScalarField& rhoE() const
            {
                return fineFieldLevel_.rhoE();
            }

            //- Return rho value for patches on fine level
            virtual const scalarField& patchRhoE
            (
                const label patchNo
            ) const
            {
                return fineFieldLevel_.patchRhoE(patchNo);
            }

            //- Access to helper data
            virtual const label& level() const
            {
                return level_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

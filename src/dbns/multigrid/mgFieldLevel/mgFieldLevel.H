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
    Foam::mgFieldLevel

Description
    Virtual base class for a level of multigrid hierarchy.

SourceFiles
    mgFieldLevel.C

\*---------------------------------------------------------------------------*/

#ifndef mgFieldLevel_H
#define mgFieldLevel_H

#include "typeInfo.H"
#include "fvMesh.H"
#include "primitiveFields.H"
#include "coarseMgMeshLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class mgFieldLevel Declaration
\*---------------------------------------------------------------------------*/

class mgFieldLevel
{
public:

        //- Runtime type information
        virtual const word& type() const = 0;

        //- Access to child array
        virtual const labelList& child() const = 0;

        //- Access to fine level
        virtual const mgFieldLevel& fineLevel() const = 0;

        //- Is this the finest level?
        virtual bool finest() const = 0;


    //- Destructor

        virtual ~mgFieldLevel()
        {}


    // Member Functions


        //- Access

            //- Return number of cells
            virtual label nCells() const = 0;

            //- Return number of internal faces
            virtual label nInternalFaces() const = 0;

            //- Return number of patches
            virtual label nPatches() const = 0;


        // Primitive variables

            //- Access to p field
            virtual const scalarField& pVar() const = 0;

            //- Access to U field
            virtual const vectorField& UVar() const = 0;

            //- Access to T field
            virtual const scalarField& TVar() const = 0;

            //- Access to fine p field
            virtual const volScalarField& p() const = 0;

            //- Access to fine U field
            virtual const volVectorField& U() const = 0;

            //- Access to fine T field
            virtual const volScalarField& T() const = 0;

            //- Access to p boundary field
            virtual const scalarField& patchP(const label patchNo) const = 0;

            //- Access to U field
            virtual const vectorField& patchU(const label patchNo) const = 0;

            //- Access to T field
            virtual const scalarField& patchT(const label patchNo) const = 0;

            //- Access to Cv field
            virtual const scalarField& CvVar() const = 0;

            //- Access to fine Cv field
            virtual const volScalarField& Cv() const = 0;

            //- Access to Cv boundary field
            virtual const scalarField& patchCv(const label patchNo) const = 0;

            //- Access to R field
            virtual const scalarField& RVar() const = 0;

            //- Access to fine R field
            virtual const volScalarField& R() const = 0;

            //- Access to R boundary field
            virtual const scalarField& patchR(const label patchNo) const = 0;

        // Fine fluxes

            //- Access to fine rhoFlux field
            virtual surfaceScalarField& rhoFlux() const = 0;

            //- Access to fine rhoUFlux field
            virtual surfaceVectorField& rhoUFlux() const = 0;

            //- Access to fine rhoEFlux field
            virtual surfaceScalarField& rhoEFlux() const = 0;


        // Residuals

            //- Access to rhoRes field
            virtual const scalarField& rhoResVar() const = 0;

            //- Access to fine rhoRes field
            virtual const volScalarField& rhoRes() const = 0;

            //- Access to rhoRes boundary field
            virtual const scalarField& patchRhoRes
            (
                const label patchNo
            ) const = 0;

            //- Access to rhoURes field
            virtual const vectorField& rhoUResVar() const = 0;

            //- Access to fine rhoRes field
            virtual const volVectorField& rhoURes() const = 0;

            //- Access to rhoRes boundary field
            virtual const vectorField& patchRhoURes
            (
                const label patchNo
            ) const = 0;

            //- Access to rhoERes field
            virtual const scalarField& rhoEResVar() const = 0;

            //- Access to fine rhoERes field
            virtual const volScalarField& rhoERes() const = 0;

            //- Access to rhoERes boundary field
            virtual const scalarField& patchRhoERes
            (
                const label patchNo
            ) const = 0;


        // Conservative variables

            //- Access to rho field
            virtual const scalarField& rhoVar() const = 0;

            //- Access to fine rho field
            virtual const volScalarField& rho() const = 0;

            //- Access to rho boundary field
            virtual const scalarField& patchRho
            (
                const label patchNo
            ) const = 0;

            //- Access to rhoURes field
            virtual const vectorField& rhoUVar() const = 0;

            //- Access to fine rhoRes field
            virtual const volVectorField& rhoU() const = 0;

            //- Access to rhoRes boundary field
            virtual const vectorField& patchRhoU
            (
                const label patchNo
            ) const = 0;

            //- Access to rhoERes field
            virtual const scalarField& rhoEVar() const = 0;

            //- Access to fine rhoERes field
            virtual const volScalarField& rhoE() const = 0;

            //- Access to rhoERes boundary field
            virtual const scalarField& patchRhoE
            (
                const label patchNo
            ) const = 0;

            //- Access to rhoERes boundary field
            virtual label const& level() const = 0;


        // Restriction and prolongation

            //- Restrict (from fine to coarse).  Call from the level
            //  where data originates
            template<class T>
            void restrict
            (
                UList<T>& coarseData,
                const UList<T>& fineData
            ) const
            {
                // Get addressing
                const labelList& c = child();

                // Restriction = summation.  Set coarse data to zero
                coarseData = pTraits<T>::zero;

                forAll (c, i)
                {
                    coarseData[c[i]] += fineData[i];
                }
            }

            //- Prolong (from  coarse to fine).  Call from the level
            //  where data originates
            template<class T>
            void prolong
            (
                const UList<T>& coarseData,
                UList<T>& fineData
            ) const
            {
                if (finest())
                {
                    FatalErrorIn
                    (
                        "void mgMeshLevel::prolong\n"
                        "(\n"
                        "    UList<T>& coarseData,\n"
                        "    const UList<T>& fineData\n"
                        ") const\n"
                    )   << "Requested prolong from finest level"
                        << abort(FatalError);
                }

                if
                (
                    coarseData.size() != nCells()
                 || fineData.size() != fineLevel().nCells()
                )
                {
                    FatalErrorIn
                    (
                        "void mgMeshLevel::prolong\n"
                        "(\n"
                        "    UList<T>& coarseData,\n"
                        "    const UList<T>& fineData\n"
                        ") const\n"
                    )    << "Incorrect data sizes "
                         << abort(FatalError);
                }

                // Get addressing
                const labelList& c = fineLevel().child();

                // Prolongation = assignment.  Initialisation not needed
                forAll (c, i)
                {
                    fineData[i] = coarseData[c[i]];
                }
            }

            //- Prolong to the finest level (from  coarse to fine)
            template<class T>
            tmp<Field<T> > prolongToFinest
            (
                const UList<T>& coarseData
            ) const
            {
                if (finest())
                {
                    // Already on finest level
                    return tmp<Field<T> >(new Field<T>(coarseData));
                }
                else if (fineLevel().finest())
                {
                    tmp<Field<T> > tresult
                    (
                        new Field<T>(fineLevel().nCells())
                    );

                    prolong(coarseData, tresult());

                    // Reached finest level, return
                    return tresult;
                }
                else
                {
                    Field<T> intermediate(fineLevel().nCells());

                    prolong(coarseData, intermediate);

                    // Recursive return
                    return fineLevel().prolongToFinest(intermediate);
                }
            }

        // Coarsening

            //- Create next level from current level
            autoPtr<mgFieldLevel> makeNextLevel
            (
                const mgMeshLevel& coarseMeshLevel
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

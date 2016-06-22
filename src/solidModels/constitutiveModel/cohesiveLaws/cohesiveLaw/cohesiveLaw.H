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
    cohesiveLaw

Description
    Cohesive law for solids.
    The cohesive properties are kept separate from the mechanical properties
    which are kept in rheologyLaw.

SourceFiles
    cohesiveLaw.C
    newCohesiveLaw.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef cohesiveLaw_H
#define cohesiveLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class cohesiveLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to stress field
        const volSymmTensorField& sigma_;

    // Private Member Functions

        //- Disallow copy construct
        cohesiveLaw(const cohesiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const cohesiveLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return sigma_.mesh();
        }

        //- Return reference to stress field
        const volSymmTensorField& sigma() const
        {
            return sigma_;
        }


public:

    //- Runtime type information
    TypeName("cohesiveLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cohesiveLaw,
            dictionary,
            (
                const word name,
                const volSymmTensorField& sigma,
                const dictionary& dict
            ),
            (name, sigma, dict)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<cohesiveLaw> New
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        cohesiveLaw
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~cohesiveLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return materials
    virtual tmp<volScalarField> materials() const = 0;

        //-  Return sigmaMax for interface given two material indicators
        virtual scalar interfaceSigmaMax(double mat1, double mat2) const
    {
        return 0;
    }

        //-  Return tauMax for interface given two material indicators
        virtual scalar interfaceTauMax(double mat1, double mat2) const
    {
        return 0;
    }

        //-  Return GIc for interface given two material indicators
        virtual scalar interfaceGIc(double mat1, double mat2) const
    {
        return 0;
    }

        //-  Return GIIc for interface given two material indicators
        virtual scalar interfaceGIIc(double mat1, double mat2) const
    {
        return 0;
    }

        //- Return opening cohesive strength - surface field
        virtual scalar sigmaMaxValue() const = 0;
        virtual tmp<surfaceScalarField> sigmaMax() const = 0;

        //- Return shearing cohesive strength - surface field
        virtual scalar tauMaxValue() const = 0;
        virtual tmp<surfaceScalarField> tauMax() const = 0;

        //- Return Mode I fracture energy - surface field
        virtual scalar GIcValue() const = 0;
        virtual tmp<surfaceScalarField> GIc() const = 0;

        //- Return Mode II fracture energy - surface field
        virtual scalar GIIcValue() const = 0;
        virtual tmp<surfaceScalarField> GIIc() const = 0;

        //- Return tractions during damage given current deltas and Gs
        virtual void damageTractions
    (
     scalar& tN, // current normal traction
     scalar& tS, // current shear traction
     const scalar deltaN, // current deltaN
     const scalar deltaS, // current deltaS
     const scalar GI, // current dissipated GI
     const scalar GII, // current dissipated GII
     const label faceID, // needed for multi-mat
     const scalarField& globalPatchMaterials
     ) const = 0;

    //- Return interface traction
        //- hmmnn this should not be here - maybe put in solidInterface
    virtual tmp<surfaceVectorField> interfaceTraction
    (
     surfaceVectorField n,
     volVectorField U,
     volTensorField gradU,
     volScalarField mu,
     volScalarField lambda
     ) const = 0;

        //- Correct the rheological model
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

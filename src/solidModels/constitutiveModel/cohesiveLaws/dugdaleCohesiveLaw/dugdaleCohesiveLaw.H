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
    dugdaleCohesiveLaw

Description
    Dugdale cohesive law
    The stress stays constant during damage for a given mode-mix
    Stress may vary if mode-mix varies.

SourceFiles
    dugdaleCohesiveLaw.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef dugdaleCohesiveLaw_H
#define dugdaleCohesiveLaw_H

#include "cohesiveLaw.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dugdaleCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class dugdaleCohesiveLaw
:
    public cohesiveLaw
{

  private:

    // Private data

        //- Criticial fracture energy mode I
        dimensionedScalar GIc_;

        //- Critical fracture energy mode II
        dimensionedScalar GIIc_;

        //- Maximum opening cohesive strength
        dimensionedScalar sigmaMax_;

        //- Maximum shearing cohesive strength
        dimensionedScalar tauMax_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dugdaleCohesiveLaw(const dugdaleCohesiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const dugdaleCohesiveLaw&);

public:

    //- Runtime type information
    TypeName("dugdale");

    // Static data members


    // Constructors

        //- Construct from dictionary
        dugdaleCohesiveLaw
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~dugdaleCohesiveLaw();


    // Member Functions

        //- return materials;
    virtual tmp<volScalarField> materials() const;

        //- Return opening cohesive strength - surface field
        virtual scalar sigmaMaxValue() const
        {
      return sigmaMax_.value();
    };
        virtual tmp<surfaceScalarField> sigmaMax() const;

        //- Return shearing cohesive strength - surface field
        virtual scalar tauMaxValue() const
        {
      return tauMax_.value();
    };
        virtual tmp<surfaceScalarField> tauMax() const;

        //-  Return Mode I fracture energy
        virtual scalar GIcValue() const
        {
      return GIc_.value();
    };
        virtual tmp<surfaceScalarField> GIc() const;

        //-  Return Mode II fracture energy
        virtual scalar GIIcValue() const
        {
      return GIIc_.value();
    };
        virtual tmp<surfaceScalarField> GIIc() const;

        //- Return normal traction (during damage) given current deltas
        //virtual scalar normalTraction(scalar deltaN, scalar deltaS) const;

        //- Return shear traction (during damage) given current deltas
        //virtual scalar shearTraction(scalar deltaN, scalar deltaS) const;

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
     ) const;

    //- Calculate tractions
    virtual tmp<surfaceVectorField> interfaceTraction
    (
     surfaceVectorField n,
     volVectorField U,
     volTensorField gradU,
     volScalarField mu,
     volScalarField lambda
     ) const;

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

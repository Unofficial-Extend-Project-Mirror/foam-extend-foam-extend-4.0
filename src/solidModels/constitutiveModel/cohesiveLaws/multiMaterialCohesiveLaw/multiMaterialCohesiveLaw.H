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
    multiMaterialCohesiveLaw

Description
    Zoned multi-material cohesive controlled by an indicator field.
    Essentially multiMaterial class for cohesive properties instead
    of cohesive properties

SourceFiles
    multiMaterialCohesiveLaw.C

Author
    Philip Cardiff UCD
    based on multiMaterial

\*---------------------------------------------------------------------------*/

#ifndef multiMaterialCohesiveLaw_H
#define multiMaterialCohesiveLaw_H

#include "cohesiveLaw.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiMaterialCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class multiMaterialCohesiveLaw
:
    public cohesiveLaw,
    public PtrList<cohesiveLaw>
{

private:

  // Private data

       //- Material indicator field
       // we will read materials again; this will allow
       // the rheologyLaw to be a single material whereas the cohesiveLaw
       // be multiple materials
       //volScalarField materials_;
       const volScalarField& materials_;

        // List of interface cohesive laws
        // number of interfaces is a triangular number
        // (this->.size()*(this->.size() - 1)/2)
        PtrList<cohesiveLaw> interfaceCohesiveLaws_;

        // indicator surface field
        // for faces not on a multi-material interface the value is equal to
        // the ID of the material in the cohesiveLaws list
        // for faces on a multi-material interface the value is equal to
        // the ID of the material in the interfaceCohesiveLaws plus
        // cohesiveLaws.size()
        //surfaceScalarField materialsSurf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterialCohesiveLaw(const multiMaterialCohesiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterialCohesiveLaw&);

        //- Calculate indicator field given index
        tmp<scalarField> indicator(const label index) const;

        //- Calculate indicator field given index
        scalar indicator(const label index, const label cellID) const;

        //- find interface ID (position in intefaceCohesiveLaws list)
        // given two material indices
        scalar interfaceID(const label mat1, const label mat2) const;


public:

    //- Runtime type information
    TypeName("multiMaterial");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiMaterialCohesiveLaw
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMaterialCohesiveLaw();


    // Member Functions

        //- return materials;
    virtual tmp<volScalarField> materials() const;

        //-  Return sigmaMax for interface given two material indicators
        virtual scalar interfaceSigmaMax(double mat1, double mat2) const;

        //-  Return tauMax for interface given two material indicators
        virtual scalar interfaceTauMax(double mat1, double mat2) const;

        //-  Return GIc for interface given two material indicators
        virtual scalar interfaceGIc(double mat1, double mat2) const;

        //-  Return GIIc for interface given two material indicators
        virtual scalar interfaceGIIc(double mat1, double mat2) const;

        //- Return opening cohesive strength - surface field
        virtual scalar sigmaMaxValue() const
        {
      notImplemented(type() + "::sigmaMaxValue()");
      return 0.0;
    };
        virtual tmp<surfaceScalarField> sigmaMax() const;

        //- Return shearing cohesive strength - surface field
        virtual scalar tauMaxValue() const
        {
      notImplemented(type() + "::tauMaxValue()");
      return 0.0;
    };
        virtual tmp<surfaceScalarField> tauMax() const;

        //-  Return Mode I fracture energy
        virtual scalar GIcValue() const
        {
      notImplemented(type() + "::GIcValue()");
      return 0.0;
    };
        virtual tmp<surfaceScalarField> GIc() const;

        //-  Return Mode II fracture energy
        virtual scalar GIIcValue() const
        {
      notImplemented(type() + "::GIIcValue()");
      return 0.0;
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
        //- hmmnn this should not be here - maybe put in solidInterface
    virtual tmp<surfaceVectorField> interfaceTraction
    (
     surfaceVectorField n,
     volVectorField U,
     volTensorField gradU,
     volScalarField mu,
     volScalarField lambda
     ) const;

        //- Correct the rheological model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

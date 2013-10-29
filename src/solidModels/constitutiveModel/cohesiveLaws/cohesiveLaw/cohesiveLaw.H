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

        //- Return normal traction (during damage) given current deltas
        // virtual scalar normalTraction(scalar deltaN, scalar deltaS) const = 0;

        //- Return shear traction (during damage) given current deltas
        // virtual scalar shearTraction(scalar deltaN, scalar deltaS) const = 0;

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

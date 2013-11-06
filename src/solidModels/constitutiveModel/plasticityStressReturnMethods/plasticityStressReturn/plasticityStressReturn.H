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
    plasticityStressReturn

Description
    Method to calcluate the plastic strain increment, returning the
    stress to the yield surface

SourceFiles
    plasticityStressReturn.C
    newPlasticityStressReturn.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef plasticityStressReturn_H
#define plasticityStressReturn_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
//#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class constitutiveModel;

/*---------------------------------------------------------------------------*\
                         Class plasticityStressReturn Declaration
\*---------------------------------------------------------------------------*/

class plasticityStressReturn
{
    // Private data

        //- Name
        const word name_;


    // Private Member Functions

        //- Disallow copy construct
        plasticityStressReturn(const plasticityStressReturn&);

        //- Disallow default bitwise assignment
        void operator=(const plasticityStressReturn&);


  //protected:

public:

    //- Runtime type information
    TypeName("plasticityStressReturn");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            plasticityStressReturn,
            dictionary,
            (
	     const word name,
	     constitutiveModel& constitutiveModel
            ),
            (name, constitutiveModel)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<plasticityStressReturn> New
        (
	 const word& name,
	 constitutiveModel& constitutiveModel
        );


    // Constructors

        //- Construct from dictionary
        plasticityStressReturn
        (
	 const word& name,
	 constitutiveModel& constitutiveModel
	 );


    // Destructor

        virtual ~plasticityStressReturn()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- If plasticity is active
        virtual bool plasticityActive() const
        {
	  return true;
	}

        //- Increment of plastic strain
        virtual const volSymmTensorField& DEpsilonP() const = 0;

        //- Update the yield stress
        //- called at end of time-step
        virtual void updateYieldStress() = 0;

        //- Correct model
        virtual void correct() = 0;

        //- Write model dictionary
        //virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    noneSolidInterface

Description
    solidInterface procedure that does nothing.
    This essentially allows you to turn off the solidInterface
    procedure.

SourceFiles
    noneSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef noneSolidInterface_H
#define noneSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noneSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class noneSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        noneSolidInterface(const noneSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const noneSolidInterface&);

public:

    //- Runtime type information
    TypeName("none");

    // Static data members


    // Constructors

        //- Construct from dictionary
        noneSolidInterface
        (
     const word& name,
     const fvMesh& mesh,
     const constitutiveModel& rheology
        );


    // Destructor

        virtual ~noneSolidInterface();


    // Member Functions

        // Correct UEqn does nothing
        virtual void correct(fvVectorMatrix& UEqn)  {}

        // Modify properties does nothing
        virtual void modifyProperties
    (
     surfaceScalarField& s
     ) const {}

        virtual void modifyProperties
    (
     surfaceSymmTensor4thOrderField& st
     ) const {}

        virtual void modifyProperties
    (
     surfaceDiagTensorField& dt
     ) const {}

        virtual void modifyProperties
    (
     surfaceScalarField& mu,
     surfaceScalarField& lambda
     ) const {}

        virtual void modifyProperties
    (
     surfaceScalarField& mu,
     surfaceScalarField& lambda,
     surfaceScalarField& threeKalpha
     ) const {}

        virtual void modifyProperties
    (
     surfaceSymmTensor4thOrderField& C,
     surfaceDiagTensorField& K
     ) const {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

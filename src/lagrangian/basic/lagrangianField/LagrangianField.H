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
    Foam::LagrangianField

Description

SourceFiles
    LagrangianFieldI.H
    LagrangianField.C
    LagrangianFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianField_H
#define LagrangianField_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LagrangianField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Cloud>
class LagrangianField
:
    public DimensionedField<Type>
{
    // Private data

        //- Reference to lagrangian cloud this field is part of
        const Cloud& cloud_;


public:

    //- Runtime type information
    TypeName("LagrangianField");


    // Constructors

        //- Construct from components
        LagrangianField
        (
            const IOobject& ioObject,
            const Cloud& cloud,
            const dimensions& dims,
            const Field<Type> field
        );

        //- Construct from IOobject
        LagrangianField
        (
            const IOobject& ioObject,
            const Cloud& cloud
        );

        //- Construct from Cloud data
        LagrangianField(const Cloud& cloud, const int typeFieldIndex);

        //- Construct as copy
        LagrangianField(const LagrangianField<Type>&);


    // Destructor

        ~LagrangianField();


    // Member Functions

        // Access

            //- Map field back into Cloud
            void map(Cloud& cloud, const int typeFieldIndex);


        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const LagrangianField&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const LagrangianField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

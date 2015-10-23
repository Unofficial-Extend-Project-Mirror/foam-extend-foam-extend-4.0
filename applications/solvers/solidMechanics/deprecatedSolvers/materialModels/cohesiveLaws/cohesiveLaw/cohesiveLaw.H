/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Virtual base class for cohesive zone model.

SourceFiles
    cohesiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveLaw_H
#define cohesiveLaw_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dimensionedTypes.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class cohesiveLaw
:
    public refCount
{

private:

    // Private data

        //- Cohesive law coefficients
        dictionary cohesiveLawCoeffs_;

        //- Fracture energy
        dimensionedScalar GIc_;

        //- Maximum cohesive strength
        dimensionedScalar sigmaMax_;


public:

    //- Runtime type information
    TypeName("cohesiveLaw");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            cohesiveLaw,
            dictionary,
            (
                const word& cohesiveLawName,
                const dictionary& dict
            ),
            (cohesiveLawName, dict)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<cohesiveLaw> New
        (
            const word& cohesiveLawName,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        cohesiveLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict
        );

        //- Construct as copy
        cohesiveLaw(const cohesiveLaw&);


        //- Construct and return a clone
        virtual autoPtr<cohesiveLaw> clone() const = 0;


    // Destructor

        virtual ~cohesiveLaw();


    // Member Functions

        //- Return cohesive law coefficients
        const dictionary& cohesiveLawCoeffs() const
        {
            return cohesiveLawCoeffs_;
        }

        //- Return reference to fracture energy
        const dimensionedScalar& GIc() const
        {
            return GIc_;
        }

        //- Return reference to maximum cohesive strength
        const dimensionedScalar& sigmaMax() const
        {
            return sigmaMax_;
        }

        //- Return reference to critical separation distance
        virtual const dimensionedScalar& deltaC() const = 0;

        //- Return current holding traction
        virtual scalar traction(scalar delta) const = 0;

        //- Write dictionary
        void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::ersViewFactor

Description
    helper class to sum up external radiation sources


SourceFiles
    ersViewFactor.C

\*---------------------------------------------------------------------------*/

#ifndef ersViewFactor_H
#define ersViewFactor_H

#include "externalRadiationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ersViewFactor
:
    public externalRadiationSource
{
    // Private data

        //- temperature at infinity
        scalar Tinf_;

        //- view factor
        scalarField F_;

        //- espison
        scalar epsilon_;

public:

    //- Runtime type information
    TypeName("viewFactor");


    // Constructors

        //- Construct from components
        explicit ersViewFactor
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );


        //- Construct from components
        explicit ersViewFactor
        (
            const word& name,
            const dictionary& dict
        );


        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            return autoPtr<externalRadiationSource>
            (
                new ersViewFactor(*this)
            );
        }

    // Destructor

        virtual ~ersViewFactor()
        {};


    // Member functions

        // Access

        const scalarField& F() const
        {
            return F_;
        };

        scalarField& F()
        {
            return F_;
        };

        const scalar& epsilon() const
        {
            return epsilon_;
        };

        scalar& epsilson()
        {
            return epsilon_;
        };


        // Mapping functions

        // Evaluation functions

            //- Return heat flux
            virtual tmp<scalarField> q(const scalarField& Tw) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

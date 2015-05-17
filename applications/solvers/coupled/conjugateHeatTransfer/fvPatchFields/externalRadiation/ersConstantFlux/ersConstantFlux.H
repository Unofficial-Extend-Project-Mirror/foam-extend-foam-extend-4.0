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
    Foam::ersConstantFlux

Description
    helper class to sum up external radiation sources


SourceFiles
    ersConstantFlux.C

\*---------------------------------------------------------------------------*/

#ifndef ersConstantFlux_H
#define ersConstantFlux_H

#include "externalRadiationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ersConstantFlux
:
    public externalRadiationSource
{
    // Private data

        //- heat flux
        scalarField q_;


public:

    //- Runtime type information
    TypeName("constantFlux");


    // Constructors

        //- Construct from components
        explicit ersConstantFlux
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );

        explicit ersConstantFlux
        (
            const word& name
        );

        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            return autoPtr<externalRadiationSource>
            (
                new ersConstantFlux(*this)
            );
        }

    // Destructor

        virtual ~ersConstantFlux()
        {};


    // Member functions

        // Access

        const scalarField& q() const
        {
            return q_;
        };

        scalarField& q()
        {
            return q_;
        };


        // Mapping functions

        // Evaluation functions

            //- Return heat flux
            virtual tmp<scalarField> q(const scalarField&) const
            {
                return q_;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

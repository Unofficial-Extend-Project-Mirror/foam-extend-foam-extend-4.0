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
    Foam::externalRadiationSource

Description
    Virtual base class for external radiation sources


SourceFiles
    externalRadiationSource.C

\*---------------------------------------------------------------------------*/

#ifndef externalRadiationSource_H
#define externalRadiationSource_H

#include "scalarField.H"
#include "tmp.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class externalRadiationSource
{
    // Private data

        //- name
        word name_;


public:

    //- Runtime type information
    TypeName("externalRadiationSource");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            externalRadiationSource,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const fvPatch& p
            ),
            (name, dict, p)
        );


    // Destructor

        virtual ~externalRadiationSource()
        {};

    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<externalRadiationSource> New
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );


        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            notImplemented("autoPt<externalRadiationSource> clone() const");
            return autoPtr<externalRadiationSource>(NULL);

            return autoPtr<externalRadiationSource>
            (
                new externalRadiationSource(*this)
            );
        };

    // Constructors

        //- Construct from components
        explicit externalRadiationSource
        (
            const word& name
        );


    // Member functions

        // Access

        const word& name() const
        {
            return name_;
        }

        // Mapping functions

        // Evaluation functions

            //- Return heat flux
            virtual tmp<scalarField> q(const scalarField& Tw) const
            {
                notImplemented
                (
                    type() +"::q(const scalarField& Tw) const"
                );

                return tmp<scalarField>(NULL);
            };


        //- Write

            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

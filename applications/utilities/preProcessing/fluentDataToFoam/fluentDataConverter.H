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
    fluentDataConverter

Description
    Fluent to OpenFOAM field data converter

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    fluentDataConverter.C

\*---------------------------------------------------------------------------*/

#ifndef fluentDataConverter_H
#define fluentDataConverter_H

#include "fvCFD.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluentDataConverter Declaration
\*---------------------------------------------------------------------------*/

class fluentDataConverter
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Zone to patch name lookup
        const wordIOList zoneToPatchName_;

        //- Field ID list
        const SLList<label>& fieldID_;

        //- Zone ID list
        const SLList<label>& zoneID_;

        //- Zone start index
        const SLList<label>& firstID_;

        //- Zone end index
        const SLList<label>& lastID_;

        //- Fluent zone data
        const SLPtrList<FieldField<Field, scalar> >& zoneData_;
        


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluentDataConverter(const fluentDataConverter&);

        //- Disallow default bitwise assignment
        void operator=(const fluentDataConverter&);


public:

    // Constructors

        //- Construct from components
        fluentDataConverter
        (
            const fvMesh& mesh,
            const SLList<label>& fieldID,
            const SLList<label>& zoneID,
            const SLList<label>& firstID,
            const SLList<label>& lastID,
            const SLPtrList<FieldField<Field, scalar> >& zoneData
        );


    // Destructor - default


    // Member Functions

        //- Convert field
        tmp<volScalarField> convertField
        (
            const word& fieldName,
            const label unitNumber,
            const dimensionedScalar& defaultValue
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

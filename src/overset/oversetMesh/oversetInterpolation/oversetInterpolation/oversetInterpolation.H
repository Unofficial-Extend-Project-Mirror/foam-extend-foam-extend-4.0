/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    oversetInterpolation

Description
    Virtual base class for overset interpolation selection.

Notes
    Weights are organized in a list of fields of scalar fields. First list
    corresponds to particular overset region. Second field corresponds to given
    acceptor cell in region. Third field corresponds to all donors for a given
    acceptor cell. For information on usage, see oversetMesh::interpolate member
    function.

Authors
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

SourceFiles
    oversetInterpolation.C
    newOversetInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef oversetInterpolation_H
#define oversetInterpolation_H

#include "dictionary.H"
#include "typeInfo.H"
#include "labelList.H"
#include "scalarField.H"
#include "vectorField.H"
#include "scalarFieldField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;

/*---------------------------------------------------------------------------*\
                    Class oversetInterpolation Declaration
\*---------------------------------------------------------------------------*/

class oversetInterpolation
{
public:

    // Public data types

        typedef List<scalarFieldField> scalarFieldFieldList;


private:

    // Private data

        //- Overset mesh reference
        const oversetMesh& overset_;


    // Private Member Functions

        //- Disallow copy construct
        oversetInterpolation(const oversetInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const oversetInterpolation&);


public:

    //- Runtime type information
    TypeName("oversetInterpolation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oversetInterpolation,
            word,
            (
                const oversetMesh& overset,
                const word& name
            ),
            (overset, name)
        );


    // Selectors

        //- Return a reference to the selected overset interpolation scheme
        static autoPtr<oversetInterpolation> New
        (
            const oversetMesh& overset,
            const word& name
        );


    // Constructors

        //- Construct from dictionary
        oversetInterpolation
        (
            const oversetMesh& overset,
            const word& name
        );


    //- Destructor
    virtual ~oversetInterpolation()
    {}


    // Member Functions

        // Access

            //- Return reference to overset mesh
            const oversetMesh& overset() const
            {
                return overset_;
            }

            //- Return weights
            virtual const scalarFieldFieldList& weights() const = 0;


        //- Update the interpolation
        virtual void update() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

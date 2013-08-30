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
    Foam::equationVariable

Description
    A generic interface template for active equation variables for the
    equationReader. Use this for objects that act as equation data sources,
    but do not store their data in a permanent location.

SourceFiles
    equationVariableI.H
    equationVariable.C
    equationVariableTemplates.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationVariable_H
#define equationVariable_H

#include "word.H"
#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class equationVariable Declaration
\*---------------------------------------------------------------------------*/

class equationVariable
{
        
public:

    // Destructor
    virtual ~equationVariable() {}

    virtual const word& name() const = 0;
    
    virtual const dimensionSet& dimensions() const = 0;
    
    virtual label lookupComponentIndex(const word) const = 0;
    
    virtual scalar evaluateScalar
    (
        const label componentIndex,
        const label cellIndex,
        const label geoIndex
    ) const = 0;

    virtual void evaluateScalarField
    (
        scalarField& result,
        const label componentIndex,
        const label geoIndex
    ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

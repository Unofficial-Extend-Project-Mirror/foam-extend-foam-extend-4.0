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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh>
void Foam::equationReader::addDataSource
(
    const DimensionedField<scalar, GeoMesh>& dfield
)
{
//    const scalarList * slist(&dfield);
    word name(dfield.name());
    dimensionSet dims(dfield.dimensions());
    addDataSource(dfield, name, dims);
}


template<class GeoMesh>
void Foam::equationReader::evaluateField
(
    const label& index,
    DimensionedField<scalar, GeoMesh>& dfield
)
{
    evaluateField(index, dfield, dfield.dimensions());
}


template<class GeoMesh>
void Foam::equationReader::linkOutput
(
    const word& eqnName,
    DimensionedField<scalar, GeoMesh>& dfield
)
{
    label index(lookup(eqnName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::linkOutput")
            << "Equation name " << eqnName << "not found."
            << abort(FatalError);
    }
    linkOutput
    (
        index,
        dfield,
        dfield.dimensions()
    );
}


template<class GeoMesh>
void Foam::equationReader::linkOutput
(
    label index,
    DimensionedField<scalar, GeoMesh>& dfield
)
{
    linkOutput
    (
        index,
        dfield,
        dfield.dimensions()
    );
}


template<class GeoMesh>
void Foam::equationReader::readEquation
(
    equation eqn,
    DimensionedField<scalar, GeoMesh>& dfield
)
{
    readEquation
    (
        eqn,
        dfield,
        dfield.dimensions()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

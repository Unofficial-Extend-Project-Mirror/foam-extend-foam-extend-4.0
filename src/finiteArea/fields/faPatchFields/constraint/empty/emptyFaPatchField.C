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

#include "emptyFaPatchField.H"
#include "faPatchFieldMapper.H"
#include "areaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>&,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper&
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{
    if (!isType<emptyFaPatch>(p))
    {
        FatalErrorIn
        (
            "emptyFaPatchField<Type>::emptyFaPatchField\n"
            "(\n"
            "    const emptyFaPatchField<Type>&,\n"
            "    const faPatch& p,\n"
            "    const DimensionedField<Type, areaMesh>& iF,\n"
            "    const faPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{
    if (typeid(p) != typeid(emptyFaPatch))
    {
        FatalIOErrorIn
        (
            "emptyFaPatchField<Type>::emptyFaPatchField\n"
            "(\n"
            "    const faPatch& p,\n"
            "    const DimensionedField<Type, areaMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )    << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>& ptf
)
:
    faPatchField<Type>
    (
        ptf.patch(),
        ptf.dimensionedInternalField(),
        Field<Type>(0)
    )
{}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf.patch(), iF, Field<Type>(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void emptyFaPatchField<Type>::updateCoeffs()
{
    // ZT, 26/06/2010, bug-fix for faMesh of zero size
    if (this->dimensionedInternalField().mesh().nFaces())
    {
        if
        (
            this->patch().faPatch::size()
          % this->dimensionedInternalField().mesh().nFaces()
        )
        {
            FatalErrorIn("emptyFaPatchField<Type>::updateCoeffs()")
                << "This mesh contains patches of type empty but is not 1D or 2D\n"
                "    by virtue of the fact that the number of faces of this\n"
                "    empty patch is not divisible by the number of cells."
                    << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

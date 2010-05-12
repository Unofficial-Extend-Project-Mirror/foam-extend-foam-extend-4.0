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

#include "fixedGradientFaPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF),
    gradient_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const fixedGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper)
{}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF),
    gradient_("gradient", dict, p.size())
{
    evaluate();
}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const fixedGradientFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf),
    gradient_(ptf.gradient_)
{}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const fixedGradientFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fixedGradientFaPatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    gradient_.autoMap(m);
}


template<class Type>
void fixedGradientFaPatchField<Type>::rmap
(
    const faPatchField<Type>& ptf,
    const labelList& addr
)
{
    faPatchField<Type>::rmap(ptf, addr);

    const fixedGradientFaPatchField<Type>& fgptf =
        refCast<const fixedGradientFaPatchField<Type> >(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}


template<class Type>
void fixedGradientFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    faPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >(new Field<Type>(this->size(), pTraits<Type>::one));
}


template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> >
fixedGradientFaPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> >
fixedGradientFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradient();
}


template<class Type>
void fixedGradientFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

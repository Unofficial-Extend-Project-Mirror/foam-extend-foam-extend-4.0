/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::timeVaryingUniformFixedValueFaPatchField

Description
    A time-varying form of a uniform fixed value finite area
    boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingUniformFixedValue;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

Note
    This class is derived directly from a fixedValue patch rather than from
    a uniformFixedValue patch.

See Also
    Foam::interpolationTable and Foam::fixedValueFaPatchField

SourceFiles
    timeVaryingUniformFixedValueFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingUniformFixedValueFaPatchField_H
#define timeVaryingUniformFixedValueFaPatchField_H

#include "fixedValueFaPatchField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class timeVaryingUniformFixedValueFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingUniformFixedValueFaPatchField
:
    public fixedValueFaPatchField<Type>
{
    // Private data

        //- The time series being used, including the bounding treatment
        interpolationTable<Type> timeSeries_;


public:

    //- Runtime type information
    TypeName("timeVaryingUniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingUniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingUniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >
            (
                new timeVaryingUniformFixedValueFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type> >
            (
                new timeVaryingUniformFixedValueFaPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the time series used
            const interpolationTable<Type>& timeSeries() const
            {
                return timeSeries_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "timeVaryingUniformFixedValueFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Vector-matrix multiplication operations for a block matrix

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::BlockLduMatrix<Type>::H(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    // Create result
    tmp<Field<Type> > tresult
    (
        new Field<Type>(lduAddr().size(), pTraits<Type>::zero)
    );
    Field<Type>& result = tresult();

    if (this->thereIsUpper())
    {
        const unallocLabelList& l = lduAddr().lowerAddr();
        const unallocLabelList& u = lduAddr().upperAddr();
        const TypeCoeffField& Upper = this->upper();

        // Create multiplication function object
        typename BlockCoeff<Type>::multiply mult;

        // Lower multiplication

        if (symmetric())
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& activeUpper = Upper.asScalar();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(activeUpper[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& activeUpper = Upper.asLinear();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(activeUpper[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& activeUpper = Upper.asSquare();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // Use transpose upper coefficient
                    result[u[coeffI]] -=
                        mult(activeUpper[coeffI].T(), x[l[coeffI]]);
                }
            }
        }
        else // Asymmetric matrix
        {
            const TypeCoeffField& Lower = this->lower();

            if (Lower.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& activeLower = Lower.asScalar();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& activeLower = Lower.asLinear();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& activeLower = Lower.asSquare();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
        }


        // Upper multiplication

        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(activeUpper[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(activeUpper[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(activeUpper[coeffI], x[u[coeffI]]);
            }
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::BlockLduMatrix<Type>::faceH(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const unallocLabelList& u = lduAddr().upperAddr();
    const unallocLabelList& l = lduAddr().lowerAddr();

    // Create result
    tmp<Field<Type> > tresult(new Field<Type>(u.size(), pTraits<Type>::zero));
    Field<Type>& result = tresult();

    if (this->thereIsUpper())
    {
        const TypeCoeffField& Upper = this->upper();

        // Create multiplication function object
        typename BlockCoeff<Type>::multiply mult;

        // Lower multiplication

        if (symmetric())
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& activeUpper = Upper.asScalar();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // This can be optimised with a subtraction.
                    // Currently not done for clarity.  HJ, 31/Oct/2007
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeUpper[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& activeUpper = Upper.asLinear();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // This can be optimised with a subtraction.
                    // Currently not done for clarity.  HJ, 31/Oct/2007
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeUpper[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& activeUpper = Upper.asSquare();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // Use transpose upper coefficient
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeUpper[coeffI].T(), x[l[coeffI]]);
                }
            }
        }
        else // Asymmetric matrix
        {
            const TypeCoeffField& Lower = this->lower();

            if (Lower.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& activeUpper = Upper.asScalar();
                const scalarTypeField& activeLower = Lower.asScalar();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& activeUpper = Upper.asLinear();
                const linearTypeField& activeLower = Lower.asLinear();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& activeUpper = Upper.asSquare();
                const squareTypeField& activeLower = Lower.asSquare();

                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(activeUpper[coeffI], x[u[coeffI]])
                    - mult(activeLower[coeffI], x[l[coeffI]]);
                }
            }
        }
    }

    return tresult;
}


// ************************************************************************* //

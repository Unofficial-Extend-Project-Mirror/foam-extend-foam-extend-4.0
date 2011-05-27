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

Description
    RBF interpolation templates

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::RBFInterpolation::interpolate
(
    const Field<Type>& ctrlField
) const
{
    // HJ and FB (05 Jan 2009)
    // Collect the values from ALL control points to all CPUs
    // Then, each CPU will do interpolation only on local allPoints_

    if (ctrlField.size() != controlPoints_.size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > RBFInterpolation::interpolate\n"
            "(\n"
            "    const Field<Type>& ctrlField\n"
            ") const"
        )   << "Incorrect size of source field.  Size = " << ctrlField.size()
            << " nControlPoints = " << controlPoints_.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(allPoints_.size(), pTraits<Type>::zero)
    );

    Field<Type>& result = tresult();

    // FB 21-12-2008
    // 1) Calculate alpha and beta coefficients using the Inverse
    // 2) Calculate displacements of internal nodes using RBF values,
    //    alpha's and beta's
    // 3) Return displacements using tresult()

    const label nControlPoints = controlPoints_.size();
    const scalarSquareMatrix& mat = this->B();

    // Determine interpolation coefficients
    Field<Type> alpha(nControlPoints, pTraits<Type>::zero);
    Field<Type> beta(4, pTraits<Type>::zero);

    for (label row = 0; row < nControlPoints; row++)
    {
        for (label col = 0; col < nControlPoints; col++)
        {
            alpha[row] += mat[row][col]*ctrlField[col];
        }
    }

    if (polynomials_)
    {
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 4;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                beta[row - nControlPoints] += mat[row][col]*ctrlField[col];
            }
        }
    }

    // Evaluation
    scalar t;

    // Algorithmic improvement, Matteo Lombardi.  21/Mar/2011

    forAll (allPoints_, flPoint)
    {
        // Cut-off function to justify neglecting outer boundary points
        t = (Foam::mag(allPoints_[flPoint] - focalPoint_) - innerRadius_)/
            (outerRadius_ - innerRadius_);

        if (t >= 1)
        {
            // Increment is zero: w = 0
            result[flPoint] = 0*result[flPoint];
        }
        else
        {
            // Full calculation of weights
            scalarField weights =
                RBF_->weights(controlPoints_, allPoints_[flPoint]);

            forAll (controlPoints_, i)
            {
                result[flPoint] += weights[i]*alpha[i];
            }

            if (polynomials_)
            {
                result[flPoint] +=
                    beta[0]
                  + beta[1]*allPoints_[flPoint].x()
                  + beta[2]*allPoints_[flPoint].y()
                  + beta[3]*allPoints_[flPoint].z();
            }

            scalar w;

            if (t <= 0)
            {
                w = 1.0;
            }
            else
            {
                w = 1 - sqr(t)*(3-2*t);
            }

            result[flPoint] = w*result[flPoint];
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


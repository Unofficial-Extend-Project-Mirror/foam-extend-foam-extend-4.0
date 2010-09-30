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
    EigenSolver

\*---------------------------------------------------------------------------*/

#include "EigenSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::EigenSolver<T>::checkMatrix(const SquareMatrix<T>& mtx) const
{
    // Check symmetry
    T assymetry = pTraits<T>::zero;

    for (label i = 0; i < mtx.m(); i++)
    {
        for (label j = i; j < mtx.n(); j++)
        {
            assymetry += mag(mtx[i][j] - mtx[j][i]);
        }
    }

    if (assymetry > mtx.m()*SMALL)
    {
        FatalErrorIn
        (
            "void EigenSolver<T>::checkMatrix(const SquareMatrix<T>& m) const"
        )   << "Matrix is not symmetric.  Assymetry = " << assymetry
            << abort(FatalError);
    }
}


template<class T>
void Foam::EigenSolver<T>::factorise(const SquareMatrix<T>& mtx)
{
    // Copy the original matrix into scratch space
    SquareMatrix<T> a = mtx;

    label n = a.m();

    // Create and initialise the matrix to hold eigen-vectors in columns
    SquareMatrix<T> v(n, pTraits<T>::zero);

    for (label ip = 0; ip < n; ip++)
    {
        v[ip][ip] = pTraits<T>::one;
    }

    // Create and initialise scratch vectors
    List<T> b(n);
    List<T>& d = values_;

    // Initialise b and d to the diagonal of a
    for (label ip = 0; ip < n; ip++)
    {
        b[ip] = a[ip][ip];
        d[ip] = a[ip][ip];
    }

    List<T> z(n, pTraits<T>::zero);

    label nRotations = 0;

    // Main iteration loop
    for (label i = 0; i <= 50; i++)
    {
        T sm = pTraits<T>::zero;
        for (label ip = 0; ip < n - 1; ip++)
        {
            for (label iq = ip + 1; iq < n; iq++)
            {
                sm += mag(a[ip][iq]);
            }
        }

        if (sm < VSMALL)
        {
            // Normal return, relying on quadratic convergence

            // Copy eigen-vectors.  Note that v stores vectors in columns
            for (label ip = 0; ip < n; ip++)
            {
                for (label iq = 0; iq < n; iq++)
                {
                    vectors_[iq][ip] = v[ip][iq];
                }
            }

            return;
        }

        // Calculate treshold
        T treshold = pTraits<T>::zero;

        if (i < 4)
        {
            // Operation on first 3 sweeps
            treshold = 0.2*sm/(n*n);
        }

        T theta, tau, g, h, t, c, s;

        for (label ip = 0; ip < n - 1; ip++)
        {
            for (label iq = ip + 1; iq < n; iq++)
            {
                g = 100.0*mag(a[ip][iq]);

                // After four sweeps, skip the rotation if the off-diagonal
                // element is small
                if
                (
                    i > 4
                 && mag(d[ip] + g) == mag(d[ip])
                 && mag(d[iq] + g) == mag(d[iq])
                )
                {
                    a[ip][iq] = pTraits<T>::zero;
                }
                else if (mag(a[ip][iq]) > treshold)
                {
                    h = d[iq] - d[ip];

                    if (mag(h) + g == mag(h))
                    {
                        t = a[ip][iq]/h;
                    }
                    else
                    {
                        theta = 0.5*h/(a[ip][iq]);

                        t = 1.0/(mag(theta) + Foam::sqrt(1.0 + sqr(theta)));

                        if (theta < 0)
                        {
                            t = -t;
                        }
                    }

                    c = 1.0/sqrt(1.0 + sqr(t));
                    s = t*c;
                    tau = s/(1.0 + c);
                    h = t*a[ip][iq];

                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;

                    a[ip][iq] = pTraits<T>::zero;

                    // Rotations 0 <= j < p
                    for (label j = 0; j < ip; j++)
                    {
                        rotate(a, s, tau, j, ip, j, iq);
                    }

                    // Rotations p < j < q
                    for (label j = ip + 1; j < iq; j++)
                    {
                        rotate(a, s, tau, ip, j, j, iq);
                    }

                    // Rotations q < j < n
                    for (label j = iq + 1; j < n; j++)
                    {
                        rotate(a, s, tau, ip, j, iq, j);
                    }

                    for (label j = 0; j < n; j++)
                    {
                        rotate(v, s, tau, j, ip, j, iq);
                    }

                    // Increment number of rotations
                    nRotations++;
                }
            }
        }

        // Update d with the sum of t a_pq and clear z
        for (label ip = 0; ip < n; ip++)
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = pTraits<T>::zero;
        }

    } // End of main iteration loop

    FatalErrorIn
    (
        "void EigenSolver<T>::factorise(const SquareMatrix<T>& mtx)"
    )
        << "Maximum number of iterations exceeded"
        << abort(FatalError);
}


template<class T>
inline void Foam::EigenSolver<T>::rotate
(
    SquareMatrix<T>& a,
    const T s,
    const T tau,
    const label i,
    const label j,
    const label k,
    const label l
) const
{
    T g = a[i][j];
    T h = a[k][l];

    a[i][j] = g - s*(h + g*tau);
    a[k][l] = h + s*(g - h*tau);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class T>
Foam::EigenSolver<T>::EigenSolver(const SquareMatrix<T>& mtx)
:
    values_(mtx.m(), pTraits<T>::zero),
    vectors_(mtx.m())
{
    // Size the eigen vectors
    forAll (vectors_, rowI)
    {
        vectors_[rowI].setSize(mtx.m());
        vectors_[rowI] = pTraits<T>::zero;
    }

    this->checkMatrix(mtx);

    factorise(mtx);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::EigenSolver<T>::eigenValue(const label n) const
{
    return values_[n];
}


template<class T>
const Foam::List<T>& Foam::EigenSolver<T>::eigenVector(const label n) const
{
    return vectors_[n];
}


// ************************************************************************* //

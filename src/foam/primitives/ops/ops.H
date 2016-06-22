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

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.

    The information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

Note:
    Format of this file has considerably changed to remove cpp pre-processor
    definition in order to help Swig with parsing in.  Implemented by
    Alexey Petrov; merged by HJ, 14/Oct/2010

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T1, class T2>
class eqOp2
{ public: void operator()(T1& x, const T2& y) const { x = y; } };

template<class T>
class eqOp
{ public: void operator()(T& x, const T& y) const { x = y; } };

template<class T1, class T2>
class plusEqOp2
{ public: void operator()(T1& x, const T2& y) const { x += y; } };

template<class T>
 class plusEqOp
{ public: void operator()(T& x, const T& y) const { x += y; } };

template<class T1, class T2>
 class minusEqOp2
{ public: void operator()(T1& x, const T2& y) const { x -= y; } };

template<class T>
class minusEqOp
{ public: void operator()(T& x, const T& y) const { x -= y; } };

template<class T1, class T2>
class multiplyEqOp2
{ public: void operator()(T1& x, const T2& y) const { x *= y; } };

template<class T>
class multiplyEqOp
{ public: void operator()(T& x, const T& y) const { x *= y; } };

template<class T1, class T2>
class divideEqOp2
{ public: void operator()(T1& x, const T2& y) const { x /= y; } };

template<class T>
class divideEqOp
{ public: void operator()(T& x, const T& y) const { x /= y; } };

template<class T1, class T2>
class eqMagOp2
{ public: void operator()(T1& x, const T2& y) const { x = mag(y); } };

template<class T>
class eqMagOp
{ public: void operator()(T& x, const T& y) const { x = mag(y); } };

template<class T>
class eqSignOp
{ public: void operator()(T& x, const T& y) const { x = sign(y); } };

template<class T1, class T2>
class plusEqMagSqrOp2
{ public: void operator()(T1& x, const T2& y) const { x += magSqr(y); } };

template<class T>
class plusEqMagSqrOp
{ public: void operator()(T& x, const T& y) const { x += magSqr(y); } };

template<class T1, class T2>
class maxEqOp2
{ public: void operator()(T1& x, const T2& y) const { x = max(x, y); } };

template<class T>
class maxEqOp
{ public: void operator()(T& x, const T& y) const { x = max(x, y); } };

template<class T1, class T2>
class minEqOp2
{ public: void operator()(T1& x, const T2& y) const { x = min(x, y); } };

template<class T>
class minEqOp
{ public: void operator()(T& x, const T& y) const { x = min(x, y); } };

template<class T1, class T2>
class andEqOp2
{ public: void operator()(T1& x, const T2& y) const { x = (x && y); } };

template<class T>
class andEqOp
{ public: void operator()(T& x, const T& y) const { x = (x && y); } };

template<class T1, class T2>
class orEqOp2
{ public: void operator()(T1& x, const T2& y) const { x = (x || y); } };

template<class T>
class orEqOp
{ public: void operator()(T& x, const T& y) const { x = (x || y); } };

template<class T1, class T2>
class eqMinusOp2
{ public: void operator()(T1& x, const T2& y) const { x = -y; } };

template<class T>
class eqMinusOp
{ public: void operator()(T& x, const T& y) const { x = -y; } };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class T1, class T2>
class sumOp3
{ public: T operator()(const T1& x, const T2& y) const { return x + y; } };

template<class T1, class T2>
class sumOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x + y; } };

template<class T>
class sumOp
{ public: T operator()(const T& x, const T& y) const { return x + y; } };

template<class T, class T1, class T2>
class plusOp3
{ public: T operator()(const T1& x, const T2& y) const { return x + y; } };

template<class T1, class T2>
class plusOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x + y; } };

template<class T>
class plusOp
{ public: T operator()(const T& x, const T& y) const { return x + y; } };

template<class T, class T1, class T2>
class minusOp3
{ public: T operator()(const T1& x, const T2& y) const { return x - y; } };

template<class T1, class T2>
class minusOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x - y; } };

template<class T>
class minusOp
{ public: T operator()(const T& x, const T& y) const { return x - y; } };

template<class T, class T1, class T2>
class multiplyOp3
{ public: T operator()(const T1& x, const T2& y) const { return x * y; } };

template<class T1, class T2>
class multiplyOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x * y; } };

template<class T>
class multiplyOp
{ public: T operator()(const T& x, const T& y) const { return x * y; } };

template<class T, class T1, class T2>
class divideOp3
{ public: T operator()(const T1& x, const T2& y) const { return x / y; } };

template<class T1, class T2>
class divideOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x / y; } };

template<class T>
class divideOp
{ public: T operator()(const T& x, const T& y) const { return x / y; } };

template<class T, class T1, class T2>
class cmptMultiplyOp3
{ public: T operator()(const T1& x, const T2& y) const { return cmptMultiply(x, y); } };

template<class T1, class T2>
class cmptMultiplyOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return cmptMultiply(x, y); } };

template<class T>
class cmptMultiplyOp
{ public: T operator()(const T& x, const T& y) const { return cmptMultiply(x, y); } };

template<class T, class T1, class T2>
class cmptDivideOp3
{ public: T operator()(const T1& x, const T2& y) const { return cmptDivide(x, y); } };

template<class T1, class T2>
class cmptDivideOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return cmptDivide(x, y); } };

template<class T>
class cmptDivideOp
{ public: T operator()(const T& x, const T& y) const { return cmptDivide(x, y); } };

template<class T, class T1, class T2>
class stabiliseOp3
{ public: T operator()(const T1& x, const T2& y) const { return stabilise(x, y); } };

template<class T1, class T2>
class stabiliseOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return stabilise(x, y); } };

template<class T>
class stabiliseOp
{ public: T operator()(const T& x, const T& y) const { return stabilise(x, y); } };

template<class T, class T1, class T2>
class maxOp3
{ public: T operator()(const T1& x, const T2& y) const { return max(x, y); } };

template<class T1, class T2>
class maxOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return max(x, y); } };

template<class T>
class maxOp
{ public: T operator()(const T& x, const T& y) const { return max(x, y); } };

template<class T, class T1, class T2>
class minOp3
{ public: T operator()(const T1& x, const T2& y) const { return min(x, y); } };

template<class T1, class T2>
class minOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return min(x, y); } };

template<class T>
class minOp
{ public: T operator()(const T& x, const T& y) const { return min(x, y); } };

template<class T, class T1, class T2>
class minModOp3
{ public: T operator()(const T1& x, const T2& y) const { return minMod(x, y); } };

template<class T1, class T2>
class minModOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return minMod(x, y); } };

template<class T>
class minModOp
{ public: T operator()(const T& x, const T& y) const { return minMod(x, y); } };

template<class T, class T1, class T2>
class andOp3
{ public: T operator()(const T1& x, const T2& y) const { return x && y; } };

template<class T1, class T2>
class andOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x && y; } };

template<class T>
class andOp
{ public: T operator()(const T& x, const T& y) const { return x && y; } };

template<class T, class T1, class T2>
class orOp3
{ public: T operator()(const T1& x, const T2& y) const { return x || y; } };

template<class T1, class T2>
class orOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x || y; } };

template<class T>
class orOp
{ public: T operator()(const T& x, const T& y) const { return x || y; } };

template<class T, class T1, class T2>
class eqEqOp3
{ public: T operator()(const T1& x, const T2& y) const { return x == y; } };

template<class T1, class T2>
class eqEqOp2
{ public: T1 operator()(const T1& x, const T2& y) const { return x == y; } };

template<class T>
class eqEqOp
{ public: T operator()(const T& x, const T& y) const { return x == y; } };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

InNamespace
    Foam

Description
    3D symmetric tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef symmTransform_H
#define symmTransform_H

#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar transform(const symmTensor&, const scalar s)
{
    return s;
}


template<class Cmpt>
inline Vector<Cmpt> transform(const symmTensor& stt, const Vector<Cmpt>& v)
{
    return stt & v;
}


template<class Cmpt>
inline Tensor<Cmpt> transform(const symmTensor& stt, const Tensor<Cmpt>& t)
{
    //return stt & t & stt.T();
    return Tensor<Cmpt>
    (
        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xx()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.xy()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.xz(),

        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xy()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.yy()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.yz(),

        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xz()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.yz()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.zz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xx()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.xy()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.xz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xy()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.yy()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.yz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xz()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.yz()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.zz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xx()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.xy()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.xz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xy()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.yy()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.yz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xz()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.yz()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.zz()
    );
}


template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const symmTensor& stt,
    const SphericalTensor<Cmpt>& st
)
{
    return st;
}


template<class Cmpt>
inline SymmTensor4thOrder<Cmpt> transform
(
 const symmTensor& stt,
 const SymmTensor4thOrder<Cmpt>& C
 )
{
    //- represent fourth order tensors in 6x6 notation.  Rotation is given by
    //- C_rotated_af = A_ba * C_cd * A_ef
    //- where A is a function of stt

    const scalar s = ::sqrt(2);
    const scalar A[6][6] =
    {
        { stt.xx()*stt.xx(), stt.xy()*stt.xy(), stt.xz()*stt.xz(), s*stt.xx()*stt.xy(), s*stt.xy()*stt.xz(), s*stt.xz()*stt.xx() },
        { stt.xy()*stt.xy(), stt.yy()*stt.yy(), stt.yz()*stt.yz(), s*stt.xy()*stt.yy(), s*stt.yy()*stt.yz(), s*stt.yz()*stt.xy() },
        { stt.xz()*stt.xz(), stt.yz()*stt.yz(), stt.zz()*stt.zz(), s*stt.xz()*stt.yz(), s*stt.yz()*stt.zz(), s*stt.zz()*stt.xz() },
        { s*stt.xx()*stt.xy(), s*stt.xy()*stt.yy(), s*stt.xz()*stt.yz(),
          (stt.xy()*stt.xy()+stt.xx()*stt.yy()), (stt.xz()*stt.yy()+stt.xy()*stt.yz()), (stt.xx()*stt.yz()+stt.xz()*stt.xy()) },
        { s*stt.xy()*stt.xz(), s*stt.yy()*stt.yz(), s*stt.yz()*stt.zz(),
          (stt.yy()*stt.xz()+stt.xy()*stt.yz()), (stt.yz()*stt.yz()+stt.yy()*stt.zz()), (stt.xy()*stt.zz()+stt.yz()*stt.xz()) },
        { s*stt.xz()*stt.xx(), s*stt.yz()*stt.xy(), s*stt.zz()*stt.xz(),
          (stt.yz()*stt.xx()+stt.xz()*stt.xy()), (stt.zz()*stt.xy()+stt.yz()*stt.xz()), (stt.xz()*stt.xz()+stt.zz()*stt.xx()) }
    };

    return symmTensor4thOrder
    (
        // xxxx
        A[0][0] * C.xxxx() * A[0][0] +
        A[1][0] * C.xxyy() * A[0][0] +
        A[2][0] * C.xxzz() * A[0][0] +
        A[0][0] * C.xxyy() * A[1][0] +
        A[1][0] * C.yyyy() * A[1][0] +
        A[2][0] * C.yyzz() * A[1][0] +
        A[0][0] * C.xxzz() * A[2][0] +
        A[1][0] * C.yyzz() * A[2][0] +
        A[2][0] * C.zzzz() * A[2][0] +
        A[3][0] * C.xyxy() * A[3][0] +
        A[4][0] * C.yzyz() * A[4][0] +
        A[5][0] * C.zxzx() * A[5][0],

        // xxyy
        A[0][0] * C.xxxx() * A[0][1] +
        A[1][0] * C.xxyy() * A[0][1] +
        A[2][0] * C.xxzz() * A[0][1] +
        A[0][0] * C.xxyy() * A[1][1] +
        A[1][0] * C.yyyy() * A[1][1] +
        A[2][0] * C.yyzz() * A[1][1] +
        A[0][0] * C.xxzz() * A[2][1] +
        A[1][0] * C.yyzz() * A[2][1] +
        A[2][0] * C.zzzz() * A[2][1] +
        A[3][0] * C.xyxy() * A[3][1] +
        A[4][0] * C.yzyz() * A[4][1] +
        A[5][0] * C.zxzx() * A[5][1],

        // xzz
        A[0][0] * C.xxxx() * A[0][2] +
        A[1][0] * C.xxyy() * A[0][2] +
        A[2][0] * C.xxzz() * A[0][2] +
        A[0][0] * C.xxyy() * A[1][2] +
        A[1][0] * C.yyyy() * A[1][2] +
        A[2][0] * C.yyzz() * A[1][2] +
        A[0][0] * C.xxzz() * A[2][2] +
        A[1][0] * C.yyzz() * A[2][2] +
        A[2][0] * C.zzzz() * A[2][2] +
        A[3][0] * C.xyxy() * A[3][2] +
        A[4][0] * C.yzyz() * A[4][2] +
        A[5][0] * C.zxzx() * A[5][2],

        // yyyy
        A[0][1] * C.xxxx() * A[0][1] +
        A[1][1] * C.xxyy() * A[0][1] +
        A[2][1] * C.xxzz() * A[0][1] +
        A[0][1] * C.xxyy() * A[1][1] +
        A[1][1] * C.yyyy() * A[1][1] +
        A[2][1] * C.yyzz() * A[1][1] +
        A[0][1] * C.xxzz() * A[2][1] +
        A[1][1] * C.yyzz() * A[2][1] +
        A[2][1] * C.zzzz() * A[2][1] +
        A[3][1] * C.xyxy() * A[3][1] +
        A[4][1] * C.yzyz() * A[4][1] +
        A[5][1] * C.zxzx() * A[5][1],

        // yyzz
        A[0][1] * C.xxxx() * A[0][2] +
        A[1][1] * C.xxyy() * A[0][2] +
        A[2][1] * C.xxzz() * A[0][2] +
        A[0][1] * C.xxyy() * A[1][2] +
        A[1][1] * C.yyyy() * A[1][2] +
        A[2][1] * C.yyzz() * A[1][2] +
        A[0][1] * C.xxzz() * A[2][2] +
        A[1][1] * C.yyzz() * A[2][2] +
        A[2][1] * C.zzzz() * A[2][2] +
        A[3][1] * C.xyxy() * A[3][2] +
        A[4][1] * C.yzyz() * A[4][2] +
        A[5][1] * C.zxzx() * A[5][2],

        // zzzz
        A[0][2] * C.xxxx() * A[0][2] +
        A[1][2] * C.xxyy() * A[0][2] +
        A[2][2] * C.xxzz() * A[0][2] +
        A[0][2] * C.xxyy() * A[1][2] +
        A[1][2] * C.yyyy() * A[1][2] +
        A[2][2] * C.yyzz() * A[1][2] +
        A[0][2] * C.xxzz() * A[2][2] +
        A[1][2] * C.yyzz() * A[2][2] +
        A[2][2] * C.zzzz() * A[2][2] +
        A[3][2] * C.xyxy() * A[3][2] +
        A[4][2] * C.yzyz() * A[4][2] +
        A[5][2] * C.zxzx() * A[5][2],

        // xyxy
        A[0][3] * C.xxxx() * A[0][3] +
        A[1][3] * C.xxyy() * A[0][3] +
        A[2][3] * C.xxzz() * A[0][3] +
        A[0][3] * C.xxyy() * A[1][3] +
        A[1][3] * C.yyyy() * A[1][3] +
        A[2][3] * C.yyzz() * A[1][3] +
        A[0][3] * C.xxzz() * A[2][3] +
        A[1][3] * C.yyzz() * A[2][3] +
        A[2][3] * C.zzzz() * A[2][3] +
        A[3][3] * C.xyxy() * A[3][3] +
        A[4][3] * C.yzyz() * A[4][3] +
        A[5][3] * C.zxzx() * A[5][3],

        // yzyz
        A[0][4] * C.xxxx() * A[0][4] +
        A[1][4] * C.xxyy() * A[0][4] +
        A[2][4] * C.xxzz() * A[0][4] +
        A[0][4] * C.xxyy() * A[1][4] +
        A[1][4] * C.yyyy() * A[1][4] +
        A[2][4] * C.yyzz() * A[1][4] +
        A[0][4] * C.xxzz() * A[2][4] +
        A[1][4] * C.yyzz() * A[2][4] +
        A[2][4] * C.zzzz() * A[2][4] +
        A[3][4] * C.xyxy() * A[3][4] +
        A[4][4] * C.yzyz() * A[4][4] +
        A[5][4] * C.zxzx() * A[5][4],

        // zxzx
        A[0][5] * C.xxxx() * A[0][5] +
        A[1][5] * C.xxyy() * A[0][5] +
        A[2][5] * C.xxzz() * A[0][5] +
        A[0][5] * C.xxyy() * A[1][5] +
        A[1][5] * C.yyyy() * A[1][5] +
        A[2][5] * C.yyzz() * A[1][5] +
        A[0][5] * C.xxzz() * A[2][5] +
        A[1][5] * C.yyzz() * A[2][5] +
        A[2][5] * C.zzzz() * A[2][5] +
        A[3][5] * C.xyxy() * A[3][5] +
        A[4][5] * C.yzyz() * A[4][5] +
        A[5][5] * C.zxzx() * A[5][5]
    );
}


template<class Cmpt>
inline DiagTensor<Cmpt> transform
(
    const symmTensor& stt,
    const DiagTensor<Cmpt>& dt
 )
{
    return dt;
}


template<class Cmpt>
inline SymmTensor<Cmpt> transform
(
    const symmTensor& stt,
    const SymmTensor<Cmpt>& st
)
{
    return SymmTensor<Cmpt>
    (
        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xx()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.xy()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.xz(),

        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xy()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.yy()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.yz(),

        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xz()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.yz()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.zz(),

        (stt.xy()*st.xx() + stt.yy()*st.xy() + stt.yz()*st.xz())*stt.xy()
      + (stt.xy()*st.xy() + stt.yy()*st.yy() + stt.yz()*st.yz())*stt.yy()
      + (stt.xy()*st.xz() + stt.yy()*st.yz() + stt.yz()*st.zz())*stt.yz(),

        (stt.xy()*st.xx() + stt.yy()*st.xy() + stt.yz()*st.xz())*stt.xz()
      + (stt.xy()*st.xy() + stt.yy()*st.yy() + stt.yz()*st.yz())*stt.yz()
      + (stt.xy()*st.xz() + stt.yy()*st.yz() + stt.yz()*st.zz())*stt.zz(),

        (stt.xz()*st.xx() + stt.yz()*st.xy() + stt.zz()*st.xz())*stt.xz()
      + (stt.xz()*st.xy() + stt.yz()*st.yy() + stt.zz()*st.yz())*stt.yz()
      + (stt.xz()*st.xz() + stt.yz()*st.yz() + stt.zz()*st.zz())*stt.zz()
    );
}


template<>
inline sphericalTensor transformMask<sphericalTensor>(const symmTensor& st)
{
    return sph(st);
}


template<>
inline symmTensor transformMask<symmTensor>(const symmTensor& st)
{
    return st;
}


template<>
inline symmTensor4thOrder transformMask<symmTensor4thOrder>
(
    const symmTensor& st
)
{
    notImplemented
    (
        "template<>\n"
        "inline symmTensor4thOrder transformMask<symmTensor4thOrder>"
        "(const symmTensor& st)"
    );

    return symmTensor4thOrder::zero;
}


template<>
inline diagTensor transformMask<diagTensor>(const symmTensor& st)
{
    return diagTensor(st.xx(), st.yy(), st.zz());
}


template<>
inline tensor transformMask<tensor>(const symmTensor& st)
{
    return st;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

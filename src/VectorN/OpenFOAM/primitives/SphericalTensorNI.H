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

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Cmpt, int length>
const char* const SphericalTensorN<Cmpt, length>::typeName =
    ("sphericalTensor" + name(length)).c_str();

template <class Cmpt, int length>
const SphericalTensorN<Cmpt, length> SphericalTensorN<Cmpt, length>::zero(0);

template <class Cmpt, int length>
const SphericalTensorN<Cmpt, length> SphericalTensorN<Cmpt, length>::one(1);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>::SphericalTensorN()
{}


// Construct given VectorSpace
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>::SphericalTensorN
(
    const VectorSpace<SphericalTensorN<Cmpt, length>, Cmpt, 1>& vs
)
:
    VectorSpace<SphericalTensorN<Cmpt, length>, Cmpt, 1>(vs)
{}


//- Construct from component
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>::SphericalTensorN(const Cmpt& tx)
{
    this->v_[0] = tx;
}


// Construct from Istream
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>::SphericalTensorN(Istream& is)
:
    VectorSpace<SphericalTensorN<Cmpt, length>, Cmpt, 1>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return diagonal tensor diagonal
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length> SphericalTensorN<Cmpt, length>::diag() const
{
    return *this;
}

//- Return spherical tensor transpose
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length> SphericalTensorN<Cmpt, length>::T() const
{
    return *this;
}


//- Transform the spherical tensor
//- The components are assumed to be individual scalars
//- i.e. transform has no effect
template<class Cmpt, int length>
inline SphericalTensorN<Cmpt, length> transform
(
    const tensor& tt,
    const SphericalTensorN<Cmpt, length>& v
)
{
    return v;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //


//- Addition of SphericalTensorN and SphericalTensorN
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt,length>
operator+(const SphericalTensorN<Cmpt,length>& st1, const SphericalTensorN<Cmpt,length>& st2)
{
    return SphericalTensorN<Cmpt, length>(st1.v_[0] + st2.v_[0]);
}


//- Subtraction of SphericalTensorN and SphericalTensorN
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt,length>
operator-(const SphericalTensorN<Cmpt,length>& st1, const SphericalTensorN<Cmpt,length>& st2)
{
    return SphericalTensorN<Cmpt, length>(st1.v_[0] - st2.v_[0]);
}


//- Inner-product between spherical tensor and spherical tensor
template <class Cmpt, int length>
inline typename
innerProduct<SphericalTensorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >::type
operator&(const SphericalTensorN<Cmpt, length>& st1, const SphericalTensorN<Cmpt, length>& st2)
{
    return SphericalTensorN<Cmpt, length>(st1.v_[0]*st2.v_[0]);
}


//- Inner-product between a spherical tensor and a vector
template <class Cmpt, int length>
inline typename
innerProduct<SphericalTensorN<Cmpt, length>, VectorN<Cmpt, length> >::type
operator&(const SphericalTensorN<Cmpt, length>& st, const VectorN<Cmpt, length>& v)
{
    const Cmpt& s = st.v_[0];
    VectorN<Cmpt, length> res;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::opSV(res, s, v, multiplyOp<Cmpt>());
    return res;
}


//- Inner-product between a vector and a spherical tensor
template <class Cmpt, int length>
inline typename
innerProduct<VectorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >::type
operator&(const VectorN<Cmpt, length>& v, const SphericalTensorN<Cmpt, length>& st)
{
    const Cmpt& s = st.v_[0];
    VectorN<Cmpt, length> res;
    VectorSpaceOps<VectorN<Cmpt, length>::nComponents,0>::opVS(res, v, s, multiplyOp<Cmpt>());
    return res;
}


//- Product of a scalar and a spherical tensor
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>
operator*(const scalar s, const SphericalTensorN<Cmpt, length>& st)
{
    return SphericalTensorN<Cmpt, length>(s*st.v_[0]);
}


//- Division of a scalar by a spherical tensor
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length>
operator/(const scalar s, const SphericalTensorN<Cmpt, length>& st)
{
    return SphericalTensorN<Cmpt, length>(s/st.v_[0]);
}


//- Inner Product of a VectorN by an inverse SphericalTensorN
template <class Cmpt, int length>
inline VectorN<Cmpt,length>
operator/(const VectorN<Cmpt,length>& v, const SphericalTensorN<Cmpt,length>& st)
{
    return v/st.v_[0];
}


//- Inner Product of a SphericalTensorN and an inverse SphericalTensorN
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt,length>
operator/(const SphericalTensorN<Cmpt,length>& st1, const SphericalTensorN<Cmpt,length>& st2)
{
    return SphericalTensorN<Cmpt, length>(st1.v_[0]/st2.v_[0]);
}


//- Return the inverse of a spherical tensor
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length> inv(const SphericalTensorN<Cmpt, length>& st)
{
    return SphericalTensorN<Cmpt, length>(pTraits<Cmpt>::one/st.v_[0]);
}


//- Return tensor diagonal
template <class Cmpt, int length>
inline SphericalTensorN<Cmpt, length> diag(const SphericalTensorN<Cmpt, length>& st)
{
    return st;
}

//- Return the component sum
template <class Cmpt, int length>
inline Cmpt sum(const SphericalTensorN<Cmpt, length>& st)
{
    return SphericalTensorN<Cmpt, length>::rowLength*st.v_[0];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt, int length>
class outerProduct<Cmpt, SphericalTensorN<Cmpt, length> >
{
public:

    typedef SphericalTensorN<Cmpt, length> type;
};

template<class Cmpt, int length>
class outerProduct<SphericalTensorN<Cmpt, length>, Cmpt>
{
public:

    typedef SphericalTensorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<SphericalTensorN<Cmpt, length>, VectorN<Cmpt, length> >
{
public:

    typedef VectorN<Cmpt, length> type;
};

template<class Cmpt, int length>
class innerProduct<VectorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >
{
public:

    typedef VectorN<Cmpt, length> type;
};

template<class Cmpt, int length>
class innerProduct<SphericalTensorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >
{
public:

    typedef SphericalTensorN<Cmpt, length> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

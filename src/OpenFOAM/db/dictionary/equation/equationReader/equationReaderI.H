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

inline const equationList& equationReader::eqns() const
{
    return eqns_;
}

inline const UPtrList<const dictionary>& equationReader::dictSources() const
{
    return dictSources_;
}


inline UPtrList<const dictionary>& equationReader::dictSources()
{
    return dictSources_;
}


inline const PtrList<word>& equationReader::dictLookups() const
{
    return dictLookups_;
}


inline const
    UPtrList<const scalarList>& equationReader::externalScalarLists() const
{
    return externalScalarLists_;
}


inline UPtrList<const scalarList>& equationReader::externalScalarLists()
{
    return externalScalarLists_;
}


inline const
    PtrList<dimensionSet>& equationReader::externalScalarListDimensions() const
{
    return externalScalarListDimensions_;
}


inline PtrList<dimensionSet>& equationReader::externalScalarListDimensions()
{
    return externalScalarListDimensions_;
}


inline const wordList& equationReader::externalScalarListNames() const
{
    return externalScalarListNames_;
}


inline wordList& equationReader::externalScalarListNames()
{
    return externalScalarListNames_;
}


inline const labelList& equationReader::externalScalarListIndex() const
{
    return externalScalarListIndex_;
}


inline labelList& equationReader::externalScalarListIndex()
{
    return externalScalarListIndex_;
}


inline const
    UPtrList<const dimensionedScalar>& equationReader::externalDScalars() const
{
    return externalDScalars_;
}


inline UPtrList<const dimensionedScalar>& equationReader::externalDScalars()
{
    return externalDScalars_;
}


inline const UPtrList<const scalar>& equationReader::externalScalars() const
{
    return externalScalars_;
}


inline UPtrList<const scalar>& equationReader::externalScalars()
{
    return externalScalars_;
}


inline const wordList& equationReader::externalScalarNames() const
{
    return externalScalarNames_;
}


inline wordList& equationReader::externalScalarNames()
{
    return externalScalarNames_;
}



inline const
    PtrList<dimensionSet>& equationReader::externalScalarDimensions() const
{
    return externalScalarDimensions_;
}


inline PtrList<dimensionSet>& equationReader::externalScalarDimensions()
{
    return externalScalarDimensions_;
}


inline const
    UPtrList<scalar>& equationReader::outputScalars() const
{
    return outputScalars_;
}


inline UPtrList<scalar>& equationReader::outputScalars()
{
    return outputScalars_;
}


inline const
    wordList& equationReader::outputScalarNames() const
{
    return outputScalarNames_;
}


inline
    wordList& equationReader::outputScalarNames()
{
    return outputScalarNames_;
}


inline const
    PtrList<dimensionSet>& equationReader::outputScalarDimensions() const
{
    return outputScalarDimensions_;
}


inline
    PtrList<dimensionSet>& equationReader::outputScalarDimensions()
{
    return outputScalarDimensions_;
}


inline const
    UPtrList< scalarList>& equationReader::outputScalarLists() const
{
    return outputScalarLists_;
}

//- Access to pointers to output scalarLists
inline UPtrList<scalarList>& equationReader::outputScalarLists()
{
    return outputScalarLists_;
}

//- Const access to dimensions associated with outputScalarLists
inline const
    PtrList<dimensionSet>& equationReader::outputScalarListDimensions() const
{
    return outputScalarListDimensions_;
}

//- Access to dimensions associated with outputScalarLists
inline PtrList<dimensionSet>& equationReader::outputScalarListDimensions()
{
    return outputScalarListDimensions_;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

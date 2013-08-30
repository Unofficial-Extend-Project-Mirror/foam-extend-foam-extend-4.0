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

inline const word& equation::name() const
{
    return equationName_;
}


inline label equation::size() const
{
    return ops_.size();
}


inline void equation::setSize(const label newSize) const
{
    ops_.setSize(newSize);
}


inline autoPtr<equationOperation> equation::set
(
    const label elementIndex,
    equationOperation * newOperation
) const
{
    return ops_.set(elementIndex, newOperation);
}


inline word& equation::name()
{
    return equationName_;
}


inline const string& equation::rawText() const
{
    return rawText_;
}


inline void equation::setRawText(const string& newRawText) const
{
    rawText_ = newRawText;
}


inline const dimensionedScalar& equation::lastResult() const
{
    return lastResult_;
}


inline void equation::setLastResult(const word& newName) const
{
    lastResult_.name() = newName;
}


inline void equation::setLastResult(const dimensionSet& newDims) const
{
    lastResult_.dimensions().reset(newDims);
}


inline void equation::setLastResult(const scalar& newScalar) const
{
    lastResult_.value() = newScalar;
}


inline const dimensionSet& equation::overrideDimensions() const
{
    return overrideDimensions_;
}

inline void equation::setOverrideDimensions
(
    const dimensionSet& newDims
) const
{
    overrideDimensions_.reset(newDims);
}


inline const bool& equation::changeDimensions() const
{
    return changeDimensions_;
}


inline void equation::setChangeDimensions(bool newFlag) const
{
    changeDimensions_ = newFlag;
}


inline const labelList& equation::maxFieldSizes() const
{
    return maxFieldSizes_;
}


inline void equation::setMaxFieldSizes
(
    const labelList& newSizes
) const
{
    maxFieldSizes_ = newSizes;
}

// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

inline equationOperation& equation::operator[](const label equationIndex) const
{
    return ops_[equationIndex];
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

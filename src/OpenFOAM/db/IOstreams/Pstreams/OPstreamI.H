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


\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

template<class T>
inline void OPstream::writeToBuffer(const T& t)
{
    //(T&)(buf_[bufPosition_]) = t;
    //bufPosition_ += sizeof(T);

    writeToBuffer(&t, sizeof(T), sizeof(T));
}


inline void OPstream::writeToBuffer(const char& c)
{
    if (size_t(buf_.size()) < bufPosition_ + 1U)
    {
        enlargeBuffer(1);
    }

    buf_[bufPosition_] = c;
    bufPosition_ ++;
}


inline void OPstream::writeToBuffer
(
    const void* data,
    size_t count,
    size_t align
)
{
    label oldPos = bufPosition_;

    if (align > 1)
    {
        // Align bufPosition. Pads bufPosition_-oldPos characters.
        bufPosition_ = align + ((bufPosition_ - 1) & ~(align - 1));
    }

    if (size_t(buf_.size()) < bufPosition_ + count)
    {
        enlargeBuffer(bufPosition_ - oldPos + count);
    }

    register char* bufPtr = &buf_[bufPosition_];
    register const char* dataPtr = reinterpret_cast<const char*>(data);
    register size_t i = count;
    while (i--) *bufPtr++ = *dataPtr++;

    bufPosition_ += count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidParticle::trackData::trackData
(
    solidParticleCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& nuInterp,
    const vector& g
)
:
    spc_(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    nuInterp_(nuInterp),
    g_(g)
{}


inline Foam::solidParticle::solidParticle
(
    const Cloud<solidParticle>& c,
    const vector& position,
    const label celli,
    const scalar d,
    const vector& U
)
:
    Particle<solidParticle>(c, position, celli),
    d_(d),
    U_(U)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::solidParticleCloud& Foam::solidParticle::trackData::spc()
{
    return spc_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::solidParticle::trackData::rhoInterp() const
{
    return rhoInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::solidParticle::trackData::UInterp() const
{
    return UInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::solidParticle::trackData::nuInterp() const
{
    return nuInterp_;
}

inline const Foam::vector& Foam::solidParticle::trackData::g() const
{
    return g_;
}


inline Foam::scalar Foam::solidParticle::d() const
{
    return d_;
}

inline Foam::scalar Foam::solidParticle::wallImpactDistance(const vector&) const
{
    return 0.5*d_;
}

inline const Foam::vector& Foam::solidParticle::U() const
{
    return U_;
}


// ************************************************************************* //

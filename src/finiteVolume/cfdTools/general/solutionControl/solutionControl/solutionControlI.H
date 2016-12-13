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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary& Foam::solutionControl::dict() const
{
    return mesh_.solutionDict().subDict(algorithmName_);
}


inline Foam::label Foam::solutionControl::corr() const
{
    return corr_;
}


inline Foam::label Foam::solutionControl::corrNonOrtho() const
{
    return corrNonOrtho_;
}


inline Foam::label Foam::solutionControl::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline bool Foam::solutionControl::finalNonOrthogonalIter() const
{
    return corrNonOrtho_ == nNonOrthCorr_ + 1;
}


inline bool Foam::solutionControl::momentumPredictor() const
{
    return momentumPredictor_;
}


inline bool Foam::solutionControl::transonic() const
{
    return transonic_;
}


inline bool Foam::solutionControl::consistent() const
{
    return consistent_;
}


inline bool Foam::solutionControl::correctNonOrthogonal()
{
    corrNonOrtho_++;

    if (debug)
    {
        Info<< algorithmName_ << " correctNonOrthogonal: corrNonOrtho = "
            << corrNonOrtho_ << endl;
    }

    if (corrNonOrtho_ <= nNonOrthCorr_ + 1)
    {
        return true;
    }
    else
    {
        corrNonOrtho_ = 0;
        return false;
    }
}


// ************************************************************************* //

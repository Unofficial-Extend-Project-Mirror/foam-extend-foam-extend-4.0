# -----------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.0
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     doxyFilt-top.awk
#
# Description
#     Only output the first /* ... */ comment section found in the file
#     Use @cond / @endcond to suppress documenting all classes/variables
#     - This is useful for application files in which only the first
#       block documents the application itself.
#
# -----------------------------------------------------------------------------
BEGIN {
    state = 0
}

# a '/*' at the beginning of a line starts a comment block
/^ *\/\*/ {
   state++
}

# check first line
# either started with a comment or skip documentation for the whole file
FNR == 1 {
   if (!state)
   {
      print "//! @cond OpenFOAMIgnoreAppDoxygen"
      state = 2
   }
}

# a '*/' ends the comment block
# skip documentation for rest of the file
/\*\// {
    if (state == 1)
    {
        print
        print "//! @cond OpenFOAMIgnoreAppDoxygen"
    }
    state = 2
    next
}

# print everything within the first comment block
{
    if (state)
    {
        print
    }
    next
}

END {
    if (state == 2)
    {
        print "//! @endcond OpenFOAMIgnoreAppDoxygen"
    }
}

# -----------------------------------------------------------------------------

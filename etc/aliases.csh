#----------------------------------*-sh-*--------------------------------------
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
#     aliases.csh
#
# Description
#     Aliases for working with FOAM
#     Sourced from FOAM-??/etc/cshrc and/or ~/.cshrc
#
#------------------------------------------------------------------------------

# Change compiled version aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wm64 'setenv WM_ARCH_OPTION 64; source $WM_PROJECT_DIR/etc/cshrc'
alias wm32 'setenv WM_ARCH_OPTION 32; source $WM_PROJECT_DIR/etc/cshrc'
alias wmSP 'setenv WM_PRECISION_OPTION SP; source $WM_PROJECT_DIR/etc/cshrc'
alias wmDP 'setenv WM_PRECISION_OPTION DP; source $WM_PROJECT_DIR/etc/cshrc'

# Toggle wmakeScheduler on/off
#  - also need to set WM_HOSTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wmSchedON 'setenv WM_SCHEDULER $WM_PROJECT_DIR/wmake/wmakeScheduler'
alias wmSchedOFF 'unsetenv WM_SCHEDULER'

# Change directory aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~
alias src 'cd $FOAM_SRC'
alias lib 'cd $FOAM_LIB'
alias run 'cd $FOAM_RUN'
alias foam 'cd $WM_PROJECT_DIR'
alias foamsrc 'cd $FOAM_SRC/$WM_PROJECT'
alias foamfv 'cd $FOAM_SRC/finiteVolume'
alias app 'cd $FOAM_APP'
alias util 'cd $FOAM_UTILITIES'
alias sol 'cd $FOAM_SOLVERS'
alias tut 'cd $FOAM_TUTORIALS'
alias foam3rdParty 'cd $WM_THIRD_PARTY_DIR'
alias work 'cd $WM_PROJECT_USER_DIR/work'

# -----------------------------------------------------------------------------

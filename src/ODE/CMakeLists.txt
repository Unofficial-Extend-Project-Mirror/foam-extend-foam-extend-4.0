# --------------------------------------------------------------------------
#   ========                 |
#   \      /  F ield         | foam-extend: Open Source CFD
#    \    /   O peration     | Version:     4.1
#     \  /    A nd           | Web:         http://www.foam-extend.org
#      \/     M anipulation  | For copyright notice see file Copyright
# --------------------------------------------------------------------------
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
# Description
#     CMakeLists.txt file for libraries and applications
#
# Author
#     Henrik Rusche, Wikki GmbH, 2017. All rights reserved
#
#
# --------------------------------------------------------------------------

set(ODESolvers ODESolvers)
list(APPEND SOURCES
  ${ODESolvers}/ODESolver/ODESolver.C
  ${ODESolvers}/ODESolver/newODESolver.C
  ${ODESolvers}/Euler/Euler.C
  ${ODESolvers}/RK/RK.C
  ${ODESolvers}/KRR4/KRR4.C
  ${ODESolvers}/SIBS/SIBS.C
  ${ODESolvers}/SIBS/SIMPR.C
  ${ODESolvers}/SIBS/polyExtrapolate.C
)

set(translationODE translationODE)
list(APPEND SOURCES
  ${translationODE}/translationODE.C
)

set(sixDOF sixDOF)
list(APPEND SOURCES
  ${sixDOF}/finiteRotation/finiteRotation.C
  ${sixDOF}/sixDOFqODE/sixDOFqODE.C
  ${sixDOF}/sixDOFODE/constraints/rotationalConstraints/rotationalConstraint/rotationalConstraint.C
  ${sixDOF}/sixDOFODE/constraints/rotationalConstraints/constantAngularAcceleration/constantAngularAcceleration.C
  ${sixDOF}/sixDOFODE/constraints/translationalConstraints/translationalConstraint/translationalConstraint.C
  ${sixDOF}/sixDOFODE/constraints/translationalConstraints/constantTranslationalAcceleration/constantTranslationalAcceleration.C
  ${sixDOF}/sixDOFODE/constraints/translationalConstraints/periodicOscillation/periodicOscillation.C

  ${sixDOF}/sixDOFODE/restraints/translationalRestraints/translationalRestraint/translationalRestraint.C
  ${sixDOF}/sixDOFODE/restraints/translationalRestraints/linearSpringDamper/linearSpringDamper.C
  ${sixDOF}/sixDOFODE/restraints/rotationalRestraints/rotationalRestraint/rotationalRestraint.C
  ${sixDOF}/sixDOFODE/restraints/rotationalRestraints/angularDamper/angularDamper.C

  ${sixDOF}/sixDOFODE/sixDOFODE.C
  ${sixDOF}/sixDOFODE/newSixDOFODE.C
  ${sixDOF}/quaternionSixDOF/quaternionSixDOF.C
  ${sixDOF}/geometricSixDOF/geometricSixDOF.C

  ${sixDOF}/sixDOFBodies/sixDOFBodies.C
)

add_foam_library(ODE SHARED ${SOURCES})

target_link_libraries(ODE PUBLIC meshTools)

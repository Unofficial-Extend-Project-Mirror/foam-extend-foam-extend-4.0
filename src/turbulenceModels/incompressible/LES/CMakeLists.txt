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

list(APPEND SOURCES
  vanDriestDelta/vanDriestDelta.C
  LESModel/LESModel.C
  GenEddyVisc/GenEddyVisc.C
  GenSGSStress/GenSGSStress.C
  laminar/laminar.C
  SpalartAllmaras/SpalartAllmaras.C
  SpalartAllmarasDDES/SpalartAllmarasDDES.C
  SpalartAllmarasIDDES/SpalartAllmarasIDDES.C
  SpalartAllmarasIDDES/IDDESDelta/IDDESDelta.C
  oneEqEddy/oneEqEddy.C
  dynOneEqEddy/dynOneEqEddy.C
  locDynOneEqEddy/locDynOneEqEddy.C
  Smagorinsky/Smagorinsky.C
  dynSmagorinsky/dynSmagorinsky.C
  LRRDiffStress/LRRDiffStress.C
  DeardorffDiffStress/DeardorffDiffStress.C
  spectEddyVisc/spectEddyVisc.C
  scaleSimilarity/scaleSimilarity.C
  mixedSmagorinsky/mixedSmagorinsky.C
  dynMixedSmagorinsky/dynMixedSmagorinsky.C
)

# Smagorinsky2/Smagorinsky2.C

list(APPEND SOURCES
  kOmegaSSTSAS/kOmegaSSTSAS.C
)

#  Wall functions
set(wallFunctions derivedFvPatchFields/wallFunctions)

list(APPEND SOURCES
  derivedFvPatchFields/decayingTurbulence/decayingVorton.C
  derivedFvPatchFields/decayingTurbulence/decayingTurbulenceFvPatchVectorField.C
)

set(nuSgsWallFunctions ${wallFunctions}/nuSgsWallFunctions)
list(APPEND SOURCES
  ${nuSgsWallFunctions}/nuSgsWallFunction/nuSgsWallFunctionFvPatchScalarField.C
)

add_foam_library(incompressibleLESModels SHARED ${SOURCES})

target_link_libraries(incompressibleLESModels PUBLIC incompressibleTurbulenceModel LESdeltas LESfilters)

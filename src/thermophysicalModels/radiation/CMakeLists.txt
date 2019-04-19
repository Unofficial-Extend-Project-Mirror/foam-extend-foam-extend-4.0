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

#  Radiation constants
list(APPEND SOURCES
  radiationConstants/radiationConstants.C
)

#  Radiation model
list(APPEND SOURCES
  radiationModel/radiationModel/radiationModel.C
  radiationModel/radiationModel/newRadiationModel.C
  radiationModel/noRadiation/noRadiation.C
  radiationModel/P1/P1.C
  radiationModel/fvDOM/fvDOM/fvDOM.C
  radiationModel/fvDOM/radiativeIntensityRay/radiativeIntensityRay.C
  radiationModel/fvDOM/blackBodyEmission/blackBodyEmission.C
  radiationModel/fvDOM/absorptionCoeffs/absorptionCoeffs.C
  radiationModel/viewFactor/viewFactor.C
)

#  Scatter model
list(APPEND SOURCES
  submodels/scatterModel/scatterModel/scatterModel.C
  submodels/scatterModel/scatterModel/newScatterModel.C
  submodels/scatterModel/constantScatter/constantScatter.C
)


#  Absorption/Emission model
list(APPEND SOURCES
  submodels/absorptionEmissionModel/absorptionEmissionModel/absorptionEmissionModel.C
  submodels/absorptionEmissionModel/absorptionEmissionModel/newAbsorptionEmissionModel.C
  submodels/absorptionEmissionModel/noAbsorptionEmission/noAbsorptionEmission.C
  submodels/absorptionEmissionModel/constantAbsorptionEmission/constantAbsorptionEmission.C
  submodels/absorptionEmissionModel/binaryAbsorptionEmission/binaryAbsorptionEmission.C
  submodels/absorptionEmissionModel/greyMeanAbsorptionEmission/greyMeanAbsorptionEmission.C
  submodels/absorptionEmissionModel/wideBandAbsorptionEmission/wideBandAbsorptionEmission.C
)


#  Boundary conditions
list(APPEND SOURCES
  derivedFvPatchFields/MarshakRadiation/MarshakRadiationMixedFvPatchScalarField.C
  derivedFvPatchFields/MarshakRadiationFixedT/MarshakRadiationFixedTMixedFvPatchScalarField.C
  derivedFvPatchFields/greyDiffusiveRadiation/greyDiffusiveRadiationMixedFvPatchScalarField.C
  derivedFvPatchFields/wideBandDiffusiveRadiation/wideBandDiffusiveRadiationMixedFvPatchScalarField.C
  derivedFvPatchFields/wideBandSpecularRadiation/wideBandSpecularRadiationMixedFvPatchScalarField.C
  derivedFvPatchFields/greyDiffusiveViewFactor/greyDiffusiveViewFactorFixedValueFvPatchScalarField.C
)

add_foam_library(radiation SHARED ${SOURCES})

target_link_libraries(radiation PUBLIC basicThermophysicalModels)

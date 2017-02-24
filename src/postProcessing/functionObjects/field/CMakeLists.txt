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
#     Free Software Foundation, either version 3 of the License, or {at your
#     option} any later version.
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
  fieldAverage/fieldAverage/fieldAverage.C
  fieldAverage/fieldAverageItem/fieldAverageItem.C
  fieldAverage/fieldAverageItem/fieldAverageItemIO.C
  fieldAverage/fieldAverageFunctionObject/fieldAverageFunctionObject.C
  fieldMinMax/fieldMinMax.C
  fieldMinMax/fieldMinMaxFunctionObject.C
  minMaxField/minMaxField.C
  maxFieldCell/maxFieldCell.C
  fieldValues/fieldValue/fieldValue.C
  fieldValues/faceSource/faceSource.C
  fieldValues/faceSource/faceSourceFunctionObject.C
  fieldValues/cellSource/cellSource.C
  fieldValues/cellSource/cellSourceFunctionObject.C
  readFields/readFields.C
  readFields/readFieldsFunctionObject.C
  surfaceInterpolateFields/surfaceInterpolateFieldsFunctionObject.C
  surfaceInterpolateFields/surfaceInterpolateFields.C
)

add_foam_library(fieldFunctionObjects SHARED ${SOURCES})

target_link_libraries(fieldFunctionObjects PUBLIC finiteVolume)

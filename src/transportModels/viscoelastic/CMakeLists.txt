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
  viscoelasticModel/viscoelasticModel.C
  viscoelasticLaws/viscoelasticLaw/viscoelasticLaw.C
  viscoelasticLaws/viscoelasticLaw/newViscoelasticLaw.C
  viscoelasticLaws/LPTT/LPTT.C
  viscoelasticLaws/EPTT/EPTT.C
  viscoelasticLaws/Oldroyd-B/Oldroyd_B.C
  viscoelasticLaws/linearMaxwell/linearMaxwell.C
  viscoelasticLaws/UCM/UCM.C
  viscoelasticLaws/Giesekus/Giesekus.C
  viscoelasticLaws/FENE-CR/FENE_CR.C
  viscoelasticLaws/FENE-P/FENE_P.C
  viscoelasticLaws/XPP_SE/XPP_SE.C
  viscoelasticLaws/XPP_DE/XPP_DE.C
  viscoelasticLaws/DCPP/DCPP.C
  viscoelasticLaws/Feta-PTT/Feta_PTT.C
  viscoelasticLaws/Leonov/Leonov.C
  viscoelasticLaws/WhiteMetzner/WhiteMetznerLarson/WhiteMetznerLarson.C
  viscoelasticLaws/WhiteMetzner/WhiteMetznerCross/WhiteMetznerCross.C
  viscoelasticLaws/WhiteMetzner/WhiteMetznerCarreauYasuda/WhiteMetznerCarreauYasuda.C
  viscoelasticLaws/S_MDCPP/S_MDCPP.C
  viscoelasticLaws/multiMode/multiMode.C
)

add_foam_library(viscoelasticTransportModels SHARED ${SOURCES})

target_link_libraries(viscoelasticTransportModels PUBLIC finiteVolume)

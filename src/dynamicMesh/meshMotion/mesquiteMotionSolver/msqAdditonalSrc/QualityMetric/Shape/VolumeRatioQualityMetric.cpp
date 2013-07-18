/* *****************************************************************
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2011 Oliver Borm

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    oli.borm@web.de

  ***************************************************************** */
/*!
  \file   VolumeRatioQualityMetric.cpp

  \brief Evaluates the Volume Expansion Ratio metric for two- and
  three-diminsional elements.
  \author Oliver Borm
  \date   2011-01-15
*/

#include "VolumeRatioQualityMetric.hpp"
#include <math.h>
#include "MsqMeshEntity.hpp"
#include "PatchData.hpp"

#include <vector>
using std::vector;

using namespace Mesquite;

std::string VolumeRatioQualityMetric::get_name() const
  { return "VolumeRatio"; }

int VolumeRatioQualityMetric::get_negate_flag() const
  { return -1; }

//!Evaluate expansion volume ratio on ``element''
bool VolumeRatioQualityMetric::evaluate(PatchData &pd,
                                             size_t elem_index,
                                             double &fval,
                                             MsqError &err)
{
  MsqMeshEntity& ownElement = pd.element_by_index( elem_index );

  double volOwner = 1;
  double volNei = 1;
  double volRatio = 1;

  fval = 1;

//   print_patch_data(pd);

//   std::cout << "elem_index: " << elem_index << std::endl;

  std::vector<size_t> neighbourElements;
  pd.get_element_to_element_indices(elem_index,neighbourElements,err);
  MSQ_ERRZERO(err);

//   std::cout << "neighbourElements: ";
//   for (size_t cellI=0; cellI < neighbourElements.size(); cellI++)
//   {
//     std::cout << neighbourElements[cellI] << " ";
//   }
//   std::cout << std::endl;

//   std::vector<size_t> neighbourElementsNew;
//   EntityTopology entity = ownElement.get_element_type();
//   size_t connectDimension = TopologyInfo::dimension(entity)-1;
//   // 1 = connected by lines for 2D elements; 2 = connected by faces for 3D elements
//   pd.get_adjacent_entities_via_n_dim(connectDimension,elem_index,neighbourElementsNew,err);
//   MSQ_ERRZERO(err);
//
//   std::cout << "neighbourElementsNew: ";
//   for (size_t cellI=0; cellI < neighbourElementsNew.size(); cellI++)
//   {
//     std::cout << neighbourElementsNew[cellI] << " ";
//   }
//   std::cout << std::endl;

  // compute positive volume of owner cell
//   volOwner = ownElement.compute_unsigned_area(pd, err);
  volOwner = ownElement.compute_signed_area(pd, err);
//   MSQ_ERRRTN(err);

  for (size_t cellI=0; cellI < neighbourElements.size(); cellI++)
  {
    MsqMeshEntity& neiElement = pd.element_by_index(neighbourElements[cellI]);
//     volNei = neiElement.compute_unsigned_area(pd, err);
    volNei = neiElement.compute_signed_area(pd, err);
//   MSQ_ERRRTN(err);
    volRatio = volOwner/(volNei+1e-15); //+SMALL
    if (fabs(volRatio) < 1.0) {volRatio = 1.0/(volRatio+1e-15);}
    fval = std::max(fval,volRatio);
  }

  return true;
}

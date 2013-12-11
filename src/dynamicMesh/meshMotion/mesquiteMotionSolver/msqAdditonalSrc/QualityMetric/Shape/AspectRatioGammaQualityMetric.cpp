/* *****************************************************************
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000
    with Sandia Corporation, the U.S. Government retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#include "AspectRatioGammaQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "MsqMeshEntity.hpp"
#include "PatchData.hpp"

#include <vector>
using std::vector;

using namespace Mesquite;

const double fourDivRootThree = 4.0/sqrt(3.0);
const double twelveDivRootTwo = 12.0/sqrt(2.0);
const double sixDivRootTwo = 6.0/sqrt(2.0);


std::string AspectRatioGammaQualityMetric::get_name() const
  { return "AspectRatioGamma"; }

int AspectRatioGammaQualityMetric::get_negate_flag() const
  { return 1; }

//note that we can define this metric for other element types?
//!Evaluate aspect ratio gamma on ``element''
bool AspectRatioGammaQualityMetric::evaluate(PatchData &pd,
                                             size_t elem_index,
                                             double &fval,
                                             MsqError &err)
{
  MsqMeshEntity& element = pd.element_by_index( elem_index );
  EntityTopology entity = element.get_element_type();
  double vol;
  Vector3D cross, normal(0,0,0);
  fval = MSQ_MAX_CAP;

  //get element's nodes
  vert.clear();
  pd.get_element_vertex_coordinates(elem_index, vert, err);  MSQ_ERRZERO(err);

  switch(entity)
  {
    case TRIANGLE:
      //area
      cross = (vert[1] - vert[0]) * (vert[2] - vert[0]);
      vol= cross.length() / 2.0;
      if (vol < MSQ_MIN)
        return false;

      if (pd.domain_set()) { // need domain to check for inverted elements
        pd.get_domain_normal_at_corner( elem_index, 0, normal, err ); MSQ_ERRZERO(err);
        if ((cross % normal) < -MSQ_MIN)
          return false;
      }

      // sum squares of edge lengths
      fval = ((vert[1] - vert[0]).length_squared()
            + (vert[2] - vert[0]).length_squared()
            + (vert[1] - vert[2]).length_squared());
      // average
      fval /= 3.0;

      //normalize to equil. and div by area
      fval /= vol * fourDivRootThree;

      break;

    case TETRAHEDRON:
      vol = (vert[1] - vert[0]) % ((vert[2] - vert[0]) * (vert[3] - vert[0])) / 6.0;
      if (vol < MSQ_MIN)  // zero for degenerate and negative for inverted
        return false;

      // sum squares of edge lengths
      fval = (vert[1] - vert[0]).length_squared()
           + (vert[2] - vert[0]).length_squared()
           + (vert[3] - vert[0]).length_squared()
           + (vert[2] - vert[1]).length_squared()
           + (vert[3] - vert[1]).length_squared()
           + (vert[3] - vert[2]).length_squared();
        // average
      fval /= 6.0;

      fval *= sqrt(fval);
      //normalize to equil. and div by volume
      fval /= vol * twelveDivRootTwo;

      break;

    case PRISM:
      // ideal shape is a prism with equilateral triangle base and equal length height
      // this volume is always positive with this method
      vol = element.compute_signed_area(pd, err);
      if (vol< MSQ_MIN)  // zero for degenerate and negative for inverted
        return false;

      fval = (vert[1] - vert[0]).length_squared()
           + (vert[2] - vert[0]).length_squared()
           + (vert[1] - vert[2]).length_squared()
           + (vert[4] - vert[3]).length_squared()
           + (vert[5] - vert[3]).length_squared()
           + (vert[4] - vert[5]).length_squared()
           + (vert[5] - vert[2]).length_squared()
           + (vert[3] - vert[0]).length_squared()
           + (vert[4] - vert[1]).length_squared();

      // average L2 edge length
      fval /= 9.0;

      fval *= sqrt(fval);
      //normalize to equil. and div by volume
      fval /= vol * fourDivRootThree;

      break;

    case HEXAHEDRON:
      // cube is an ideal shape
      // this volume is always positive with this method
      vol = element.compute_unsigned_area(pd, err);;
      if (vol< MSQ_MIN)  // zero for degenerate and negative for inverted
        return false;

      fval = (vert[1] - vert[0]).length_squared()
           + (vert[2] - vert[1]).length_squared()
           + (vert[3] - vert[2]).length_squared()
           + (vert[3] - vert[0]).length_squared()
           + (vert[5] - vert[4]).length_squared()
           + (vert[6] - vert[5]).length_squared()
           + (vert[7] - vert[6]).length_squared()
           + (vert[7] - vert[4]).length_squared()
           + (vert[6] - vert[2]).length_squared()
           + (vert[5] - vert[1]).length_squared()
           + (vert[4] - vert[0]).length_squared()
           + (vert[7] - vert[3]).length_squared();

      // average L2 edge length
      fval /= 12.0;

      fval *= sqrt(fval);
      //normalize to equil. and div by volume
      fval /= vol;

      break;

    case PYRAMID:
      // Johnson Pyramid with equilateral triangle sides and square base
      // this volume is always positive with this method
      vol = element.compute_unsigned_area(pd, err);
      if (vol< MSQ_MIN)  // zero for degenerate and negative for inverted
        return false;

      fval = (vert[1] - vert[0]).length_squared()
           + (vert[2] - vert[1]).length_squared()
           + (vert[3] - vert[2]).length_squared()
           + (vert[0] - vert[3]).length_squared()
           + (vert[4] - vert[0]).length_squared()
           + (vert[4] - vert[1]).length_squared()
           + (vert[4] - vert[2]).length_squared()
           + (vert[4] - vert[3]).length_squared();

      // average L2 edge length
      fval /= 8.0;

      fval *= sqrt(fval);
      //normalize to equil. and div by volume
      fval /= vol * sixDivRootTwo;

     break ;

    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                     "Entity type %d is not valid for Aspect Ratio Gamma\n",
                     (int)entity);
      return false;
  };

  return true;
}

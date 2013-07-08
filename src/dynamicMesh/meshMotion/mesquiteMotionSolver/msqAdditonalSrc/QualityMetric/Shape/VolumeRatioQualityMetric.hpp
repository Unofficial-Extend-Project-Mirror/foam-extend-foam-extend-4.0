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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file VolumeRatioQualityMetric.hpp
  \brief
  Header file for the Mesquite::VolumeRatioQualityMetric class

  \author Oliver Borm
  \date   2011-01-17
 */


#ifndef VolumeRatioQualityMetric_hpp
#define VolumeRatioQualityMetric_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "ElementQM.hpp"
namespace MESQUITE_NS
{
     /*! \class VolumeRatioQualityMetric
       \brief Object for computing the maximum volume ratio of
       an element, by checking all neighbour elements.
     */
   class VolumeRatioQualityMetric : public ElementQM
   {
   public:     
     VolumeRatioQualityMetric() {}

       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~VolumeRatioQualityMetric()
        {}

     virtual std::string get_name() const;

     int get_negate_flag() const;

     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );

   private:
      std::vector<Vector3D> vert;
   };


} //namespace


#endif // VolumeRatioQualityMetric_hpp



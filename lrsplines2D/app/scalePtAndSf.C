/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

//#define DEBUG

using namespace Go;
using std::vector;
using std::string;

int main(int argc, char *argv[])
{
  if (argc != 6) 
    {
      std::cout << "Parameters: insurf(.g2), inpoints(.g2) outsurf(.g2) outpoints(.g2) scaling_factor" << std::endl;
      return 1;
    }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream sfout(argv[3]);
  std::ofstream ptsout(argv[4]);
  double scale = atof(argv[5]);

   // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(sfin);
  }
  catch (...)
    {
      std::cout << "ERROR: Input object not recognized. Exiting" << std::endl;
      return 1;
    }

  shared_ptr<GeomObject> geom_obj;
   try {
    geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
    geom_obj->read(sfin);
  }
  catch (...)
    {
      std::cout << "ERROR: Input surface could not be read. Exiting" << std::endl;
      return 1;
    }
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cout << "ERROR: Input file contains no surface" << std::endl;
      return 1;
    }

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  shared_ptr<BoundedSurface> bdsf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  RectDomain dom;
  if (bdsf.get())
    dom = bdsf->underlyingSurface()->containingDomain();
  else
    dom = sf->containingDomain();
  double umin = dom.umin()*scale;
  double umax = dom.umax()*scale;
  double vmin = dom.vmin()*scale;
  double vmax = dom.vmax()*scale;
  sf->setParameterDomain(umin, umax, vmin, vmax);

  double *data = points.rawData();
  int nmb = points.numPoints();

  int ki;
  for (ki=0; ki<nmb; ++ki)
    {
      data[3*ki] *= scale;
      data[3*ki+1] *= scale;
    }

  sf->writeStandardHeader(sfout);
  sf->write(sfout);

  points.writeStandardHeader(ptsout);
  points.write(ptsout);
}

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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace Go;
using std::vector;



int main(int argc, char *argv[])
{
  if (argc != 5)
  {
      std::cout << "Usage: input_(lr_)spline_sf.g2, upar, vpar number of derivatives " << std::endl;
      return -1;
  }

  std::ifstream filein(argv[1]); // Input LRSplineSurface.
  double upar = atof(argv[2]);
  double vpar = atof(argv[3]);
  int nder = atoi(argv[4]);

  Go::ObjectHeader header;
  filein >> header;

  LRSplineSurface lrspline_sf;
  lrspline_sf.read(filein);

  std::cout << "Done reading lrspline_sf." << std::endl;

  Element2D *elem = lrspline_sf.coveringElement(upar, vpar);
  const vector<LRBSpline2D*>& Bfunctions = elem->getSupport();

  Point zero(lrspline_sf.dimension());
  zero.setValue(0.0);
	       
  for (size_t ki=0; ki<Bfunctions.size(); ++ki)
    {
      vector<double> der1((nder+1)*(nder+2));
      size_t kj=0;
      for (int ka=0; ka<=nder; ++ka)
	for (int kb=0; kb<=nder; ++kb, ++kj)
	  der1[kj] = Bfunctions[ki]->evalBasisFunction(upar, vpar, ka, kb,
						       false, false);

      vector<Point> der2((nder+1)*(nder+2)/2, zero);
      Bfunctions[ki]->evalder_add(upar, vpar, nder, &der2[0], false, false);
      int stop_break = 1;
    }
}


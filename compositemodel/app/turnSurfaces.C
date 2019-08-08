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

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/errormacros.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file on (g22/g2) format, output file, file type" << std::endl;
    exit(-1);
  }

  // Read input arguments
  char* file1(argv[1]);

  std::ofstream file2(argv[2]);
  ALWAYS_ERROR_IF(file2.bad(), "Bad or no output filename");

  int type = atoi(argv[3]);

  if (type == 1)
    {
      std::ifstream is(file1);
      IGESconverter conv;
      conv.readgo(is);
      vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();
      int nmbgeom = (int)gogeom.size();
      for (int i=0; i<nmbgeom; i++)
	{
	  if (gogeom[i].get() == 0)
	    continue;
	  shared_ptr<GeomObject> lg = gogeom[i];

	  shared_ptr<ParamSurface> sf =
	    dynamic_pointer_cast<ParamSurface, GeomObject>(lg);

	  sf->swapParameterDirection();

	  sf->writeStandardHeader(file2);
	  sf->write(file2);
	}
    }
  else
    {
      shared_ptr<SurfaceModel> shell;

      CompositeModelFileHandler fileread;
      shared_ptr<Body> body = fileread.readBody(file1);
      if (!body.get())
	{
	  shell = fileread.readShell(file1);
	  if (!shell.get())
	    exit(1);
	}
      else
	{
	  shell = body->getOuterShell();
	}

      int nmb = shell->nmbEntities();
      vector<shared_ptr<ParamSurface> > all_sfs;
      for (int ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> surf = shell->getSurface(ki);
	  surf->swapParameterDirection();
	  all_sfs.push_back(surf);
	}

      tpTolerances tptol = shell->getTolerances();
      shared_ptr<SurfaceModel> shell2(new SurfaceModel(tptol.gap, tptol.gap,
						       tptol.neighbour,
						       tptol.kink, tptol.bend,
						       all_sfs));
      CompositeModelFileHandler filewrite;
      filewrite.writeStart(file2);
      filewrite.writeHeader("Turned orientation", file2);
      filewrite.writeSurfModel(*shell2, file2);
      filewrite.writeEnd(file2);
    }
}




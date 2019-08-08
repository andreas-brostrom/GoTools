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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 6) {
    std::cout << "Input parameters : Input file(g2/g22), type of input file, nest stimplify step (0/1), output file, output g22 (0/1)" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::string file1(argv[1]);

  int file_type = atoi(argv[2]);
  int simplify2 = atoi(argv[3]);
  std::ofstream file2(argv[4]);
  bool g22 = atoi(argv[5]);

  shared_ptr<SurfaceModel> sfmodel;
  if (file_type == 1)
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.0001; // 0.001;
      double neighbour = 0.001; // 0.01;
      double kink = 0.01;
      double approxtol = 0.001;
      
      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);
      
      std::ifstream is(file1);
      CompositeModel *model = factory.createFromG2(is);
      
      sfmodel =  
	shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
    }
  else
    {
      CompositeModelFileHandler filehandler;
      sfmodel = filehandler.readShell(file1.c_str());
    }

  if (sfmodel)
    {
      int degree = 3;

      if (simplify2 == 2)
	{
	  try {
	    SurfaceModelUtils::simplifySurfaceModel2(sfmodel, degree);
	  }
	  catch (...)
	    {
	      ;
	    }
	}
      else if (simplify2 == 1)
	{
	  try {
	    SurfaceModelUtils::simplifySurfaceModel(sfmodel, degree);
	  }
	  catch (...)
	    {
	      ;
	    }
	}
      else
	sfmodel->simplifyShell();
      


      if (g22)
	{
	  CompositeModelFileHandler filehandler;
	  filehandler.writeStart(file2);
	  filehandler.writeHeader("Simplified model", file2);
	  filehandler.writeSurfModel(*sfmodel, file2);
	  filehandler.writeEnd(file2);
	}
      else
	{
	  int nmb = sfmodel->nmbEntities();
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);
	    surf->writeStandardHeader(file2);
	    surf->write(file2);
	    }
	}
    }
}


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
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.0001; //0.001;
  double neighbour = 0.001; //0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  vector<shared_ptr<ftSurface> > faces = sfmodel->allFaces();

  // Identify faces with same underlying surfaces
  vector<vector<shared_ptr<ftSurface> > > face_grps;
  vector<shared_ptr<ParamSurface> > under_sfs;
  SurfaceModelUtils::sameUnderlyingSurf(faces, gap, kink, face_grps, under_sfs);

  // Extend underlying surfaces 
  BoundingBox tot = model->boundingBox();
  Point high1 = tot.high();
  Point low1 = tot.low();
  double l1 = high1.dist(low1);
  for (size_t ki=0; ki<under_sfs.size(); ++ki)
    {
      BoundingBox curr = under_sfs[ki]->boundingBox();
      Point high2 = curr.high();
      Point low2 = curr.low();
      double l2 = high2.dist(low2);
      double len = std::max(2.0*l1/l2, l1/4.0);
      ElementarySurface *elem = under_sfs[ki]->elementarySurface();
      if (elem)
	{
	  shared_ptr<ElementarySurface> elem2(elem->clone());
	  elem2->enlarge(len, len, len, len);
	  under_sfs[ki] = elem2;
	}
    }

  // Fetch remaining underlying surfaces
  for (size_t ki=0; ki<faces.size(); ++ki)
    {
      // Check if the current face belongs to a group
      size_t kj, kr;
      for (kj=0; kj<face_grps.size(); ++kj)
	{
	  for (kr=0; kr<face_grps[kj].size(); ++kr)
	    {
	      if (face_grps[kj][kr].get() == faces[ki].get())
		break;
	    }
	  if (kr < face_grps[kj].size())
	    break;
	}
      if (face_grps.size() == 0 || kj == face_grps.size())
	{
	  // Make face group with a single face
	  vector<shared_ptr<ftSurface> > single_face;
	  single_face.push_back(faces[ki]);
	  face_grps.push_back(single_face);
	  shared_ptr<ParamSurface> surf = faces[ki]->surface();
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bd_sf.get())
	    under_sfs.push_back(bd_sf->underlyingSurface());
	  else
	    under_sfs.push_back(surf);
	}
    }

  // Write
  std::ofstream out_file("under_sfs.g2");
  for (size_t ki=0; ki<under_sfs.size(); ++ki)
    {
      under_sfs[ki]->writeStandardHeader(out_file);
      under_sfs[ki]->write(out_file);
    }
}


 

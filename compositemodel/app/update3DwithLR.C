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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: surface in (g2), point cloud (.g2), surface out.g2, tol, maxiter" << std::endl;
    return -1;
  }

  int ki, kj, kh;

  std::ifstream surfin(argv[1]);
  std::ifstream pntin(argv[2]);
  std::ofstream fileout(argv[3]);
  double AEPSGE = atof(argv[4]);
  int max_iter = atoi(argv[5]);

  double knot_diff_tol = 1.0e-10;

  // Read surface
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  ObjectHeader header;
  header.read(surfin);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(surfin);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);

  // Represent as tensor product splines
  shared_ptr<SplineSurface> tmp_spline;
  shared_ptr<SplineSurface> tmp_spline2;
  SplineSurface *base = sf->getSplineSurface();
  if (!base)
    {
      tmp_spline = shared_ptr<SplineSurface>(sf->asSplineSurface());
      base = tmp_spline.get();
    }
  
  if (base->rational())
    {
      // Approximate with non-rational surface
      shared_ptr<SplineSurface> base2(base->clone());
      vector<shared_ptr<ParamSurface> > base_sfs;
      base_sfs.push_back(base2);
      OffsetSurfaceStatus status = 
	OffsetSurfaceUtils::offsetSurfaceSet(base_sfs, 0.0, 0.1*AEPSGE, 
					     tmp_spline2);
      base = tmp_spline2.get();
    }

  if (base->numCoefs_u() == base->order_u() && base->order_u() > 2)
    {
      vector<double> new_knots(2);
      new_knots[0] = (2.0*base->startparam_u() + base->endparam_u())/3.0;
      new_knots[1] = (base->startparam_u() + 2.0*base->endparam_u())/3.0;
      base->insertKnot_u(new_knots);
    }
  if (base->numCoefs_v() == base->order_v() && base->order_v() > 2)
    {
      vector<double> new_knots(2);
      new_knots[0] = (2.0*base->startparam_v() + base->endparam_v())/3.0;
      new_knots[1] = (base->startparam_v() + 2.0*base->endparam_v())/3.0;
      base->insertKnot_v(new_knots);
    }

  // Read parameterized points (u, v, x, y, z)
  int nmb_pts;
  int dim=3, del=5;
  pntin >> nmb_pts;
  vector<double> data(del*nmb_pts);
  for (ki=0; ki<nmb_pts; ++ki)
    pntin >> data[del*ki] >> data[del*ki+1] >> data[del*ki+2] 
	  >> data[del*ki+3] >> data[del*ki+4];

  // For each dimension
  int order = 4; 
  vector<double> uknots(2*order);
  vector<double> vknots(2*order);
  for (ki=0; ki<order; ++ki)
    {
      uknots[ki] = base->startparam_u();
      uknots[order+ki] = base->endparam_u();
      vknots[ki] = base->startparam_v();
      vknots[order+ki] = base->endparam_v();
    }

  vector<shared_ptr<SplineSurface> > diff_sfs(3);

  // Extract initial knot values
  vector<double> init_knots_u, init_knots_v;
  base->basis_u().knotsSimple(init_knots_u);
  base->basis_v().knotsSimple(init_knots_v);

  // Remove first and last knot to keep only the inner ones
  init_knots_u.pop_back();
  init_knots_u.erase(init_knots_u.begin(), init_knots_u.begin()+1);
  init_knots_v.pop_back();
  init_knots_v.erase(init_knots_v.begin(), init_knots_v.begin()+1);

  for (kj=0; kj<dim; ++kj)
    {
      // Fetch data points
      vector<double> data2(3*nmb_pts);
      for (ki=0; ki<nmb_pts; ++ki)
	{
	  data2[3*ki] = data[del*ki];
	  data2[3*ki+1] = data[del*ki+1];
	  data2[3*ki+2] = data[del*ki+kj+2];
	}

      // Create approximation machinery
      LRSurfApprox approx(order, uknots, order, vknots, data2, 1, AEPSGE, 
			  true, 0.0, false, false);
      approx.setFixCorner(false /*true*/);
      approx.setUseMBA(true);

      // Add knot line information from base surface
      if (init_knots_u.size() > 0 || init_knots_v.size() > 0)
	approx.setInitialKnots(init_knots_u, init_knots_v);

      approx.setVerbose(true);

      double maxdist, avdist, avdist_total; // will be set below
      int nmb_out_eps;        // will be set below
      shared_ptr<LRSplineSurface> surf = 
	approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, max_iter);

      std::cout << "No. elements: " << surf->numElements();
      std::cout << ", maxdist= " << maxdist << "avdist= " << avdist_total;
      std::cout << ", avdist(out)= " << avdist;
      std::cout << ", nmb out= " << nmb_out_eps << std::endl;

      diff_sfs[kj] = shared_ptr<SplineSurface>(surf->asSplineSurface());
      int stop_break = 1;
    }

  // Ensure that the difference surface have the same knot vectors
  GeometryTools::unifySurfaceSplineSpace(diff_sfs, knot_diff_tol);

  // Create 3D difference surface
  int ncoef = diff_sfs[0]->numCoefs_u()*diff_sfs[0]->numCoefs_v();
  vector<double> coefs(dim*ncoef);
  for (kj=0; kj<dim; ++kj)
    {
      vector<double>::const_iterator it = diff_sfs[kj]->coefs_begin();
      for (ki=0; it != diff_sfs[kj]->coefs_end(); ++it, ++ki)
	coefs[ki*dim+kj] = (*it);
    }
  shared_ptr<SplineSurface> diff3D(new SplineSurface(diff_sfs[0]->numCoefs_u(),
						     diff_sfs[0]->numCoefs_v(),
						     diff_sfs[0]->order_u(),
						     diff_sfs[0]->order_v(),
						     diff_sfs[0]->basis_u().begin(),
						     diff_sfs[0]->basis_v().begin(),
						     coefs.begin(), dim));

  // Create 3D unified surface
  shared_ptr<SplineSurface> unified_sf = 
    GeometryTools::surfaceSum(*base, 1.0, *diff3D, 1.0, knot_diff_tol);

  unified_sf->writeStandardHeader(fileout);
  unified_sf->write(fileout);
}

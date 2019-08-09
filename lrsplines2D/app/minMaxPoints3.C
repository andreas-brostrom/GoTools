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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <fstream>
#include <iostream>

using namespace Go;
using std::vector;

void extractInnerCurves(vector<vector<shared_ptr<ParamCurve> > >& cvs,
			vector<BoundingBox>& bbox,
			double eps)
{
  size_t ki, kj;
  for (ki=0; ki<cvs.size(); )
    {
      for (kj=ki+1; kj<cvs.size(); )
	{
	  if (!(bbox[ki].overlaps(bbox[kj])))
	    {
	      ++kj;
	      continue;  // No possibility of one curve lying inside the other
	    }

	  double size1 = bbox[ki].low().dist(bbox[ki].high());
	  double size2 = bbox[kj].low().dist(bbox[kj].high());
	  size_t ix1 = (size1 < size2) ? ki : kj;
	  size_t ix2 = (ix1 == kj) ? ki : kj;

	  int ka = 0;
	  for (ka=0; ka<2; ++ka)
	    {
	      // Check if curve ix1 lies inside curve ix2
	      // Compute one point on the curve
	      Point pos = cvs[ix1][0]->point(0.5*(cvs[ix1][0]->startparam()+
						  cvs[ix1][0]->endparam()));
	      shared_ptr<CurveLoop> cvloop(new CurveLoop(cvs[ix2], eps, false));
	      CurveBoundedDomain cvdom(cvloop);

	      Vector2D pos2(pos[0], pos[1]);
	      bool inside = cvdom.isInDomain(pos2, eps);
	      if (inside)
		{
		  break;
		}
	      std::swap(ix1, ix2);
	    }

	  if (ka < 2)
	    {
	      // An "inside" curve is found. Remove the other one
	      cvs.erase(cvs.begin()+ix2);
	      bbox.erase(bbox.begin()+ix2);
	      if (ix2 == ki)
		break;
	    }
	  else
	    ++kj;
	}
      if (kj == cvs.size())
	++ki;
    }
}



int main(int argc, char* argv[] )
{

  if (argc != 5)
    {
      std::cout << "Usage: " " infile curves, infile surface, outfile, tolerance" << std::endl;
      exit(1);
    }

  // Open input files
  std::ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  std::ifstream is2(argv[2]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");

  // Open outfile
  std::ofstream os(argv[3]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  double tol = atof(argv[4]);
  double eps = 1.0e-6;

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ObjectHeader header;
  shared_ptr<GeomObject> geom_obj;
  try {
  header.read(is2);
  geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
  geom_obj->read(is2);
  }
  catch(...)
    {
      std::cout << "WARNING: No surface found" << std::endl;
      exit(0);
    }
  
  shared_ptr<ParamSurface> surf = 
    dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  
  // Fetch curve loop
  CurveLoop loop = SurfaceTools::outerBoundarySfLoop(surf, eps);

  // Fetch parameter curves
  int nmb = loop.size();
  vector<shared_ptr<ParamCurve> > par_loop(nmb);
  for (int ka=0; ka<nmb; ++ka)
    {
      if (loop[ka]->dimension() == 2)
	par_loop[ka] = loop[ka];
      else
	{
	  shared_ptr<CurveOnSurface> sf_cv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loop[ka]);
	  if (!sf_cv.get())
	    {
	      std::cout << "Parameter loop not available" << std::endl;
	      exit(1);
	    }

	  if (!sf_cv->hasParameterCurve())
	    {
	      sf_cv->ensureParCrvExistence(tol);
	    }
	  par_loop[ka] = sf_cv->parameterCurve();
	  if (!par_loop[ka].get())
	    {
	      std::cout << "Parameter loop not available" << std::endl;
	      exit(1);
	    }
	}
    }

  // Read curves
  vector<vector<shared_ptr<ParamCurve> > > cvs;
  vector<BoundingBox> bbox;
  while (!is.eof())
    {
      // Read curves from file
      shared_ptr<SplineCurve> cv(new SplineCurve());
      is >> header;
      cv->read(is);

      if (cv->dimension() != 2)
	{
	  std::cout << "Curve dimension different from 2" << std::endl;
	  exit(1);
	}

      // Compute distance between curve endpoints
      Point pos1 = cv->ParamCurve::point(cv->startparam());
      Point pos2 = cv->ParamCurve::point(cv->endparam());

      // Assign bounding boxes
      BoundingBox bb = cv->boundingBox();
      vector<shared_ptr<ParamCurve> > contour_loop;
      contour_loop.push_back(cv);
      double dist0 = pos1.dist(pos2);
      if (dist0 > eps)
	{
	  // Close the curve loop
	  int ind1, ind2;
	  double par1, par2, dist1, dist2;
	  Point pt1, pt2;
	  loop.closestParPoint(pos1, ind1, par1, pt1, dist1);
	  loop.closestParPoint(pos2, ind2, par2, pt2, dist2);

	  if (dist0 < std::min(dist1, dist2))
	    {
	      // Close original loop
	      shared_ptr<ParamCurve> cv2(new SplineCurve(pos1, pos2));
	      BoundingBox bb2 = cv2->boundingBox();
	      contour_loop.push_back(cv2);
	      bb.addUnionWith(bb2);
	      cvs.push_back(contour_loop);
	      bbox.push_back(bb);
	    }
	  else
	    {
	      // Include part of surface loop in the curve loop
	      // Compute curve lengths
	      double len0 = cv->estimatedCurveLength();
	      int nmb = loop.size();
	      if (ind1 > ind2 || (ind1 == ind2 && par1 > par2))
		{
		  std::swap(ind1, ind2);
		  std::swap(par1, par2);
		}
	      double len1 = 
		par_loop[ind1]->estimatedCurveLength(par1, 
						     (ind2 != ind1) ? 
						     loop[ind1]->endparam() : par2);
	      for (int ka=ind1+1; ka<ind2; ++ka)
		len1 += par_loop[ka]->estimatedCurveLength();
	      if (ind1 != ind2)
		len1 += par_loop[ind2]->estimatedCurveLength(loop[ind2]->startparam(),
							     par2);

	      double len2 = 
		par_loop[ind2]->estimatedCurveLength(par2, loop[ind2]->endparam());
	      for (int ka=(ind2+1)%nmb; ka!=ind1; ka=(ka+1)%nmb)
		len2 += par_loop[ka]->estimatedCurveLength();
	      len2 += par_loop[ind1]->estimatedCurveLength(loop[ind1]->startparam(),
							   par1);
	      
	      vector<shared_ptr<ParamCurve> > contour_loop2;
	      contour_loop2.push_back(cv);
	      BoundingBox bb3 = bb;
	      if (len2 >= len0 + len1)
		{
		  for (int ka=ind1; ka<=ind2; ++ka)
		    {
		      double t1 = (ka == ind1) ? par1 : loop[ka]->startparam();
		      double t2 = (ka == ind2) ? par2 : loop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> cv2 = (ka != ind1 && ka != ind2) ?
			par_loop[ka] : 
			shared_ptr<ParamCurve>(par_loop[ka]->subCurve(t1, t2));
		      BoundingBox bb2 = cv2->boundingBox();
		      contour_loop.push_back(cv2);
		      bb.addUnionWith(bb2);
		    }
		  cvs.push_back(contour_loop);
		  bbox.push_back(bb);
		}

	      if (len1 >= len0 + len2)
		{
		  int ka, kb;
		  int nmb2 = nmb - (ind2-ind1);
		  for (ka=ind2, kb=0; kb<=nmb2; ka=(ka+1)%nmb, ++kb)
		    {
		      double t1 = (ka == ind2) ? par2 : loop[ka]->startparam();
		      double t2 = (ka == ind1) ? par1 : loop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> cv2 = (ka != ind1 && ka != ind2) ?
			par_loop[ka] : 
			shared_ptr<ParamCurve>(par_loop[ka]->subCurve(t1, t2));
		      BoundingBox bb2 = cv2->boundingBox();
		      contour_loop2.push_back(cv2);
		      bb3.addUnionWith(bb2);
		    }
		  cvs.push_back(contour_loop2);
		  bbox.push_back(bb3);
		}
	    }
	}
      else
	{
	  cvs.push_back(contour_loop);
	  bbox.push_back(bb);
	}

      Utils::eatwhite(is);
    }
  
  // Ensure consistent direction of curves in loop
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      if (cvs[ki].size() == 1)
	continue;

      Point pos1 = cvs[ki][0]->point(cvs[ki][0]->startparam());
      Point pos2 = cvs[ki][0]->point(cvs[ki][0]->endparam());
      Point pos3 = cvs[ki][1]->point(cvs[ki][1]->startparam());
      if (pos3.dist(pos1) < pos3.dist(pos2))
	cvs[ki][0]->reverseParameterDirection();
    }

  // For each combination of curves, check if one curve can possibly lie 
  // inside the other. If so, check if it is the case. Curves that are
  // external to another curve is removed.
  extractInnerCurves(cvs, bbox, eps);

  // Write result to file
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      for (size_t kj=0; kj<cvs[ki].size(); ++kj)
	{
	  cvs[ki][kj]->writeStandardHeader(os);
	  cvs[ki][kj]->write(os);
	}
    }

}


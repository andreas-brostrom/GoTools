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

#ifndef _CUTCELLQUAL_H_
#define  _CUTCELLQUAL_H_

#include "GoTools/utils/Point.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveBoundedDomain.h"

namespace Go
{
  class CutCellQuad
  {
  public:
    // Constructor
    CutCellQuad(std::vector<shared_ptr<ParamCurve> >& bd_curves,
		double tol);

    // Define quadrature information
    void setQuadratureInfo(std::vector<double>& quadpar,
			   double min_cell_size)
    {
      quadpar_ = quadpar;
      min_cell_size_ = min_cell_size;
    }
    
    // Check cell status
    int cellStat(const Point& ll, const Point& ur);

    // Compute quadratures
    void quadrature(const Point& ll, const Point& ur,
		    std::vector<std::vector<double> >& quadraturepoints,
		    std::vector<std::vector<shared_ptr<ParamCurve> > >& unresolved_cells,
		    int stat = -1);

  private:
    //shared_ptr<CurveBoundedDomain> domain_;
    double tol_;
    double angtol_;
    std::vector<shared_ptr<CurveLoop> > loops_;
    std::vector<double> quadpar_;
    double min_cell_size_;

    // Sort boundary curves
    void sortBoundary(std::vector<shared_ptr<ParamCurve> >& bd_curves,
		      std::vector<shared_ptr<CurveLoop> >& bd_loops);

    void
      createCutCell(std::vector<shared_ptr<CurveLoop> >& bd_loops,
		    std::vector<shared_ptr<SplineCurve> >& cell_cvs,
		    std::vector<std::vector<shared_ptr<CurveLoop> > >& trim_loops,
    		    bool test_inside);

    
    void
      sortLoops(std::vector<std::vector<shared_ptr<ParamCurve> > >& loop_cvs,
		std::vector<shared_ptr<ParamCurve> >& cell_cvs,
		std::vector<std::vector<shared_ptr<CurveLoop> > >& trim_loops,
		bool test_inside);

    void
      quadraturePoints(std::vector<shared_ptr<CurveLoop> >& cell_loops,
		       std::vector<std::vector<double> >& quadraturepoints,
		       std::vector<std::vector<shared_ptr<ParamCurve> > >& unresolved_cells);
    
    int heightDirection(shared_ptr<CurveLoop> cell_loop,
			const RectDomain& domain);
    
    void
      computeQuadraturePoints(shared_ptr<CurveBoundedDomain> cvdom,
			      const RectDomain& domain,
			      int dir,
			      std::vector<double>& quadraturepoints);
    
    void splitPars(std::vector<shared_ptr<CurveLoop> >& cell_loops,
		   shared_ptr<CurveBoundedDomain> cvdom,
		   const RectDomain& domain, std::vector<double>& splitpar,
		   int& dir);
    
    void defineSplits1(const std::vector<Point>& corner,
		       shared_ptr<CurveBoundedDomain> cvdom,
		       const RectDomain& domain,
		       std::vector<double>& splitpar, int& dir,
		       std::vector<double>& candpar);

    void defineSplits2(const std::vector<Point>& turnpts,
		       shared_ptr<CurveBoundedDomain> cvdom,
		       const RectDomain& domain,
		       std::vector<double>& splitpar, int& dir,
		       std::vector<double>& candpar);

    void fetchTurningPoints(shared_ptr<ParamCurve> cv,
			    std::vector<Point>& turnpts);

    bool checkSplits(std::vector<double>& splitpar);

  };
} // end namespace Go

#endif

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

#ifndef _MODELUPDATE_H
#define _MODELUPDATE_H

#include "GoTools/utils/Point.h"
#include <vector>

namespace Go
{

  /// Modify brep model with respect to point information

  class SurfaceModel;
  class ParamSurface;
  class SplineSurface;
  class CurveOnSurface;
  class ftEdgeBase;
  class ftEdge;
  class ftSurface;

  class ModelPointData
  {
    friend class ModelUpdate;

  public:
    ModelPointData(Point& pos)
    {
      pos_ = pos;
      has_dist_ = false;
      sf_ix_ = -1;
    }

    ModelPointData(Point& pos, double dist)
    {
      pos_ = pos;
      has_dist_ = true;
      dist_ = dist;
      sf_ix_ = -1;
    }

    Point getPos()
    {
      return pos_;
    }

    Point getPar()
    {
      return Point(upar_, vpar_);
    }

    double surfIx()
    {
      return sf_ix_;
    }
    
    bool hasDist()
    {
      return has_dist_;
    }

    double getDist()
    {
      return dist_;
    }

    void setDist(double dist)
    {
      dist_ = dist;
    }

    void setParameterInfo(int sf_ix, double upar, double vpar)
    {
      sf_ix_ = sf_ix;
      upar_ = upar;
      vpar_ = vpar;
    }

  private:
    Point pos_;
    bool has_dist_;
    double dist_;
    int sf_ix_;
    double upar_, vpar_;

  };

  struct SurfPointInfo
  {
    int nmb_pts_;
    double max_dist_;
    double min_dist_;
    double av_dist_;

    SurfPointInfo()
    {
      nmb_pts_ = 0;
      max_dist_ = -HUGE;
      min_dist_ = HUGE;
      av_dist_ = 0.0;
    }
  };

class ModelUpdate
{
 public:
  /// Constructor
  ModelUpdate(shared_ptr<SurfaceModel> sfmodel);

  /// Destructor
  ~ModelUpdate();

  /// Set point information. 
  void setPointSet(std::vector<Point>& point_set);

  /// Set point information including a distance to be added in
  /// the normal direction
  void setPointSet(std::vector<std::pair<Point,double> >& point_set,
		   bool negate, bool remove=true);

  /// Set maximum iteration level in surface update
  void setMaxIter(int max_iter)
  {
    max_iter_ = max_iter;
  }

  /// Update the model based on given point information
  shared_ptr<SurfaceModel> updateModel();

  /// Update the model based on input information
  shared_ptr<SurfaceModel> updateModel(std::vector<ModelPointData>& data, 
				       bool negate, bool remove=true);

 private:

  shared_ptr<SurfaceModel> model_;
  int max_iter_;
  std::vector<std::pair<double,double> > sf_size_;
  double sf_side_mean_;
  std::vector<ModelPointData> points_;
  std::vector<SurfPointInfo> sf_info_;
  bool negate_;
  bool remove_;
  std::vector<shared_ptr<ParamSurface> > mod_sfs_;
  std::vector<ftEdgeBase*> smooth_edgs_;
  std::vector<ftEdgeBase*> smooth_edgs2_;
  std::vector<shared_ptr<ParamSurface> > trim_sfs_;
  double maxdist_, avdist_;
  
  // DEBUG storage
  std::vector<Point> all_points_;

  void doUpdateModel();

  void projectPoints();

  void updateSingleSurf(int sf_ix, Point& midpar, shared_ptr<ParamSurface>& base);

  void trimWithAdjacent(std::vector<Point>& midpar);

  shared_ptr<ParamSurface> getUpdatedSurf(shared_ptr<ParamSurface> surf,
					  shared_ptr<ftSurface> face,
					  std::vector<Point>& par,
					  std::vector<Point>& diff,
					  double max_dist,
					  shared_ptr<ParamSurface>& base_sf);

  void splitCommonBase(std::vector<int>& same_face_ix,
		       shared_ptr<ParamSurface> base_sf,
		       std::pair<bool,bool>& base_closed);

  double meanSurfaceSide();

  shared_ptr<SplineSurface> updateSurfaceSize(SplineSurface *surf,
					      shared_ptr<ftSurface> face);

  void
    identifyCommonBaseSfs(std::vector<std::vector<shared_ptr<ftSurface> > >& same_faces,
			  std::vector<std::vector<int> >& same_faces_ix,
			  std::vector<shared_ptr<ParamSurface> >& base_sfs,
			  std::vector<std::pair<bool, bool> >& base_closed);

  shared_ptr<ParamSurface>
    constructCommonBase(std::vector<shared_ptr<ftSurface> >& same_faces);

  void 
    fetchAdjacent(int face_ix, std::vector<int>& adj_ix,
		  std::vector<std::pair<int, shared_ptr<ftEdge> > >& smooth,
		  std::vector<std::pair<int, int> >& smooth_ix);

  void 
    splitIntersectionCurves(std::vector<int>& ix2,
			    std::vector<std::pair<shared_ptr<CurveOnSurface>,
			    shared_ptr<CurveOnSurface> > >& int_seg,
			    std::vector<int>& int_type,
			    int& nmb_curr, double epsge, double minlen);
  void 
   cleanIntersectionCurves(shared_ptr<ParamSurface> surf,
			   std::vector<int>& ix2,
			   std::vector<std::pair<shared_ptr<CurveOnSurface>,
			   shared_ptr<CurveOnSurface> > >& int_seg,
			   double epsge, double minlen, double minlen2);

  void getTangentialInt(shared_ptr<ParamSurface> surf1, 
			shared_ptr<ParamSurface> surf2, 
			std::vector<shared_ptr<CurveOnSurface> >& int_segments1, 
			std::vector<shared_ptr<CurveOnSurface> >& int_segments2, 
			double epsge);
};
}

#endif // _MODELUPDATE_H

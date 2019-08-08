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

#include "GoTools/compositemodel/ModelUpdate.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/topology/FaceConnectivityUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/Domain.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "sisl.h"
#include <fstream>

#define DEBUG

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;

//==========================================================================
ModelUpdate::ModelUpdate(shared_ptr<SurfaceModel> sfmodel)
  : max_iter_(2), negate_(false), remove_(true), maxdist_(0.0), avdist_(0.0)
//==========================================================================
{
  model_ = sfmodel;

  // Compute mean estimated surface size
  sf_side_mean_ = meanSurfaceSide();
}

//==========================================================================
ModelUpdate::~ModelUpdate()
//==========================================================================
{
}

//==========================================================================
void ModelUpdate::setPointSet(vector<Point>& point_set)
//==========================================================================
{
  //points_.resize(point_set.size());
  points_.reserve(point_set.size());
  for (size_t ki=0; ki<point_set.size(); ++ki)
    points_.push_back(ModelPointData(point_set[ki]));
  //points_[ki] = ModelPointData(point_set[ki]);
}

//==========================================================================
void ModelUpdate::setPointSet(vector<std::pair<Point,double> >& point_set,
			      bool negate, bool remove)
//==========================================================================
{
  //points_.resize(point_set.size());
  points_.reserve(point_set.size());
  for (size_t ki=0; ki<point_set.size(); ++ki)
    points_.push_back(ModelPointData(point_set[ki].first, point_set[ki].second));
  //points_[ki] = ModelPointData(point_set[ki].first, point_set[ki].second);
  negate_ = negate;
  remove_ = remove;
}

//==========================================================================
shared_ptr<SurfaceModel> ModelUpdate::updateModel()
//==========================================================================
{
  // Parameterize point set based by projecting the points onto the model
  projectPoints();

  doUpdateModel();

  tpTolerances tptol = model_->getTolerances();
  shared_ptr<SurfaceModel> updated_model(new SurfaceModel(tptol.gap, tptol.gap,
							  tptol.neighbour,
							  tptol.kink,
							  tptol.bend, 
							  trim_sfs_));
#ifdef DEBUG
  if (updated_model->nmbBoundaries() > 0)
    {
      std::cout << "Open shell." << std::endl;
      vector<shared_ptr<ftEdge> > bd_edgs = 
	updated_model->getBoundaryEdges();
      std::ofstream of_bd("updated_nontwin_bd.g2");
      for (size_t kr=0; kr<bd_edgs.size(); ++kr)
	{
	  shared_ptr<ParamCurve> tmp_cv = bd_edgs[kr]->geomCurve();
	  shared_ptr<ParamCurve> tmp_cv2(tmp_cv->geometryCurve());
	  tmp_cv2->writeStandardHeader(of_bd);
	  tmp_cv2->write(of_bd);
	}
    }
#endif

  return updated_model;
}

//==========================================================================
shared_ptr<SurfaceModel> 
ModelUpdate::updateModel(vector<ModelPointData>& data, bool negate,
			 bool remove)
//==========================================================================
{
  // Store input information
  points_ = data;
  negate_ = negate;
  remove_ = remove;

  int nmb_sfs = model_->nmbEntities();
  int nmb_points = (int)points_.size();
  sf_info_.resize(nmb_sfs);
  for (int ka=0; ka<nmb_points; ++ka)
    {
      int sf_ix = points_[ka].surfIx();
      double pt_dist = points_[ka].hasDist() ? points_[ka].getDist() : 0.0;
      if (remove_ == false)
	{
	  if (negate_ && pt_dist > 0.0)
	    pt_dist = 0.0;
	  else if (negate_ == false && pt_dist < 0)
	    pt_dist = 0.0;
	}
      sf_info_[sf_ix].nmb_pts_++;
      sf_info_[sf_ix].max_dist_ = std::max(sf_info_[sf_ix].max_dist_, pt_dist);
      sf_info_[sf_ix].min_dist_ = std::min(sf_info_[sf_ix].min_dist_, pt_dist);
      sf_info_[sf_ix].av_dist_ += fabs(pt_dist);
    }  
  for (int ka=0; ka<nmb_sfs; ++ka)
    {
      if (sf_info_[ka].nmb_pts_ > 0)
	sf_info_[ka].av_dist_ /= (double)sf_info_[ka].nmb_pts_;
    }

#ifdef DEBUG
  std::ofstream info("input_info.txt");
  double maxtotal = -HUGE;
  double mintotal = HUGE;
  double avtotal = 0.0;
  for (size_t ki=0; ki<sf_info_.size(); ++ki)
    {
      info << ki << " " << sf_info_[ki].nmb_pts_ << " " << sf_info_[ki].max_dist_
	   << " " << sf_info_[ki].min_dist_ << " " << sf_info_[ki].av_dist_ << std::endl;
      maxtotal = std::max(maxtotal, sf_info_[ki].max_dist_);
      mintotal = std::min(mintotal, sf_info_[ki].min_dist_);
      avtotal += (sf_info_[ki].nmb_pts_*sf_info_[ki].av_dist_);
    }
  avtotal /= (double)nmb_points;

  info << "Number of points: " << nmb_points << std::endl;
  info << "Maximum positive distance: " << maxtotal << std::endl;
  info << "Maximum negative distance: " << mintotal << std::endl;
  info << "Average distance: " << avtotal << std::endl;
#endif

  doUpdateModel();

  tpTolerances tptol = model_->getTolerances();
  shared_ptr<SurfaceModel> updated_model(new SurfaceModel(tptol.gap, tptol.gap,
							  tptol.neighbour,
							  tptol.kink,
							  tptol.bend, 
							  trim_sfs_));
#ifdef DEBUG
  if (updated_model->nmbBoundaries() > 0)
    {
      std::cout << "Open shell." << std::endl;
      vector<shared_ptr<ftEdge> > bd_edgs = 
	updated_model->getBoundaryEdges();
      std::ofstream of_bd("updated_nontwin_bd.g2");
      for (size_t kr=0; kr<bd_edgs.size(); ++kr)
	{
	  shared_ptr<ParamCurve> tmp_cv = bd_edgs[kr]->geomCurve();
	  shared_ptr<ParamCurve> tmp_cv2(tmp_cv->geometryCurve());
	  tmp_cv2->writeStandardHeader(of_bd);
	  tmp_cv2->write(of_bd);
	}
    }
#endif

  return updated_model;
}

//==========================================================================
void ModelUpdate::doUpdateModel()
//==========================================================================
{
  // Create storage for updated surfaces
  int nmb = model_->nmbEntities();
  mod_sfs_.resize(nmb);
  tpTolerances tptol = model_->getTolerances();

#ifdef DEBUG
  std::ofstream of("mod_sfs.g2");
#endif

  // Identify surfaces having the same underlying base surface
  vector<vector<shared_ptr<ftSurface> > > same_faces;
  vector<shared_ptr<ParamSurface> > base_sfs;
  vector<vector<int> > same_faces_ix;
  vector<pair<bool, bool> > base_closed;
  identifyCommonBaseSfs(same_faces, same_faces_ix, base_sfs, base_closed);
 
  // Update each surface
  vector<Point> midparam(nmb);
  maxdist_ = -HUGE;
  avdist_ = 0.0;
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> curr_base;
      int curr_ix = -1;
      for (size_t kr=0; kr<same_faces_ix.size(); ++kr)
	{
	  for (size_t kh=0; kh<same_faces_ix[kr].size(); ++kh)
	    if (same_faces_ix[kr][kh] == ki)
	      {
		curr_base = base_sfs[kr];
		curr_ix = (int)kr;
		break;
	      }
	  if (curr_ix >= 0)
	    break;
	}

      Point midpar;
      updateSingleSurf(ki, midpar, curr_base);
      midparam[ki] = midpar;
      if (curr_ix >= 0)
	{
	  for (int ka=0; ka<ki; ++ka)
	    if (mod_sfs_[ka].get() == base_sfs[curr_ix].get())
	      mod_sfs_[ka] = curr_base;

	  base_sfs[curr_ix] = curr_base;
	}

#ifdef DEBUG
      if (all_points_.size() > 0)
	{
	  std::ofstream of1("curr_data_pts.g2");
	  of1 << "400 1 0 4 255 0 0 255" << std::endl;
	  of1 << all_points_.size() << std::endl;
	  for (size_t kr=0; kr<all_points_.size(); ++kr)
	    of1 << all_points_[kr] << std::endl;
	}

      if (mod_sfs_[ki].get())
	{
	  mod_sfs_[ki]->writeStandardHeader(of);
	  mod_sfs_[ki]->write(of);
	}
#endif
    }
  avdist_ /= points_.size();

#ifdef DEBUG
  if (all_points_.size() > 0)
    {
      std::ofstream of2("data_pts.g2");
      of2 << "400 1 0 4 255 0 0 255" << std::endl;
      of2 << all_points_.size() << std::endl;
      for (size_t kr=0; kr<all_points_.size(); ++kr)
	of2 << all_points_[kr] << std::endl;
    }
#endif

#ifdef DEBUG
  std::ofstream info("approx_info.txt");
  info << "Number of points: " << points_.size() << std::endl;
  info << "Maximum distance: " << maxdist_ << std::endl;
  info << "Average distance: " << avdist_ << std::endl;
#endif

  for (size_t kr=0; kr<base_sfs.size(); ++kr)
    {
      if (base_closed[kr].first || base_closed[kr].second)
	{
	  shared_ptr<SplineSurface> base(base_sfs[kr]->asSplineSurface());
	  if (base_closed[kr].first)
	    {
	      Point dummy;
	      GeometryTools::averageCoefsAtSeam(base, 1, true);
	    }
	  if (base_closed[kr].second)
	    {
	      Point dummy;
	      GeometryTools::averageCoefsAtSeam(base, 2, true);
	    }
	  base_sfs[kr] = base;
	  for (size_t kh=0; kh<same_faces_ix[kr].size(); ++kh)
	    mod_sfs_[same_faces_ix[kr][kh]] = base;
	}
    }
  // //Split common base surfaces according to face information and
  // //ensure continuity across the seam
  // for (size_t kr=0; kr<same_faces_ix.size(); ++kr)
  //   splitCommonBase(same_faces_ix[kr], base_sfs[kr], base_closed[kr]);

#ifdef DEBUG
  std::ofstream pre_trim("pre_trim.g2");
  for (size_t kr=0; kr<mod_sfs_.size(); ++kr)
    {
      mod_sfs_[kr]->writeStandardHeader(pre_trim);
      mod_sfs_[kr]->write(pre_trim);
    }
#endif

  trimWithAdjacent(midparam);

  // Check for existence of trimmed surfaces and split closed surfaces
  for (size_t kr=0; kr<trim_sfs_.size(); )
    {
      if (!trim_sfs_[kr].get())
	trim_sfs_.erase(trim_sfs_.begin()+kr);
      else
	++kr;
    }

  int stop_break = 1;
}

//==========================================================================
void ModelUpdate::projectPoints()
//==========================================================================
{
#ifdef DEBUG
  vector<Point> proj_pts;
  std::ofstream of0("closest_point.txt");
#endif

  int nmb_points = (int)points_.size();
  bool use_surface_model = true; //false;
  
  if (use_surface_model)
    {
      // Project the points onto the model and remember surface identification
      // and parameter values
      int nmb_sfs = model_->nmbEntities();
      sf_info_.resize(nmb_sfs);
     for (size_t ki=0; ki<points_.size(); ++ki)
	{
	  Point clo_pnt;
	  double clo_par[2];
	  double dist;
	  int sf_ix;
	  Point pos = points_[ki].getPos();
	  model_->closestPoint(pos, clo_pnt, sf_ix,
			       clo_par, dist);
	  points_[ki].setParameterInfo(sf_ix, clo_par[0], clo_par[1]);

	  double pt_dist = points_[ki].hasDist() ? points_[ki].getDist() : 0.0;
	  if (remove_ == false)
	    {
	      if (negate_ && pt_dist > 0.0)
		pt_dist = 0.0;
	      else if (negate_ == false && pt_dist < 0)
		pt_dist = 0.0;
	    }
	  sf_info_[sf_ix].nmb_pts_++;
	  sf_info_[sf_ix].max_dist_ = std::max(sf_info_[sf_ix].max_dist_, pt_dist);
	  sf_info_[sf_ix].min_dist_ = std::min(sf_info_[sf_ix].min_dist_, pt_dist);
	  sf_info_[sf_ix].av_dist_ += pt_dist;
#ifdef DEBUG
	  of0 << sf_ix << " " << clo_par[0] << " " << clo_par[1] << " " 
	      << clo_pnt << " " << points_[ki].getDist() << std::endl;
	  proj_pts.push_back(clo_pnt);
#endif
	  int stop_break = 1;
	}

      for (int ka=0; ka<nmb_sfs; ++ka)
	{
	  if (sf_info_[ka].nmb_pts_ > 0)
	    sf_info_[ka].av_dist_ /= (double)sf_info_[ka].nmb_pts_;
	}
    }
  else
    {
      // Project pointset using registration functionality
      int nmb_sfs = model_->nmbEntities();
      vector<shared_ptr<GeomObject> > surfaces(nmb_sfs);
      for (int ka=0; ka<nmb_sfs; ++ka)
	surfaces[ka] = model_->getSurface(ka);

      // Preparations
      shared_ptr<boxStructuring::BoundingBoxStructure> structure = 
	preProcessClosestVectors(surfaces, 200.0);//, &of_status_filename);
      
      vector<vector<double> > rotation(3);
      for (int ka=0; ka<3; ++ka)
	{
	  rotation[ka].resize(3, 0.0);
	  rotation[ka][ka] = 1.0;
	}
      Point translation(0.0, 0.0, 0.0);
      pair<vector<vector<double> >, Point> currTransformation(rotation,
							      translation);

      // Assemble points
      vector<float> pts(3*points_.size());
      for (int kb=0; kb<nmb_points; ++kb)
	{
	  Point curr = points_[kb].getPos();
	  for (int ka=0; ka<3; ++ka)
	    pts[3*kb+ka] = curr[ka];
	}

      // Closest point computations
      vector<float> signed_dists;
      signed_dists = closestSignedDistanceSfParams(pts, structure,
						   currTransformation.first, 
						   currTransformation.second);

      // Extract results and compute statistic
      
      int result_size = (int)signed_dists.size()/nmb_points;
      if (result_size != 4)
	THROW("Projection failed");

      sf_info_.resize(nmb_sfs);
      for (int ka=0; ka<nmb_points; ++ka)
	{
	  int sf_ix = (int)floor(signed_dists[4*ka+1]+0.1);
	  double upar = signed_dists[4*ka+2];
	  double vpar = signed_dists[4*ka+3];
	  points_[ka].setParameterInfo(sf_ix, upar, vpar);

	  double pt_dist = points_[ka].hasDist() ? points_[ka].getDist() : 0.0;
	  sf_info_[sf_ix].nmb_pts_++;
	  sf_info_[sf_ix].max_dist_ = std::max(sf_info_[sf_ix].max_dist_, pt_dist);
	  sf_info_[sf_ix].min_dist_ = std::min(sf_info_[sf_ix].min_dist_, pt_dist);
	  sf_info_[sf_ix].av_dist_ += pt_dist;

#ifdef DEBUG
	  Point clo_pnt = model_->getSurface(sf_ix)->point(upar, vpar);
	  of0 << sf_ix << " " << upar << " " << vpar << " " 
	      << clo_pnt << " " << points_[ka].getDist() << std::endl;
	  proj_pts.push_back(clo_pnt);
#endif
	}

      for (int ka=0; ka<nmb_sfs; ++ka)
	{
	  if (sf_info_[ka].nmb_pts_ > 0)
	    sf_info_[ka].av_dist_ /= (double)sf_info_[ka].nmb_pts_;
	}
    }

#ifdef DEBUG
  std::ofstream of("proj_points.g2");
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << proj_pts.size() << std::endl;
  for (size_t ki=0; ki<proj_pts.size(); ++ki)
    of << proj_pts[ki] << std::endl;
#endif
#ifdef DEBUG
  std::ofstream info("input_info.txt");
  double maxtotal = -HUGE;
  double mintotal = HUGE;
  double avtotal = 0.0;
  for (size_t ki=0; ki<sf_info_.size(); ++ki)
    {
      info << ki << " " << sf_info_[ki].nmb_pts_ << " " << sf_info_[ki].max_dist_
	   << " " << sf_info_[ki].min_dist_ << " " << sf_info_[ki].av_dist_ << std::endl;
      maxtotal = std::max(maxtotal, sf_info_[ki].max_dist_);
      mintotal = std::min(mintotal, sf_info_[ki].min_dist_);
      avtotal += (sf_info_[ki].nmb_pts_*sf_info_[ki].av_dist_);
    }
  avtotal /= (double)nmb_points;

  info << "Number of points: " << nmb_points << std::endl;
  info << "Maximum positive distance: " << maxtotal << std::endl;
  info << "Maximum negative distance: " << mintotal << std::endl;
  info << "Average distance: " << avtotal << std::endl;
#endif

  int stop_break = 2;
}

//==========================================================================
 void ModelUpdate::updateSingleSurf(int sf_ix, Point& midpar,
				    shared_ptr<ParamSurface>& base)
//==========================================================================
{
  // Fetch surface
  shared_ptr<ParamSurface> surf = model_->getSurface(sf_ix);

  // Fetch internal point in surface
  double u_in, v_in;
  Point pt_in = surf->getInternalPoint(u_in, v_in);
  Point par_in = Point(u_in, v_in);
  Point midpoint = pt_in;
  midpar = par_in;

  // Extract corresponding point data and identify point closest
  // to the internal surface point (in the parameter domain)
  vector<Point> pos;
  vector<Point> diff;
  vector<Point> par;
  double par_dist = HUGE;
  double max_dist = 0.0;
  for (size_t ki=0; ki<points_.size(); ++ki)
    {
      if (points_[ki].surfIx() == sf_ix)
	{
	  Point sfpar = points_[ki].getPar();
	  par.push_back(sfpar);
	  if (sfpar.dist(par_in) < par_dist)
	    {
	      par_dist = sfpar.dist(par_in);
	      midpar = sfpar;
	      midpoint = surf->point(sfpar[0], sfpar[1]);
	    }

	  int sgn = (negate_) ? -1 : 1;
	  if (points_[ki].hasDist())
	    {
	      double dist = points_[ki].getDist();
	      vector<Point> sfpts;
	      if (base.get())
		sfpts = base->point(sfpar[0], sfpar[1], 1);
	      else
		sfpts = surf->point(sfpar[0], sfpar[1], 1);
	      Point sfnorm = sfpts[1]%sfpts[2];
	      double len = sfnorm.normalize_checked();
	      Point pt_diff;
	      if (remove_ || sgn*dist > 0.0)
		pt_diff = sgn*dist*sfnorm;
	      else
		pt_diff = Point(0.0, 0.0, 0.0);
	      Point offpt = sfpts[0] + pt_diff;
	      pos.push_back(offpt);
	      diff.push_back(offpt - sfpts[0]);
	      max_dist = std::max(max_dist, offpt.dist(sfpts[0]));
	    }
	  else
	    {
	      Point tmp_pt = points_[ki].getPos();
	      pos.push_back(tmp_pt);
	      Point sfpt;
	      if (base.get())
		sfpt = base->point(sfpar[0], sfpar[1]);
	      else
		sfpt = surf->point(sfpar[0], sfpar[1]);
	      diff.push_back(tmp_pt - sfpt);
	      max_dist = std::max(max_dist, tmp_pt.dist(sfpt));
	    }
	}
    }

  if (base.get())
    {
      // Project point onto corresponding base surface
      double u_close, v_close, d_close;
      Point close;
      base->closestPoint(midpoint, u_close, v_close, close, d_close,
			 model_->getTolerances().gap, NULL, midpar.begin());
      midpar = Point(u_close, v_close);
    }


#ifdef DEBUG
  std::ofstream of1("off_pts.g2");
  of1 << "400 1 0 4 255 0 0 255" << std::endl;
  of1 << pos.size() << std::endl;
  for (size_t ki=0; ki<pos.size(); ++ki)
    of1 << pos[ki] << std::endl;

  std::ofstream of2("curr_sf.g2");
  surf->writeStandardHeader(of2);
  surf->write(of2);
#endif

#ifdef DEBUG
  all_points_.insert(all_points_.end(), pos.begin(), pos.end());
#endif

  std::ofstream of3("parpts.raw");
  for (size_t ki=0; ki<diff.size(); ++ki)
    of3 << par[ki] << " " << diff[ki] << std::endl;
  
  mod_sfs_[sf_ix] = getUpdatedSurf(surf, model_->getFace(sf_ix),
				   par, diff, max_dist, base);
  // vector<double> data;
  // data.reserve(5*pos.size());
  // for (size_t ki=0; ki<pos.size(); ++ki)
  //   {
  //     data.insert(data.end(), par[ki].begin(), par[ki].end());
  //     data.insert(data.end(), diff[ki].begin(), diff[ki].end());
  //   }
#ifdef DEBUG
  if (mod_sfs_[sf_ix].get())
    {
      std::ofstream of3("curr_mod_sf.g2");
      mod_sfs_[sf_ix]->writeStandardHeader(of3);
      mod_sfs_[sf_ix]->write(of3);
    }
#endif
  int stop_break = 1;
}

//==========================================================================
void ModelUpdate::trimWithAdjacent(vector<Point>& midpar)
//==========================================================================
{
#ifdef DEBUG
  std::ofstream of("all_trim_sfs.g2");
#endif

  int nmb = model_->nmbEntities();
  trim_sfs_.resize(nmb);

  // For each surface, assemble information about adjacent surfaces
  vector<vector<int> > adjacent(nmb);
  vector<vector<pair<int, shared_ptr<ftEdge> > > > smoothedges(nmb);
  vector<pair<int, int> > adj_ix_smooth;
  for (int ki=0; ki<nmb; ++ki)
    {
      // Check if the surface is trimmed already
      int ka;
      for (ka=0; ka<ki; ++ka)
	if (mod_sfs_[ka].get() == mod_sfs_[ki].get())
	  break;
      if (ka < ki)
	continue;  // Treated as same underlying surface

      fetchAdjacent(ki, adjacent[ki], smoothedges[ki],
		    adj_ix_smooth);
    }

  // Prioritize surfaces such that the one with fewest neighbours are
  // trimmed first
  vector<int> prio(nmb);
  for (int ki=0; ki<nmb; ++ki)
    prio[ki] = ki;

  for (int ki=0; ki<nmb; ++ki)
    {
      for (int kj=ki+1; kj<nmb; ++kj)
	if (adjacent[prio[kj]].size() < adjacent[prio[ki]].size())
	  std::swap(prio[ki], prio[kj]);
    }
      
  // Compute trimming curves
  tpTolerances tptol = model_->getTolerances();
  double epsge = tptol.gap;
  vector<pair<pair<int,shared_ptr<CurveOnSurface> >,
	      pair<int,shared_ptr<CurveOnSurface> > > > trim_curves;
  for (int ki=0; ki<nmb; ++ki)
    {
      // Check if the surface corresponds to another one that is
      // selected for trimming
      int ka;
      for (ka=0; ka<prio[ki]; ++ka)
	if (mod_sfs_[ka].get() == mod_sfs_[prio[ki]].get())
	  break;
      if (ka < prio[ki])
	continue;  // Treated as same underlying surface

      shared_ptr<ParamSurface> surf = mod_sfs_[prio[ki]];
 #ifdef DEBUG
  std::ofstream of2("trim_curves.g2");
  surf->writeStandardHeader(of2);
  surf->write(of2);
#endif
#ifdef DEBUG
  std::ofstream adj_sf("curr_adj_sfs.g2");
  surf->writeStandardHeader(adj_sf);
  surf->write(adj_sf);
#endif
      vector<int> ix2;
      vector<pair<shared_ptr<CurveOnSurface>, 
		  shared_ptr<CurveOnSurface> > > all_seg;
      vector<int> all_type;
      for (size_t kr=0; kr<adjacent[prio[ki]].size(); ++kr)
	{
	  vector<shared_ptr<CurveOnSurface> > int_segments1, int_segments2;
	  int adj_ix = adjacent[prio[ki]][kr];
	  shared_ptr<ParamSurface> surf2 = mod_sfs_[adj_ix];
#ifdef DEBUG
	  std::ofstream int_sf("int_sfs.g2");
	  surf->writeStandardHeader(int_sf);
	  surf->write(int_sf);
	  surf2->writeStandardHeader(int_sf);
	  surf2->write(int_sf);
	  surf2->writeStandardHeader(adj_sf);
	  surf2->write(adj_sf);
#endif
	  // Check for almost tangential intersections
	  size_t kh;
	  for (kh=0; kh<adj_ix_smooth.size(); ++kh)
	    {
	      if ((adj_ix_smooth[kh].first == prio[ki] &&
		   adj_ix_smooth[kh].second == adj_ix) ||
		  (adj_ix_smooth[kh].second == prio[ki] &&
		   adj_ix_smooth[kh].first == adj_ix))
		break;
	    }
	  if (kh < adj_ix_smooth.size())
	    {
	      getTangentialInt(surf, surf2, int_segments1, 
			       int_segments2, epsge);
	    }
	  else
	    {
	      try {
		BoundedUtils::getIntersectionCurve(surf, surf2, int_segments1, 
						   int_segments2, 0.5*epsge);
	      } catch (...) {
		THROW("Failed intersecting the two spline surfaces.");
	      }
	    }

	  for (size_t kh=0; kh<int_segments1.size(); ++kh)
	    {
	      ix2.push_back(adj_ix);
	      all_seg.push_back(make_pair(int_segments1[kh], 
					  int_segments2[kh]));
	      all_type.push_back(1);  // Currently computed intersection curve
#ifdef DEBUG
	      of2 << "100 1 0 4 255 0 0 255" << std::endl;
	      int_segments1[kh]->spaceCurve()->write(of2);
#endif
	    }
	}

      if (smoothedges[prio[ki]].size() > 0)
	{
	  // Add tangential intersection curves
	  for (size_t kr=0; kr<smoothedges[prio[ki]].size(); ++kr)
	    {
	      shared_ptr<ParamCurve> cv1 = 
		smoothedges[prio[ki]][kr].second->geomCurve();
	      if (!smoothedges[prio[ki]][kr].second->twin())
		continue;
	      shared_ptr<ParamCurve> cv2 = 
		smoothedges[prio[ki]][kr].second->twin()->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv1 =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv1);
	      shared_ptr<CurveOnSurface> sf_cv2 =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv2);
	      if (sf_cv1.get())
		cv1 = sf_cv1->spaceCurve();
	      if (sf_cv2.get())
		cv2 = sf_cv2->spaceCurve();

	      int i2 = smoothedges[prio[ki]][kr].first;

	      shared_ptr<CurveOnSurface> guide1(new CurveOnSurface(mod_sfs_[prio[ki]], cv1, false));
	      shared_ptr<CurveOnSurface> guide2(new CurveOnSurface(mod_sfs_[i2], cv2, false));

	      shared_ptr<CurveOnSurface> int_cv1, int_cv2;
	      CurveCreators::guidedIntersectionCurve(guide1, guide1->startparam(),
						     guide1->endparam(), guide2,
						     guide2->startparam(), guide2->endparam(),
						     epsge, int_cv1, int_cv2);
	      ix2.push_back(i2);
	      all_seg.push_back(make_pair(int_cv1, int_cv2));
	      all_type.push_back(2); // Guided intersection curve
#ifdef DEBUG
	      of2 << "100 1 0 4 100 100 55 255" << std::endl;
	      int_cv1->spaceCurve()->write(of2);
#endif
	    }
	}

      // Fetch previously computed trimming curves
      int nmb_curr_seg = (int)all_seg.size();
	  for (size_t kr=0; kr<trim_curves.size(); )
	{
	  if (trim_curves[kr].first.first == prio[ki] ||
	      trim_curves[kr].second.first == prio[ki])
	    {
	      if (trim_curves[kr].first.first == prio[ki])
		{
		  ix2.push_back(trim_curves[kr].second.first);
		  all_seg.push_back(make_pair(trim_curves[kr].first.second,
					      trim_curves[kr].second.second));
		  all_type.push_back(3); // Previously computed curve
		}
	      else
		{
		  ix2.push_back(trim_curves[kr].first.first);
		  all_seg.push_back(make_pair(trim_curves[kr].second.second,
					      trim_curves[kr].first.second));
		  all_type.push_back(3); 
		}
#ifdef DEBUG
	      of2 << "100 1 0 4 0 255 0 255" << std::endl;
	      all_seg[all_seg.size()-1].first->spaceCurve()->write(of2);
#endif
	      trim_curves.erase(trim_curves.begin()+kr);  // May be changed
	    }
	  else
	    ++kr;
	}

      // Split intersecting intersection curves
	  splitIntersectionCurves(ix2, all_seg, all_type, nmb_curr_seg, epsge,
				  tptol.neighbour);
 #ifdef DEBUG
    std::ofstream of3("int_seg2.g2");
    for (size_t kr=0; kr<all_seg.size(); ++kr)
      {
	of3 << "100 1 0 4 0 255 0 255" << std::endl;
	all_seg[kr].first->spaceCurve()->write(of3);
      }
#endif

      // Remove curves with a loose end
      cleanIntersectionCurves(surf, ix2, all_seg, epsge, 
			      tptol.neighbour, 5.0*tptol.neighbour);

 #ifdef DEBUG
    std::ofstream of4("int_seg3.g2");
    for (size_t kr=0; kr<all_seg.size(); ++kr)
      {
	of4 << "100 1 0 4 255 0 0 255" << std::endl;
	all_seg[kr].first->spaceCurve()->write(of4);
      }
#endif

      // // Remember trimming curves
      // for (size_t kr=0; kr<all_seg.size(); ++kr)
      // 	{
      // 	  trim_curves.push_back(make_pair(make_pair(prio[ki], all_seg[kr].first),
      // 					  make_pair(ix2[kr], all_seg[kr].second)));
      // 	}

      // // Remove the current surface as adjacent as the corresponding
      // // trimming curves are already computed and stored
      // for (int kj=0; kj<nmb; ++kj)
      // 	{
      // 	  if (kj == prio[ki])
      // 	    continue;

      // 	  for (size_t kr=0; kr<adjacent[kj].size(); )
      // 	    {
      // 	      if (adjacent[kj][kr] == prio[ki])
      // 		adjacent[kj].erase(adjacent[kj].begin()+kr);
      // 	      else
      // 		++kr;
      // 	    }
      // 	}
    // }

  // // Compute trimmed surfaces
  // for (int ki=0; ki<nmb; ++ki)
  //   {
  //     // Check if the surface is trimmed already
  //     int ka;
  //     for (ka=0; ka<prio[ki]; ++ka)
  // 	if (mod_sfs_[ka].get() == mod_sfs_[prio[ki]].get())
  // 	  break;
  //     if (ka < prio[ki])
  // 	continue;  // Treated as same underlying surface

       // shared_ptr<ParamSurface> surf = mod_sfs_[prio[ki]];

      // Collect trimming curves
    vector<shared_ptr<CurveOnSurface> > curr_trim(all_seg.size());
    for (size_t kr=0; kr<all_seg.size(); ++kr)
      curr_trim[kr] = all_seg[kr].first;
      // for (size_t kr=0; kr<trim_curves.size(); ++kr)
      // 	{
      // 	  if (trim_curves[kr].first.first == prio[ki])
      // 	    curr_trim.push_back(trim_curves[kr].first.second);
      // 	  else if (trim_curves[kr].second.first == prio[ki])
      // 	    curr_trim.push_back(trim_curves[kr].second.second);
      // 	}

      // Perform trimming
      vector<CurveLoop> loops = 
	SurfaceTools::absolutelyAllBoundarySfLoops(surf, epsge);
      shared_ptr<BoundedSurface> bd_sf = 
	shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

      vector<shared_ptr<BoundedSurface> > trimmed =
	BoundedUtils::splitWithTrimSegments(bd_sf, curr_trim, epsge);

      // Select trimmed surface entity
      // Fetch point to guide the trimming
      Point midpoint = surf->point(midpar[prio[ki]][0], midpar[prio[ki]][1]);
#ifdef DEBUG
      std::ofstream of5("curr_trim_sfs.g2");
      for (size_t kr=0; kr<trimmed.size(); ++kr)
	{
	  trimmed[kr]->writeStandardHeader(of5);
	  trimmed[kr]->write(of5);
}
      of5 << "400 1 0 4 255 0 0 255" << std::endl;
      of5 << "1" << std::endl;
      of5 << midpoint << std::endl;
#endif

      // Compute distance between an arbitrary point in the surface
      // and the trimming shell
      int ix = -1;
      if (trimmed.size() > 1)
	{
	  double u1, v1, d1=HUGE;
	  Point clo_pt1;
	  for (size_t kj=0; kj<trimmed.size(); ++kj)
	    {
	      double u2, v2, d2;
	      Point clo_pt2;
	      trimmed[kj]->closestPoint(midpoint, u2, v2, clo_pt2, d2, epsge);
	      if (d2 < d1)
		{
		  ix = (int)kj;
		  d1 = d2;
		}
	    }
	}
      if (trimmed.size() > 0 && ix >= 0)
	{
	  trim_sfs_[prio[ki]] = trimmed[ix];
#ifdef DEBUG
	  trimmed[ix]->writeStandardHeader(of);
	  trimmed[ix]->write(of);
#endif
	  // Remove curves that are not used for the selected surface
	  for (size_t kr=0; kr<all_seg.size(); )
	    {
	      Point endpts[2];
	      endpts[0] = 
		all_seg[kr].first->ParamCurve::point(all_seg[kr].first->startparam());
	      endpts[1] = 
		all_seg[kr].first->ParamCurve::point(all_seg[kr].first->endparam());
	      int ka;
	      for (ka=0; ka<2; ++ka)
		{
		  double u1, v1, d1;
		  Point clo_pt1;
		  trimmed[ix]->closestPoint(endpts[ka], u1, v1, clo_pt1, 
					    d1, epsge);
		  if (d1 > tptol.neighbour)
		    break;
		}
	      if (ka < 2)
		{
		  ix2.erase(ix2.begin()+kr);
		  all_seg.erase(all_seg.begin()+kr);
		}
	      else
		++kr;
	    }

	  // Remember trimming curves
	  for (size_t kr=0; kr<all_seg.size(); ++kr)
	    {
	      trim_curves.push_back(make_pair(make_pair(prio[ki], all_seg[kr].first),
					      make_pair(ix2[kr], all_seg[kr].second)));
	    }

	  // Remove the current surface as adjacent as the corresponding
	  // trimming curves are already computed and stored
	  for (int kj=0; kj<nmb; ++kj)
	    {
	      if (kj == prio[ki])
		continue;
	      
	      for (size_t kr=0; kr<adjacent[kj].size(); )
		{
		  if (adjacent[kj][kr] == prio[ki])
		    adjacent[kj].erase(adjacent[kj].begin()+kr);
		  else
		    ++kr;
		}
	    }
	}
    }

}

//==========================================================================
shared_ptr<ParamSurface> 
ModelUpdate::getUpdatedSurf(shared_ptr<ParamSurface> surf,
			    shared_ptr<ftSurface> face,
			    vector<Point>& par,
			    vector<Point>& diff,
			    double max_dist,
			    shared_ptr<ParamSurface>& base_sf)
//==========================================================================
{
  if (diff.size() == 0)
    return surf;  // No update required

  shared_ptr<SplineSurface> unified_sf;

  double eps = model_->getTolerances().neighbour;  // What is the appropriate value?
  double knot_diff_tol = 1.0e-10;

  shared_ptr<ParamSurface> curr_surf = (base_sf.get()) ? base_sf : surf;

  // Represent surface as tensor product splines
  shared_ptr<SplineSurface> tmp_spline;
  ElementarySurface *elem_sf = curr_surf->elementarySurface();
  
  // Fetch a possible underlying spline surface
  SplineSurface *base = NULL;
  if (elem_sf != NULL) 
    base = elem_sf->createNonRationalSpline(model_->getTolerances().gap);
  if (base == NULL)
    base = curr_surf->getSplineSurface();
  if (!base || base->rational())
    {
      // Approximate surface/underlying surface with non-rational surface
      vector<shared_ptr<ParamSurface> > base_sfs;
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(curr_surf);
      if (bd_sf.get())
	base_sfs.push_back(bd_sf->underlyingSurface());
      else
	base_sfs.push_back(surf);

      try {
	OffsetSurfaceStatus status = 
	  OffsetSurfaceUtils::offsetSurfaceSet(base_sfs, 0.0, 0.5*eps, 
					       tmp_spline);
      }
      catch (...)
	{
	  return unified_sf;
	}
      base = tmp_spline.get();
    }

  shared_ptr<SplineSurface> mod_base;
  if (!base_sf.get())
    {
      // Check that the surface size is appropriate compared to the
      // corresponding face
      mod_base = updateSurfaceSize(base, face);
      if (mod_base)
	base = mod_base.get();
    }

  if (max_dist < eps)
    {
      unified_sf = (mod_base.get()) ? mod_base :
	shared_ptr<SplineSurface>(base->clone());
      return unified_sf;
    }
      
  // Extract initial knot values
  vector<double> init_knots_u, init_knots_v;
  base->basis_u().knotsSimple(init_knots_u);
  base->basis_v().knotsSimple(init_knots_v);
  size_t nmb_init_u = init_knots_u.size();
  size_t nmb_init_v = init_knots_v.size();

  // Fetch information about the face domain
  RectDomain dom = face->surface()->containingDomain();
  vector<double> par_u(2), par_v(2);
  par_u[0] = dom.umin();
  par_u[1] = dom.umax();
  par_v[0] = dom.vmin();
  par_v[1] = dom.vmax();

  // Fetch size and shape information
  double size_u, size_v;
  base->estimateSfSize(size_u, size_v);
  DirectionCone tcone1 = base->tangentCone(true);
  DirectionCone tcone2 = base->tangentCone(false);
  int nmb_inner_u = 0;
  int nmb_inner_v = 0;

  // Add inner knots according to surface size
  double divfac1 = 1.1; //3.0;
  nmb_inner_u += (int)(size_u/(divfac1*sf_side_mean_));
  nmb_inner_v += (int)(size_v/(divfac1*sf_side_mean_));

  // Add inner knots according to shape
  double divfac2 = 1.0;
  int nmb_u = (tcone1.greaterThanPi()) ? 6 : (int)(tcone1.angle()/divfac2);
  int nmb_v = (tcone2.greaterThanPi()) ? 6 : (int)(tcone2.angle()/divfac2);
  nmb_inner_u = std::max(nmb_inner_u, nmb_u);
  nmb_inner_v = std::max(nmb_inner_v, nmb_v);

  double del_u = 
    (base->endparam_u() - base->startparam_u())/(double)(nmb_inner_u+1);
  double del_v = 
    (base->endparam_v() - base->startparam_v())/(double)(nmb_inner_v+1);

  // Modify face domain limits if they are close to an initial knot
  double start_u, end_u, start_v, end_v;
  double facu = 0.1, facv = 0.1;

  // Avoid extra knots if the base surface side is small
  if (size_u < 0.8*sf_side_mean_)
    facu = 0.75;
  if (size_v < 0.8*sf_side_mean_)
    facv = 0.75;
  
  if (par_u[0] < init_knots_u[0] + facu*del_u)
    par_u[0] = base->startparam_u();
  if (par_u[1] > init_knots_u[nmb_init_u-1] - facu*del_u)
    par_u[1] = base->endparam_u();
  if (par_v[0] < init_knots_v[0] + facv*del_v)
    par_v[0] = base->startparam_v();
  if (par_v[1] > init_knots_v[nmb_init_v-1] - facv*del_v)
    par_v[1] = base->endparam_v();

  if (par_u[0] >= init_knots_u[0] && 
      par_u[1] <= init_knots_u[nmb_init_u-1] &&
      par_v[0] >= init_knots_v[0] && 
      par_v[1] <= init_knots_v[nmb_init_v-1])
    {
      size_t kr;
      for (kr=1; kr<nmb_init_u; ++kr)
	if (init_knots_u[kr-1] <= par_u[0] && init_knots_u[kr] > par_u[0])
	  break;
      if (par_u[0] - init_knots_u[kr-1] < 0.1*del_u)
	par_u[0] = init_knots_u[kr-1];
      else if (init_knots_u[kr] - par_u[0] < 0.5*del_u)
	par_u[0] = std::max(init_knots_u[kr-1], init_knots_u[kr]-0.5*del_u);
      for (; kr<nmb_init_u; ++kr)
	if (init_knots_u[kr-1] < par_u[1] && init_knots_u[kr] >= par_u[1])
	  break;
      if (par_u[1] - init_knots_u[kr-1] < 0.5*del_u)
	par_u[1] = std::min(init_knots_u[kr-1]+0.5*del_u, init_knots_u[kr]);
      else if (init_knots_u[kr] - par_u[0] < 0.1*del_u)
	par_u[1] = init_knots_u[kr];

      for (kr=1; kr<nmb_init_v; ++kr)
	if (init_knots_v[kr-1] <= par_v[0] && init_knots_v[kr] > par_v[0])
	  break;
      if (par_v[0] - init_knots_v[kr-1] < 0.1*del_v)
	par_v[0] = init_knots_v[kr-1];
      else if (init_knots_v[kr] - par_v[0] < 0.5*del_v)
	par_v[0] = std::max(init_knots_v[kr-1], init_knots_v[kr]-0.5*del_v);
      for (; kr<nmb_init_v; ++kr)
	if (init_knots_v[kr-1] < par_v[1] && init_knots_v[kr] >= par_v[1])
	  break;
      if (par_v[1] - init_knots_v[kr-1] < 0.5*del_v)
	par_v[1] = std::min(init_knots_v[kr-1]+0.5*del_v, init_knots_v[kr]);
      else if (init_knots_v[kr] - par_v[0] < 0.1*del_v)
	par_v[1] = init_knots_v[kr];

      start_u = par_u[0];
      end_u = par_u[1];
      start_v = par_v[0];
      end_v = par_v[1];

      // Remove domain parameters identical with surface domain
      if (fabs(init_knots_u[nmb_init_u-1]-par_u[1]) < 0.1*del_u)
	{
	  par_u.pop_back();
	  end_u = init_knots_u[nmb_init_u-1];
	}
      if (fabs(par_u[0]-init_knots_u[0]) < 0.1*del_u)
	{
	  par_u.erase(par_u.begin(),par_u.begin()+1);
	  start_u = init_knots_u[0];
	}
      if (fabs(init_knots_v[nmb_init_v-1]-par_v[1]) < 0.1*del_v)
	{
	  par_v.pop_back();
	  end_v =init_knots_v[nmb_init_v-1];
	} 
      if (fabs(par_v[0]-init_knots_v[0]) < 0.1*del_v)
	{
	  par_v.erase(par_v.begin(),par_v.begin()+1);
	  start_v = init_knots_v[0];
	}
    }
  else
    {
      start_u = init_knots_u[0];
      end_u = init_knots_u[nmb_init_u-1];
      start_v = init_knots_v[0];
      end_v = init_knots_v[nmb_init_v-1];
      par_u.clear();
      par_v.clear();
    }

  // Create initial information for difference surface
  int dim = diff[0].dimension();
  int order = 4; 
  nmb_inner_u += (int)par_u.size();
  nmb_inner_v += (int)par_v.size();
  vector<double> uknots(2*order+nmb_inner_u);
  vector<double> vknots(2*order+nmb_inner_v);
  for (int kj=0; kj<order; ++kj)
    {
      uknots[kj] = base->startparam_u();
      uknots[order+nmb_inner_u+kj] = base->endparam_u();
      vknots[kj] = base->startparam_v();
      vknots[order+nmb_inner_v+kj] = base->endparam_v();
    }

  int u_ind = 0, v_ind = 0;
  if (start_u > uknots[0])
    {
      uknots[order] = start_u;
      u_ind++;
    }
  if (end_u < uknots[order+nmb_inner_u])
    uknots[order+nmb_inner_u-1] = end_u;
  if (start_v > vknots[0])
    {
      vknots[order] = start_v;
      v_ind++;
    }
  if (end_v < vknots[order+nmb_inner_v])
    vknots[order+nmb_inner_v-1] = end_v;

  // Remove knots at surface boundaries
  init_knots_u.pop_back();
  init_knots_u.erase(init_knots_u.begin(), init_knots_u.begin()+1);
  init_knots_v.pop_back();
  init_knots_v.erase(init_knots_v.begin(), init_knots_v.begin()+1);

  nmb_inner_u -= (int)par_u.size();
  nmb_inner_v -= (int)par_v.size();
  del_u = (end_u - start_u)/(double)(nmb_inner_u+1);
  del_v = (end_v - start_v)/(double)(nmb_inner_v+1);

  for (int kj=0; kj<nmb_inner_u; ++kj)
    {
      double knot = start_u + (kj+1)*del_u;
      double knot_dist = HUGE;
      double curr_knot;
      for (size_t kr=0; kr<init_knots_u.size(); ++kr)
	{
	  double curr_dist = fabs(init_knots_u[kr]-knot);
	  if (curr_dist < knot_dist)
	    {
	      knot_dist = curr_dist;
	      curr_knot = init_knots_u[kr];
	    }
	  else if (curr_dist > knot_dist)
	    break;  // Inner knots are sorted. No closer knot exist
	}
      if (knot_dist < 0.5*del_u)
	uknots[order+u_ind+kj] = curr_knot;
      else
	uknots[order+u_ind+kj] = knot;
    }

  for (int kj=0; kj<nmb_inner_v; ++kj)
    {
      double knot = start_v + (kj+1)*del_v;
      double knot_dist = HUGE;
      double curr_knot;
      for (size_t kr=0; kr<init_knots_v.size(); ++kr)
	{
	  double curr_dist = fabs(init_knots_v[kr]-knot);
	  if (curr_dist < knot_dist)
	    {
	      knot_dist = curr_dist;
	      curr_knot = init_knots_v[kr];
	    }
	  else if (curr_dist > knot_dist)
	    break;  // Inner knots are sorted. No closer knot exist
	}
      if (knot_dist < 0.5*del_v)
	vknots[order+v_ind+kj] = curr_knot;
      else
	vknots[order+v_ind+kj] = knot;
    }


  vector<shared_ptr<SplineSurface> > diff_sfs(dim);
  for (int kj=0; kj<dim; ++kj)
    {
      // Fetch data points
      vector<double> data(3*diff.size());
      for (size_t ki=0; ki<diff.size(); ++ki)
	{
	  data[3*ki] = par[ki][0];
	  data[3*ki+1] = par[ki][1];
	  data[3*ki+2] = diff[ki][kj];
	}

      // Create approximation machinery
      LRSurfApprox approx(order, uknots, order, vknots, data, 
			  1, eps, true, 0.0, false, false);
      approx.setUseMBA(true);

      // if (!remove_)
      // 	approx.addLowerConstraint(0.0);

      // Add knot line information from base surface
      if (init_knots_u.size() > 0 || init_knots_v.size() > 0)
	approx.setInitialKnots(init_knots_u, init_knots_v);

#ifdef DEBUG
      approx.setVerbose(true);
#endif

      double maxdist, avdist, avdist_total; // will be set below
      int nmb_out_eps;        // will be set below
      shared_ptr<LRSplineSurface> lr_surf = 
	approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, 
			     max_iter_);

#ifdef DEBUG
      std::cout << "No. elements: " << lr_surf->numElements();
      std::cout << ", maxdist= " << maxdist << "avdist= " << avdist_total;
      std::cout << ", avdist(out)= " << avdist;
      std::cout << ", nmb out= " << nmb_out_eps << std::endl;
#endif

      diff_sfs[kj] = shared_ptr<SplineSurface>(lr_surf->asSplineSurface());
      int stop_break = 1;
    }

  // Ensure that the difference surface have the same knot vectors
  GeometryTools::unifySurfaceSplineSpace(diff_sfs, knot_diff_tol);
  
  // Create 3D difference surface
  int ncoef = diff_sfs[0]->numCoefs_u()*diff_sfs[0]->numCoefs_v();
  vector<double> coefs(dim*ncoef);
  for (int kj=0; kj<dim; ++kj)
    {
      vector<double>::const_iterator it = diff_sfs[kj]->coefs_begin();
      size_t ki;
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

  double curr_max = 0.0;
  double curr_av = 0.0;
  for (size_t kr=0; kr<diff.size(); ++kr)
    {
      Point diff_pos = diff3D->ParamSurface::point(par[kr][0],par[kr][1]);
      double dist = diff[kr].dist(diff_pos);
      maxdist_ = std::max(maxdist_, dist);
      avdist_ += dist;
      curr_max = std::max(curr_max, dist);
      curr_av += dist;
    }
  curr_av /= (int)diff.size();

  // Create 3D unified surface
  unified_sf = 
    GeometryTools::surfaceSum(*base, 1.0, *diff3D, 1.0, knot_diff_tol);

#ifdef DEBUG
  std::ofstream of("curr_mod_sf.g2");
  unified_sf->writeStandardHeader(of);
  unified_sf->write(of);
#endif

  if (base_sf.get())
    base_sf = unified_sf;

  return unified_sf;
}

//==========================================================================
void ModelUpdate::splitCommonBase(vector<int>& same_face_ix,
				  shared_ptr<ParamSurface> base_sf,
				  pair<bool, bool>& base_closed)
//==========================================================================
{
  double eps = 1.0e-6;

  RectDomain dom_base = base_sf->containingDomain();
  SplineSurface* base_spline = base_sf->getSplineSurface();
  shared_ptr<SplineSurface> tmp_spline;
  if (base_spline == NULL)
    {
      tmp_spline = shared_ptr<SplineSurface>(base_sf->asSplineSurface());
      base_spline = tmp_spline.get();
    }

  for (size_t ki=0; ki<same_face_ix.size(); ++ki)
    {
      shared_ptr<ftSurface> face1 = model_->getFace(same_face_ix[ki]);
      RectDomain dom1 = face1->surface()->containingDomain();
      for (size_t kj=ki+1; kj<same_face_ix.size(); ++kj)
	{
	  shared_ptr<ftSurface> face2 = model_->getFace(same_face_ix[kj]);
	  RectDomain dom2 = face2->surface()->containingDomain();

	  if (!dom1.overlap(dom2, eps))
	    continue;

	  // Check minimum domain overlap
	  if (std::min(fabs(dom1.umin()-dom2.umax()), 
		       fabs(dom2.umin()-dom1.umax())) <
	      std::min(fabs(dom1.vmin()-dom2.vmax()), 
		       fabs(dom2.vmin()-dom1.vmax())))
	    {
	      // Split in the first parameter direction
	      double upar;
	      if (fabs(dom1.umin()-dom2.umax()) < fabs(dom2.umin()-dom1.umax()))
		  upar = 0.5*(dom1.umin()+dom2.umax());
	      else
		  upar = 0.5*(dom2.umin()+dom1.umax());

	      shared_ptr<SplineSurface> sub1(base_spline->subSurface(dom_base.umin(),
								     dom_base.vmin(),
								     upar,
								     dom_base.vmax()));
	      shared_ptr<SplineSurface> sub2(base_spline->subSurface(upar,
								     dom_base.vmin(),
								     dom_base.umax(),
								     dom_base.vmax()));
	      if (base_closed.first)
		{
		  // Ensure continuity across seam
		  Point dummy;
		  GeometryTools::averageBoundaryCoefs(sub2, 1, false, 
						      sub1, 0, false, 
						      false, dummy, false, 
						      dummy, false);
		}

	      // Replace modified surface
	      mod_sfs_[same_face_ix[ki]] = (dom1.umin() < dom2.umin()) ?
		sub1 : sub2;
	      mod_sfs_[same_face_ix[kj]] = (dom1.umin() < dom2.umin()) ?
		sub2 : sub1;
	    }
	  else
	    {
	      // Split in the second parameter direction
	      double vpar;
	      if (fabs(dom1.vmin()-dom2.vmax()) < fabs(dom2.vmin()-dom1.vmax()))
		  vpar = 0.5*(dom1.vmin()+dom2.vmax());
	      else
		  vpar = 0.5*(dom2.vmin()+dom1.vmax());
	      shared_ptr<SplineSurface> sub1(base_spline->subSurface(dom_base.umin(),
								     dom_base.vmin(),
								     dom_base.umax(),
								     vpar));
	      shared_ptr<SplineSurface> sub2(base_spline->subSurface(dom_base.umin(),
								     vpar,
								     dom_base.umax(),
								     dom_base.vmax()));
	      if (base_closed.second)
		{
		  // Ensure continuity across seam
		  Point dummy;
		  GeometryTools::averageBoundaryCoefs(sub2, 3, false, 
						      sub1, 2, false, 
						      false, dummy, false, 
						      dummy, false);
		}
	      // Replace modified surface
	      mod_sfs_[same_face_ix[ki]] = (dom1.vmin() < dom2.vmin()) ?
		sub1 : sub2;
	      mod_sfs_[same_face_ix[kj]] = (dom1.vmin() < dom2.vmin()) ?
		sub2 : sub1;
	    }
	}
    }
}

//==========================================================================
double ModelUpdate::meanSurfaceSide()
//==========================================================================
{
  double mean_size = 0.0;
  int nmb_sides = 0;
  int nmb = model_->nmbEntities();
  sf_size_.resize(nmb);
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> surf = model_->getSurface(ki);
      double len_u, len_v;
      shared_ptr<BoundedSurface> bd_surf = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
      if (bd_surf.get())
	{
	  RectDomain dom = surf->containingDomain();
	  double umin = dom.umin();
	  double umax = dom.umax();
	  double vmin = dom.vmin();
	  double vmax = dom.vmax();
	  bd_surf->underlyingSurface()->estimateSubSfSize(umin, umax, len_u, 
							  vmin, vmax, len_v);
	}
      else
	surf->estimateSfSize(len_u, len_v);
      sf_size_[ki] = make_pair(len_u, len_v);
      mean_size += (len_u + len_v);
      nmb_sides += 2;
    }

  if (nmb_sides > 0)
    mean_size /= (double)nmb_sides;

  return mean_size;
}

//==========================================================================
shared_ptr<SplineSurface> 
ModelUpdate::updateSurfaceSize(SplineSurface *surf,
			       shared_ptr<ftSurface> face)
//==========================================================================
{
  double epsp = 1.0e-6;
  double tol = 10.0*model_->getTolerances().neighbour;
  double frac = 0.75;
  double frac2 = 0.2;

  // Check parameter domains
  RectDomain dom1 = surf->containingDomain();
  RectDomain dom2 = face->surface()->containingDomain();

  // First check if the base surface is unnecessarily large
  bool limit = false;
  double umin = dom1.umin();
  double umax = dom1.umax();
  double vmin = dom1.vmin();
  double vmax = dom1.vmax();
  double umin2 = dom2.umin();
  double umax2 = dom2.umax();
  double vmin2 = dom2.vmin();
  double vmax2 = dom2.vmax();
  double udel = umax2 - umin2;
  double vdel = vmax2 - vmin2;
  if (umin2 - umin > frac*udel)
    {
      limit = true;
      umin = umin2 - frac*udel;
    }
  if (umax - umax2 > frac*udel)
    {
      limit = true;
      umax = umax2 + frac*udel;
    }

  if (vmin2 - vmin > frac*vdel)
    {
      limit = true;
      vmin = vmin2 - frac*vdel;
    }
  if (vmax - vmax2 > frac*vdel)
    {
      limit = true;
      vmax = vmax2 + frac*vdel;
    }

  shared_ptr<SplineSurface> mod_surf;
  if (limit)
    mod_surf = shared_ptr<SplineSurface>(surf->subSurface(umin, vmin,
							  umax, vmax));

  // Check if the base surface should be extended
  SplineSurface *tmp = (mod_surf.get()) ? mod_surf.get() : surf;
  double del1=0.0, del2=0.0, del3=0.0, del4=0.0;
  double other;
  if (umin2 - umin > epsp)
    tmp->estimateSubSfSize(umin, umin2, del1,
			   vmin2, vmax2, other);
  if (umax - umax2 > epsp)
    tmp->estimateSubSfSize(umax2, umax, del2,
			   vmin2, vmax2, other);
  if (vmin2 - vmin > epsp)
    tmp->estimateSubSfSize(umin2, umax2, other,
			   vmin, vmin2, del3);
  if (vmax - vmax2 > epsp)
    tmp->estimateSubSfSize(umin2, umax2, other,
			   vmax2, vmax, del4);

  if (del1 < tol || del2 < tol || del3 < tol || del4 < tol)
    {
      if (!mod_surf.get())
	mod_surf = shared_ptr<SplineSurface>(surf->clone());
    }

  if (del1 < tol)
    mod_surf->enlarge(frac2*udel, true, false);
  if (del2 < tol)
    mod_surf->enlarge(frac2*udel, true, true);
  if (del3 < tol)
    mod_surf->enlarge(frac2*vdel, false, false);
  if (del4 < tol)
    mod_surf->enlarge(frac2*vdel, false, true);
   
  return mod_surf;
}

//==========================================================================
void
ModelUpdate::identifyCommonBaseSfs(vector<vector<shared_ptr<ftSurface> > >& same_faces,
				   vector<vector<int> >& same_faces_ix,
				   vector<shared_ptr<ParamSurface> >& base_sfs,
				   vector<pair<bool, bool> >& base_closed)
//==========================================================================
{
  // Identify smooth edges in the original model
  vector<shared_ptr<ftSurface> > faces = model_->allFaces();
  FaceConnectivityUtils<ftEdgeBase,ftSurface> connectivity;
  vector<ftEdgeBase*> vec;
  connectivity.smoothEdges(faces, vec, model_->getTolerances().bend);

#ifdef DEBUG
  if (vec.size() > 0)
    {
      std::ofstream edgof("smoothedges.g2");
      for (size_t kr=0; kr<vec.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = vec[kr]->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sf_cv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  if (sf_cv.get())
	    cv = sf_cv->spaceCurve();
	  cv->writeStandardHeader(edgof);
	  cv->write(edgof);
	}
    }
#endif
  // Remove edges with a small, but still significant angle
  vector<ftEdgeBase*> vec2;
  for (size_t kr=0; kr<vec.size(); )
    {
      if (vec[kr]->hasConnectivityInfo())
	{
	  int status = vec[kr]->getConnectivityInfo()->BestStatus();
	  if (status > 0)
	    {
	      vec2.push_back(vec[kr]);
	      vec.erase(vec.begin()+kr);
	    }
	  else
	    ++kr;
	}
      else 
	++kr;
    }
#ifdef DEBUG
  if (vec.size() > 0)
    {
      std::ofstream edgof2("smoothedges2.g2");
      for (size_t kr=0; kr<vec.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = vec[kr]->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sf_cv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  if (sf_cv.get())
	    cv = sf_cv->spaceCurve();
	  cv->writeStandardHeader(edgof2);
	  cv->write(edgof2);
	}
    }
#endif
  smooth_edgs_ = vec;
  smooth_edgs2_ = vec2;

  // Identify surfaces having the same underlying base surface
  vector<shared_ptr<ftSurface> > all_faces = model_->allFaces();
  tpTolerances tptol = model_->getTolerances();
  SurfaceModelUtils::sameUnderlyingSurf(all_faces, tptol.gap, tptol.kink,
					same_faces, base_sfs);

  // Dismiss instances where the corresponding faces are non-adjacent
  for (size_t kr=0; kr<same_faces.size(); )
    {
      size_t kj, kh;
      for (kj=0; kj<same_faces[kr].size(); ++kj)
	{
	  for (kh=0; kh<same_faces[kr].size(); ++kh)
	    {
	      if (kj == kh)
		continue;

	      bool smooth;
	      if (same_faces[kr][kj]->isAdjacent(same_faces[kr][kh].get(), 
						 smooth))
		break;
	    }
	  if (kh == same_faces[kr].size())
	    break;  // No adjacent face is found
	}

      if (kj < same_faces[kr].size())
	{
	  same_faces.erase(same_faces.begin()+kr);
	  base_sfs.erase(base_sfs.begin()+kr);
	}
      else
	++kr;
    }

  // Recompute parameter information for points corresponding to
  // the found faces
  double eps = model_->getTolerances().gap;  // Maybe neighbour is better?
  for (size_t kj=0; kj<same_faces.size(); ++kj)
    {
      shared_ptr<ElementarySurface> elem_sf = 
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(base_sfs[kj]);
      if (!elem_sf.get())
	continue;  // Really the same underlying surface

      shared_ptr<ParamSurface> spline(elem_sf->createNonRationalSpline(eps));
      base_sfs[kj] = spline;  // To ensure same parameterization

      for (size_t kr=0; kr<same_faces[kj].size(); ++kr)
	{
#ifdef DEBUG
	  vector<Point> pos1, pos2;
#endif
	  int face_ix = model_->getIndex(same_faces[kj][kr]);
	  for (size_t ki=0; ki<points_.size(); ++ki)
	    {
	      if (points_[ki].surfIx() == face_ix)
		{
		  Point pos = points_[ki].getPos();
		  double upar, vpar, dist;
		  Point clo_pos;
		  base_sfs[kj]->closestPoint(pos, upar, vpar, clo_pos, 
					     dist, eps);
		  points_[ki].setParameterInfo(face_ix, upar, vpar);
#ifdef DEBUG
		  pos1.push_back(pos);
		  pos2.push_back(base_sfs[kj]->point(upar,vpar));
#endif
		}
	    }
#ifdef DEBUG
	  std::ofstream ofpts("reparpts.g2");
	  ofpts << "400 1 0 4 255 0 0 255 " << std::endl;
	  ofpts << pos1.size() << std::endl;
	  for (size_t ki=0; ki<pos1.size(); ++ki)
	    ofpts << pos1[ki] << std::endl;
	  ofpts << "400 1 0 4 0 255 0 255 " << std::endl;
	  ofpts << pos2.size() << std::endl;
	  for (size_t ki=0; ki<pos2.size(); ++ki)
	    ofpts << pos2[ki] << std::endl;
#endif
	  int stop_break = 1;
	}
    }

  // Remove information about smooth edges placed between faces with
  // common underlying surface
  // First collect common edges
   vector<shared_ptr<ftEdge> > common_edgs;
  for (size_t kr=0; kr<same_faces.size(); ++kr)
    {
      size_t kj, kh;
      for (kj=0; kj<same_faces[kr].size()-1; ++kj)
	{
	  for (kh=kj+1; kh<same_faces[kr].size(); ++kh)
	    {
	      vector<shared_ptr<ftEdge> > curr_common_edgs = 
		same_faces[kr][kj]->getCommonEdges(same_faces[kr][kh].get());
	      if (curr_common_edgs.size() > 0)
		common_edgs.insert(common_edgs.end(), curr_common_edgs.begin(),
				   curr_common_edgs.end());
	    }
	}
    }
  for (size_t kr=0; kr<smooth_edgs_.size(); )
    {
      size_t kh;
      for (kh=0; kh<common_edgs.size(); ++kh)
	{
	  if (common_edgs[kh].get() == smooth_edgs_[kr] ||
	      common_edgs[kh]->twin() == smooth_edgs_[kr] ||
	      common_edgs[kh].get() == smooth_edgs_[kr]->twin() ||
	      common_edgs[kh]->twin() == smooth_edgs_[kr]->twin())
	  break;
	}
      if (kh < common_edgs.size())
	{
	  smooth_edgs_.erase(smooth_edgs_.begin()+kr);
	}
      else
	++kr;
    }
#ifdef DEBUG
  if (smooth_edgs_.size() > 0)
    {
      std::ofstream edgof3("smoothedges_final.g2");
      for (size_t kr=0; kr<smooth_edgs_.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = smooth_edgs_[kr]->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sf_cv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  if (sf_cv.get())
	    cv = sf_cv->spaceCurve();
	  cv->writeStandardHeader(edgof3);
	  cv->write(edgof3);
	}
    }
#endif

  // Merge remaining underlying surfaces meeting in a smooth edge to 
  // create common base surfaces
  // First identify groups of faces
  vector<vector<shared_ptr<ftSurface> > > face_grp;
  for (size_t kr=0; kr<smooth_edgs_.size(); ++kr)
    {
      vector<shared_ptr<ftSurface> > curr_grp;
      ftFaceBase* curr;
      if (smooth_edgs_[kr]->twin())
      curr = smooth_edgs_[kr]->geomEdge()->face();
      if (curr)
	curr_grp.push_back(model_->fetchAsSharedPtr(curr));
      if (smooth_edgs_[kr]->twin())
	{
	  curr = smooth_edgs_[kr]->twin()->geomEdge()->face();
	  if (curr)
	    curr_grp.push_back(model_->fetchAsSharedPtr(curr));
	}

      // Check if any of the faces are collected already
      size_t kh, j1, j2;
      for (kh=0; kh<face_grp.size(); ++kh)
	{
	  for (j1=0; j1<curr_grp.size(); ++j1)
	    {
	      for (j2=0; j2<face_grp[kh].size(); ++j2)
		{
		  if (curr_grp[j1].get() == face_grp[kh][j2].get())
		    break;
		}
	      if (j2 < face_grp[kh].size())
		break;
	    }
	  if (j1 < curr_grp.size())
	    break;
	}
      if (kh < face_grp.size())
	{
	  for (j1=0; j1<curr_grp.size(); ++j1)
	    {
	      for (j2=0; j2<face_grp[kh].size(); ++j2)
		{
		  if (curr_grp[j1].get() == face_grp[kh][j2].get())
		    break;
		}
	      if (j2 == face_grp[kh].size())
		face_grp[kh].push_back(curr_grp[j1]);
	    }
	}
      else
	face_grp.push_back(curr_grp);
    }

  size_t nmb_bases = same_faces.size();
  for (size_t kr=0; kr<face_grp.size(); ++kr)
    {
      shared_ptr<ParamSurface> curr_base = constructCommonBase(face_grp[kr]);
      if (curr_base.get())
	{
	  same_faces.push_back(face_grp[kr]);
	  base_sfs.push_back(curr_base);
	}
    }

  // Remove information about smooth edges placed between new faces with
  // common underlying surface
  // First collect common edges
  common_edgs.clear();
  for (size_t kr=nmb_bases; kr<same_faces.size(); ++kr)
    {
      size_t kj, kh;
      for (kj=0; kj<same_faces[kr].size()-1; ++kj)
	{
	  for (kh=kj+1; kh<same_faces[kr].size(); ++kh)
	    {
	      vector<shared_ptr<ftEdge> > curr_common_edgs = 
		same_faces[kr][kj]->getCommonEdges(same_faces[kr][kh].get());
	      if (curr_common_edgs.size() > 0)
		common_edgs.insert(common_edgs.end(), curr_common_edgs.begin(),
				   curr_common_edgs.end());
	    }
	}
    }
  for (size_t kr=0; kr<smooth_edgs_.size(); )
    {
      size_t kh;
      for (kh=0; kh<common_edgs.size(); ++kh)
	{
	  if (common_edgs[kh].get() == smooth_edgs_[kr] ||
	      common_edgs[kh]->twin() == smooth_edgs_[kr] ||
	      common_edgs[kh].get() == smooth_edgs_[kr]->twin() ||
	      common_edgs[kh]->twin() == smooth_edgs_[kr]->twin())
	  break;
	}
      if (kh < common_edgs.size())
	{
	  smooth_edgs_.erase(smooth_edgs_.begin()+kr);
	}
      else
	++kr;
    }

  // Collect indices of faces corresponding to the same base surface
  same_faces_ix.resize(same_faces.size());
  for (size_t kr=0; kr<same_faces.size(); ++kr)
    {
      same_faces_ix[kr].resize(same_faces[kr].size());
      for (size_t kh=0; kh<same_faces[kr].size(); ++kh)
	same_faces_ix[kr][kh] = model_->getIndex(same_faces[kr][kh]);
    }

  // Collect information about closed base surfaces
  base_closed.resize(base_sfs.size());
  double neighbour = model_->getTolerances().neighbour;
  for (size_t kr=0; kr<base_sfs.size(); ++kr)
    {
      bool close_u, close_v;
      SurfaceTools::checkSurfaceClosed(*base_sfs[kr], close_u, close_v,
				       neighbour);
      base_closed[kr] = make_pair(close_u, close_v);
    }
      
  // Extend surfaces in non-closed direction to ensure intersections
  // between updated neighbours
  BoundingBox box = model_->boundingBox();
  double len = box.low().dist(box.high());
  len *= 0.1;
  for (size_t kr=0; kr<base_sfs.size(); ++kr)
    {
      shared_ptr<ElementarySurface> elem = 
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(base_sfs[kr]);
      shared_ptr<SplineSurface> spline = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(base_sfs[kr]);
      if (base_closed[kr].first && !base_closed[kr].second)
	{
	  if (elem.get())
	    elem->enlarge(0.0, 0.0, len, len);
	  else if (spline.get())
	    spline->enlarge(0.0, 0.0, len, len);
	}
      else if (!base_closed[kr].first && base_closed[kr].second)
	{
	  if (elem.get())
	    elem->enlarge(len, len, 0.0, 0.0);
	  else if (spline.get())
	    spline->enlarge(len, len, 0.0, 0.0);
	}
      else if (!(base_closed[kr].first || base_closed[kr].second))
	{
	  if (elem.get())
	    elem->enlarge(len, len, len, len);
	  else if (spline.get())
	    spline->enlarge(len, len, len, len);
	}
    }

#ifdef DEBUG
  std::ofstream ofbase("same_base.g2");
  for (size_t kr=0; kr<base_sfs.size(); ++kr)
    {
      base_sfs[kr]->writeStandardHeader(ofbase);
      base_sfs[kr]->write(ofbase);
    }
#endif
}

//==========================================================================
shared_ptr<ParamSurface>
ModelUpdate::constructCommonBase(vector<shared_ptr<ftSurface> >& same_faces)
//==========================================================================
{
  shared_ptr<ParamSurface> base;
  if (same_faces.size() == 0)
    return base;

  // Check consictency
  for (size_t kh=0; kh<same_faces.size(); ++kh)
    {
      ElementarySurface* elem1 = same_faces[kh]->surface()->elementarySurface();
      for (size_t kr=kh+1; kr<same_faces.size(); ++kr)
	{
	  ElementarySurface* elem2 = 
	    same_faces[kr]->surface()->elementarySurface();
	  if (elem1 && elem2 && elem1->instanceType() != elem2->instanceType())
	    return base;
	}
    }

  // Fetch face indices
  vector<int> face_ix(same_faces.size());
  for (size_t kh=0; kh<same_faces.size(); ++kh)
    face_ix[kh] = model_->getIndex(same_faces[kh]);

  // Sort the faces according to size
  for (size_t kh=0; kh<same_faces.size(); ++kh)
    {
      double sfsize1 = 
	sf_size_[face_ix[kh]].first*sf_size_[face_ix[kh]].second;
      for (size_t kr=kh+1; kr<same_faces.size(); ++kr)
	{
	  double sfsize2 = 
	    sf_size_[face_ix[kr]].first*sf_size_[face_ix[kr]].second;
	  if (sfsize2 > sfsize1)
	    {
	      std::swap(same_faces[kh], same_faces[kr]);
	      std::swap(face_ix[kh], face_ix[kr]);
	    }
	}
    }

  // Check if the smaller faces can be represented with the (possibly enlarged)
  // base surface of the largest face.
  // First find largest surface size of the remaining surfaces
  double len = 0.0;
  for (size_t kr=1; kr<same_faces.size(); ++kr)
    len = std::max(len, std::max(sf_size_[face_ix[kr]].first,
				 sf_size_[face_ix[kr]].second));

  shared_ptr<ParamSurface> surf = same_faces[0]->surface();
  ElementarySurface *elem = surf->elementarySurface();
  if (elem)
    {
      // Enlarge surface.       
      shared_ptr<ElementarySurface> elem2(elem->clone());
      elem2->enlarge(len, len, len, len);
      base = elem2;
    }
  else
    {
      shared_ptr<SplineSurface> spline(surf->asSplineSurface()); // Copy
      spline->enlarge(len, len, len, len);
      base = spline;
    }

  size_t kr;
  double eps = model_->getTolerances().gap;  // Maybe neighbour is better?
  for (kr=1; kr<same_faces.size(); ++kr)
    {
      int face_ix = model_->getIndex(same_faces[kr]);
      double eps2 = std::max(eps, std::max(sf_info_[face_ix].max_dist_,
					   fabs(sf_info_[face_ix].min_dist_)));
      // Create sampling points
      shared_ptr<ftPointSet> points(new ftPointSet());
      shared_ptr<ParamSurface> curr_sf = same_faces[kr]->surface();
      RectDomain dom = curr_sf->containingDomain();
      int nmb_sample = 10;
      vector<int> corner;
      AdaptSurface::getBoundaryData(curr_sf, dom, nmb_sample, points, corner);
      AdaptSurface::getInnerData(surf, dom, nmb_sample, points, false);

      int nmb = points->size();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  Vector3D curr = (*points)[ki]->getPoint();
	  Point pos(curr[0], curr[1], curr[2]);

	  // Project onto the proposed base surface
	  double upar, vpar, dist;
	  Point clo_pt;
	  base->closestPoint(pos, upar, vpar, clo_pt, dist, eps);
	  if (dist > eps2)
	    break;
	}
      if (ki < nmb)
	break;
    }
  if (kr < same_faces.size())
    {
      // Distance too large. Try another approach
      tpTolerances toptol = model_->getTolerances();
      shared_ptr<SurfaceModel> tmp_model(new SurfaceModel(toptol.gap, toptol.gap,
							  toptol.neighbour,
							  toptol.kink, toptol.bend,
							  same_faces));
      double dist;
      shared_ptr<ParamSurface> merged;
      try {
	merged = tmp_model->representAsOneSurface(dist, 3);
      }
      catch (...)
	{
	  base.reset();
	}

      if (dist > eps)
	base.reset();

      if (merged.get())
	{
	  shared_ptr<SplineSurface> spline_base(merged->asSplineSurface());
	  spline_base->enlarge(len, len, len, len);
	  base = spline_base;
	}
    }

  if (base.get())
    {
      // Reparameterize data points corresponding to the smaller faces
      int face_ix = model_->getIndex(same_faces[0]);
      for (size_t kr=1; kr<same_faces.size(); ++kr)
	{
	  int face_ix2 = model_->getIndex(same_faces[kr]);
	  for (size_t ki=0; ki<points_.size(); ++ki)
	    {
	      if (points_[ki].surfIx() == face_ix2)
		{
		  Point pos = points_[ki].getPos();
		  double upar, vpar, dist;
		  Point clo_pos;
		  base->closestPoint(pos, upar, vpar, clo_pos, dist, eps);
		  points_[ki].setParameterInfo(face_ix, upar, vpar);

		  // Update distance based on differerence between the
		  // previous and the new surface in the point.
		  // First check orientation
		  Point norm;
		  base->normal(norm, upar, vpar);
		  Point vec = pos - clo_pos;
		  int sgn = (norm*vec >= 0.0) ? 1 : -1;
		  double dist2 = points_[ki].getDist();
		  dist2 += (sgn*dist);
		  points_[ki].setDist(dist2);
		}
	    }
	}
    }

  return base;
}

//==========================================================================
 void 
 ModelUpdate::fetchAdjacent(int face_ix, vector<int>& adj_ix,
			    vector<pair<int, shared_ptr<ftEdge> > >& smooth,
			    vector<pair<int, int> >& smooth_ix)
//==========================================================================
{
  // Fetch associated face
  shared_ptr<ftSurface> face = model_->getFace(face_ix);
  shared_ptr<ParamSurface> surf = mod_sfs_[face_ix];
#ifdef DEBUG
  std::ofstream of("adj_info.g2");
  surf->writeStandardHeader(of);
  surf->write(of);
#endif

  vector<int> adj_ix_smooth;
  bool has_missing_twin = false;
  vector<double> missing_len;
  vector<Point> edg_point;
  vector<shared_ptr<ftEdge> > edgs = face->getAllEdges();
  for (size_t ki=0; ki<edgs.size(); ++ki)
    {
      if (edgs[ki]->twin() == NULL)
	{
	  has_missing_twin = true;
	  missing_len.push_back(edgs[ki]->estimatedCurveLength());
	  edg_point.push_back(edgs[ki]->point(0.5*(edgs[ki]->tMin()+edgs[ki]->tMax())));
	}
    }

  // Fetch neighbours
  vector<ftSurface*> neighbours;
  face->getAdjacentFaces(neighbours);

  // Check for more neighbours
  int nmb = model_->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      if ((int)ki == face_ix)
	continue;
      if (mod_sfs_[ki].get() == surf.get())
	{
	  shared_ptr<ftSurface> face2 = model_->getFace(ki);
	  vector<ftSurface*> neighbours2;
	  face2->getAdjacentFaces(neighbours2);
	  for (size_t kj=0; kj<neighbours2.size(); ++kj)
	    {
	      size_t kr;
	      for (kr=0; kr<neighbours.size(); ++kr)
		if (neighbours[kr] == neighbours2[kj])
		  break;
	      if (kr == neighbours.size())
		neighbours.push_back(neighbours2[kj]);
	    }

	  // Add also missing twin information
	  vector<shared_ptr<ftEdge> > edgs2 = face2->getAllEdges();
	  for (size_t kj=0; kj<edgs2.size(); ++kj)
	    {
	      if (edgs2[kj]->twin() == NULL)
		{
		  has_missing_twin = true;
		  missing_len.push_back(edgs2[kj]->estimatedCurveLength());
		  edg_point.push_back(edgs2[kj]->point(0.5*(edgs2[kj]->tMin()+
							    edgs2[kj]->tMax())));
		}
	    }
	}
    }

  // Fetch associated face indices
  for (size_t kj=0; kj<neighbours.size(); ++kj)
    {
      int ix = model_->getIndex(neighbours[kj]);
      if (mod_sfs_[ix].get() == surf.get())
	continue;

      // Check for a smooth connection
      bool found_smooth = false;
      vector<shared_ptr<ftEdge> > common_edgs = 
	face->getCommonEdges(neighbours[kj]);
      size_t kr, kh;
      for (kh=0; kh<common_edgs.size(); ++kh)
	{
	  for (kr=0; kr<smooth_edgs_.size(); ++kr)
	    {
	      if (common_edgs[kh].get() == smooth_edgs_[kr] ||
		  common_edgs[kh]->twin() == smooth_edgs_[kr] ||
		  common_edgs[kh].get() == smooth_edgs_[kr]->twin() ||
		  common_edgs[kh]->twin() == smooth_edgs_[kr]->twin())
		break;
	    }
	  if (kr < smooth_edgs_.size())
	    {
	      smooth.push_back(make_pair(ix, common_edgs[kh]));
	      found_smooth = true;
	    }

	  for (kr=0; kr<smooth_edgs2_.size(); ++kr)
	    {
	      if (common_edgs[kh].get() == smooth_edgs2_[kr] ||
		  common_edgs[kh]->twin() == smooth_edgs2_[kr] ||
		  common_edgs[kh].get() == smooth_edgs2_[kr]->twin() ||
		  common_edgs[kh]->twin() == smooth_edgs2_[kr]->twin())
		break;
	    }
	  if (kr < smooth_edgs2_.size())
	    {
	      adj_ix_smooth.push_back(ix);
	    }
	}
      if (!found_smooth)
	adj_ix.push_back(ix);
    }

  double eps = model_->getTolerances().gap;
  double fac = 100.0;
  double lim = fac*model_->getTolerances().neighbour;
  if (has_missing_twin)
    {
#ifdef DEBUG
      std::ofstream of2("missing.g2");
#endif
      // Extend collection of adjacent surface to avoid loosing
      // trimming curves
      for (size_t kj=0; kj<edg_point.size(); ++kj)
	{
#ifdef DEBUG
	  of2 << "400 1 0 4 255 0 0 255" << std::endl;
	  of2 << "1" << std::endl;
	  of2 << edg_point[kj] << std::endl;
#endif
	  double upar, vpar, dist;
	  int clo_ix = -1;
	  Point clo_pnt;
	  for (int ka=0; ka<nmb; ++ka)
	    {
	      if (ka == face_ix)
		continue;
	      shared_ptr<ftSurface> curr_face = model_->getFace(ka);
	      curr_face->closestPoint(edg_point[kj], upar, vpar, clo_pnt,
				      dist, eps);
	      if (dist > lim)
		continue;  // The candidate surface is not in the
	      // neigbourhood

	      // Check if the surface is found already
	      size_t kr;
	      for (kr=0; kr<adj_ix.size(); ++kr)
		if (adj_ix[kr] == ka)
		  break;

	      if (kr == adj_ix.size() && mod_sfs_[ka].get() != surf.get())
		{
		  adj_ix.push_back(ka);
#ifdef DEBUG
		  mod_sfs_[ka]->writeStandardHeader(of2);
		  mod_sfs_[ka]->write(of2);
#endif
		  
		}
	    }
	}
      int stop_break2 = 1;
    }

  // Modify indices in the case of common base surfaces
  for (size_t kj=0; kj<adj_ix.size(); ++kj)
    {
      int ka;
      for (ka=0; ka<adj_ix[kj]; ++ka)
	if (mod_sfs_[ka].get() == mod_sfs_[adj_ix[kj]].get())
	  break;
      if (ka < adj_ix[kj])
	adj_ix[kj] = ka;
    }  

  // Remove duplicates
  for (size_t kj=0; kj<adj_ix.size(); ++kj)
    for (size_t kr=kj+1; kr<adj_ix.size();)
      {
	if (adj_ix[kj] == adj_ix[kr])
	  adj_ix.erase(adj_ix.begin()+kr);
	else
	  ++kr;
      }

  // Modify indices in the case of common base surfaces
  for (size_t kj=0; kj<adj_ix_smooth.size(); ++kj)
    {
      int ka;
      for (ka=0; ka<adj_ix_smooth[kj]; ++ka)
	if (mod_sfs_[ka].get() == mod_sfs_[adj_ix_smooth[kj]].get())
	  break;
      if (ka < adj_ix_smooth[kj])
	adj_ix_smooth[kj] = ka;
    }  

  // Remove duplicates
  for (size_t kj=0; kj<adj_ix_smooth.size(); ++kj)
    for (size_t kr=kj+1; kr<adj_ix_smooth.size();)
      {
	if (adj_ix_smooth[kj] == adj_ix_smooth[kr])
	  adj_ix_smooth.erase(adj_ix_smooth.begin()+kr);
	else
	  ++kr;
      }

  for (size_t kj=0; kj<adj_ix_smooth.size(); ++kj)
    smooth_ix.push_back(make_pair(face_ix, adj_ix_smooth[kj]));

#ifdef DEBUG
  for (size_t kr=0; kr<adj_ix.size(); ++kr)
    {
      mod_sfs_[adj_ix[kr]]->writeStandardHeader(of);
      mod_sfs_[adj_ix[kr]]->write(of);
    }
#endif
  int stop_break = 1;
}

//==========================================================================
 void 
 ModelUpdate::splitIntersectionCurves(vector<int>& ix2,
				      vector<pair<shared_ptr<CurveOnSurface>,
				      shared_ptr<CurveOnSurface> > >& int_seg,
				      vector<int>& int_type,
				      int& nmb_curr, double epsge,
				      double minlen)
//==========================================================================
{
  vector<vector<double> > par(int_seg.size());
  for (int ki=0; ki<nmb_curr; ++ki)
      {
	for (int kj=ki+1; kj<int_seg.size(); ++kj)
	  {
	    vector<pair<double,double> > int_pars;
	    bool use_par_crvs = (int_seg[ki].first->hasParameterCurve() &&
				 int_seg[kj].first->hasParameterCurve());
	    if (int_type[kj] != 3)
	      {
		if (use_par_crvs)
		  intersectParamCurves(int_seg[ki].first->parameterCurve().get(), 
				       int_seg[kj].first->parameterCurve().get(), 
				       epsge, int_pars);
		else
		  intersectParamCurves(int_seg[ki].first.get(), 
				       int_seg[kj].first.get(), 
				       epsge, int_pars);
	      }
	    else
	      {
		// Perform closest point with endpoints
		double endpar[2];
		endpar[0] = int_seg[kj].first->startparam();
		endpar[1] = int_seg[kj].first->endparam();
		for (int kb=0; kb<2; ++kb)
		  {
		    Point pos = int_seg[kj].first->ParamCurve::point(endpar[kb]);
		    double clo_par, clo_dist;
		    Point clo_pt;
		    int_seg[ki].first->closestPoint(pos, 
						    int_seg[ki].first->startparam(),
						    int_seg[ki].first->endparam(),
						    clo_par, clo_pt, clo_dist);
		    if (clo_dist < 1.5*minlen)
		      {
			int_pars.push_back(make_pair(clo_par, endpar[kb]));
		      }
		    
		  }
	      }
	    if (int_pars.size() > 0)
	      {
		// Split intersection curves
		// First represent in independent vectors and add endpoints
		for (size_t kr=0; kr<int_pars.size(); ++kr)
		  {
		    if (use_par_crvs)
		      {
			// Iterate to a better position
			double par1, par2, dist;
			Point ptc1, ptc2;
			double start1 = int_seg[ki].first->startparam();
			double end1 = int_seg[ki].first->endparam();
			double start2 = int_seg[kj].first->startparam();
			double end2 = int_seg[kj].first->endparam();
			ClosestPoint::closestPtCurves(int_seg[ki].first.get(),
						      int_seg[kj].first.get(),
						      start1, end1,
						      start2, end2,
						      int_pars[kr].first, 
						      int_pars[kr].second,
						      par1, par2, dist, ptc1, ptc2);
			par[ki].push_back(par1);
			par[kj].push_back(par2);
		      }
		    else
		      {
			par[ki].push_back(int_pars[kr].first);
			par[kj].push_back(int_pars[kr].second);
		      }
		  }
	      }
	  }
      }

  for (int ki=0; ki<(int)par.size();)
      {
	if (par[ki].size() == 0)
	  {
	    ++ki;
	    continue;
	  }

	// The corresponding intersection curve is oppositely oriented
	shared_ptr<CurveOnSurface> other(int_seg[ki].second->clone());
	other->reverseParameterDirection();

	std::sort(par[ki].begin(), par[ki].end());

	// Extend with endpoints if necessary
	if (fabs(par[ki][0]-int_seg[ki].first->startparam()) < epsge)
	  par[ki][0] = int_seg[ki].first->startparam();
	else
	  par[ki].insert(par[ki].begin(), int_seg[ki].first->startparam());
	if (fabs(int_seg[ki].first->endparam()-par[ki][par[ki].size()-1]) < epsge)
	  par[ki][par[ki].size()-1] = int_seg[ki].first->endparam();
	else
	  par[ki].push_back(int_seg[ki].first->endparam());

	// Split intersection curves in found parameters
	vector<pair<shared_ptr<CurveOnSurface>, shared_ptr<CurveOnSurface> > > sub_int;
	for (size_t kr=1; kr<par[ki].size(); ++kr)
	  {
	    if (par[ki][kr]-par[ki][kr-1] > epsge)
	      {
		shared_ptr<CurveOnSurface> sub1(int_seg[ki].first->subCurve(par[ki][kr-1], par[ki][kr]));
		shared_ptr<CurveOnSurface> sub2(other->subCurve(par[ki][kr-1], par[ki][kr]));
		sub2->reverseParameterDirection();
		sub_int.push_back(make_pair(sub1,sub2));
	      }
	  }
	if (sub_int.size() > 0)
	  {
	    nmb_curr--;
	    int_seg.erase(int_seg.begin()+ki);
	    int_seg.insert(int_seg.end(), 
			   sub_int.begin(), sub_int.end());
	    int other = ix2[ki];
	    ix2.erase(ix2.begin()+ki);
	    for (size_t kr=0; kr<sub_int.size(); ++kr)
	      ix2.push_back(other);
	    par.erase(par.begin()+ki);
	  }
	else
	  ki++;
      }
#ifdef DEBUG
  std::ofstream of("split_int_cvs.g2");
  for (size_t kr=0; kr<int_seg.size(); ++kr)
    {
      of << "100 1 0 4 255 0 0 255" << std::endl;
      int_seg[kr].first->spaceCurve()->write(of);
      of << "100 1 0 4 0 255 0 255" << std::endl;
      int_seg[kr].second->spaceCurve()->write(of);
    }
#endif
}

//==========================================================================
 void 
   ModelUpdate::cleanIntersectionCurves(shared_ptr<ParamSurface> surf,
					vector<int>& ix2,
					vector<pair<shared_ptr<CurveOnSurface>,
					shared_ptr<CurveOnSurface> > >& int_seg,
					double epsge, double minlen, 
					double minlen2)
//==========================================================================
{
  // First fetch the boundary loop of the initial surface
  CurveLoop bd_loop = SurfaceTools::outerBoundarySfLoop(surf, 0.0);
  vector<shared_ptr<ParamCurve> > bd_cvs = bd_loop.getCurves();

  // Check endpoints of intersection curves for coincidence
  size_t ki, kj;
  for (ki=0; ki<int_seg.size(); )
    {
      // First check with itself
      double cv_len = int_seg[ki].first->estimatedCurveLength();
      Point pos[2];
      int_seg[ki].first->point(pos[0], int_seg[ki].first->startparam());
      int_seg[ki].first->point(pos[1], int_seg[ki].first->endparam());
      double dist = pos[0].dist(pos[1]);
      if (cv_len < epsge)
	{
	  // Remove curve
	  ix2.erase(ix2.begin()+ki);
	  int_seg.erase(int_seg.begin()+ki);
	  continue;
	}
      else if (dist < minlen || (dist < minlen2 && cv_len > minlen2))
	{
	  ++ki;
	  continue;
	}

      int ka;
      for (ka=0; ka<2; ++ka)
	{
	  // First check if the endpoint coincides with the endpoint of 
	  // another intersection curve
	  size_t kj;
	  for (kj=0; kj<int_seg.size(); kj++)
	    {
	      if (ki == kj)
		continue;
	      double endpar[2];
	      endpar[0] = int_seg[kj].first->startparam();	
	      endpar[1] = int_seg[kj].first->endparam();
	      int kb;
	      for (kb=0; kb<2; ++kb)
		{
		  Point pos2 = int_seg[kj].first->ParamCurve::point(endpar[kb]);
		  double dist2 = pos[ka].dist(pos2);
		  if (dist2 < minlen || (dist2 < minlen2 && cv_len > minlen2))
		    break;
		}
	      if (kb < 2)
		break;
	    }
	  if (kj < int_seg.size())
	    continue;

	  // Check if the current intersection curve ends up in a 
	  // boundary curve
	  size_t kr;
	  for (kr=0; kr<bd_cvs.size(); ++kr)
	    {
	      double cpar, dist;
	      Point ptc;
	      bd_cvs[kr]->closestPoint(pos[ka], bd_cvs[kr]->startparam(),
				       bd_cvs[kr]->endparam(), cpar, ptc,
				       dist);
	      if (dist < minlen /*|| (dist < minlen2 && cv_len > minlen2)*/)
		break;
	    }
	  if (kr == bd_cvs.size())
	    break;
	}

      if (ka < 2)
	{
	  // Loose end found, remove curve
	  ix2.erase(ix2.begin()+ki);
	  int_seg.erase(int_seg.begin()+ki);
	}
      else
	++ki;
    }
}

//==========================================================================
 void 
   ModelUpdate::getTangentialInt(shared_ptr<ParamSurface> surf1, 
				 shared_ptr<ParamSurface> surf2, 
				 vector<shared_ptr<CurveOnSurface> >& int_segments1, 
				 vector<shared_ptr<CurveOnSurface> >& int_segments2, 
				 double epsge)
//==========================================================================
 {
   // Check
   DirectionCone tcone1_1 = surf1->tangentCone(true);
   DirectionCone tcone1_2 = surf1->tangentCone(false);
   DirectionCone tcone2_1 = surf2->tangentCone(true);
   DirectionCone tcone2_2 = surf2->tangentCone(false);

   // Fetch boundary curves
   double deg_eps = 1.0e-12;
   CurveLoop loop1 = SurfaceTools::outerBoundarySfLoop(surf1, deg_eps);
   CurveLoop loop2 = SurfaceTools::outerBoundarySfLoop(surf2, deg_eps);

   // Intersect boundary curves
   vector<double> parvals;
   for (int ki=0; ki<loop1.size(); ++ki)
     {
       shared_ptr<ParamCurve> crv = loop1[ki];
       vector<pair<double, Point> > intpts;
       vector<int> pretop;
       vector<pair<pair<double,Point>, pair<double,Point> > > intcvs;
       intersectParamCurveSurf(crv.get(), surf2.get(), epsge, intpts,
			       pretop, intcvs);
       if (intpts.size() > 0)
	 {
	   shared_ptr<CurveOnSurface> sfcv =
	     dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
	   if (sfcv.get())
	     {
	       for (size_t kr=0; kr<intpts.size(); ++kr)
		 {
		   Point par = sfcv->faceParameter(intpts[kr].first);
		   parvals.insert(parvals.end(), par.begin(), par.end());
		   parvals.insert(parvals.end(), 
				  intpts[kr].second.begin(),
				  intpts[kr].second.end());
		 }
	     }
	 }
     }
   for (int ki=0; ki<loop2.size(); ++ki)
     {
       shared_ptr<ParamCurve> crv = loop2[ki];
       vector<pair<double, Point> > intpts;
       vector<int> pretop;
       vector<pair<pair<double,Point>, pair<double,Point> > > intcvs;
       intersectParamCurveSurf(crv.get(), surf1.get(), epsge, intpts,
			       pretop, intcvs);
       if (intpts.size() > 0)
	 {
	   shared_ptr<CurveOnSurface> sfcv =
	     dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
	   if (sfcv.get())
	     {
	       for (size_t kr=0; kr<intpts.size(); ++kr)
		 {
		   parvals.insert(parvals.end(), 
				  intpts[kr].second.begin(),
				  intpts[kr].second.end());
		   Point par = sfcv->faceParameter(intpts[kr].first);
		   parvals.insert(parvals.end(), par.begin(), par.end());
		 }
	     }
	 }
     }

   if (parvals.size() > 0)
     {
       // Try to march between pairwise intersection points
       // Prepare for calling SISL
       shared_ptr<SplineSurface> tmp_sf1, tmp_sf2;
       SplineSurface *sf1 = surf1->getSplineSurface();
       SplineSurface *sf2 = surf2->getSplineSurface();
       if (!surf1)
	 {
	   tmp_sf1 = shared_ptr<SplineSurface>(sf1->asSplineSurface());
	   sf1 = tmp_sf1.get();
	 }
       if (!surf2)
	 {
	   tmp_sf2 = shared_ptr<SplineSurface>(sf2->asSplineSurface());
	   sf2 = tmp_sf2.get();
	 }
       if (!(sf1 && sf2))
	 THROW("No underlying spline surface");

       SISLSurf* sisl_sf1 = GoSurf2SISL(*sf1);
       SISLSurf* sisl_sf2 = GoSurf2SISL(*sf2);
       double maxstep = (double)0;
       int makecurv = 2; // Make both geometric and parametric curves.
       int draw = 0;

#ifdef DEBUG
       std::ofstream of("tang_int_cvs.g2");
#endif
       for (int ka=0; ka<(int)parvals.size();)
	 {
	   int nmbvals = (int)parvals.size();
	   int kb;
	   for (kb=ka+4; kb<(int)parvals.size(); kb+=4)
	     {
	       double epar1[4];
	       double epar2[4];
	       epar1[0] = parvals[ka];
	       epar1[1] = parvals[ka+1];
	       epar1[2] = parvals[kb];
	       epar1[3] = parvals[kb+1];
	       epar2[0] = parvals[ka+2];
	       epar2[1] = parvals[ka+3];
	       epar2[2] = parvals[kb+2];
	       epar2[3] = parvals[kb+3];
	       SISLIntcurve *intcv = newIntcurve(2, 2, 2, epar1, epar2, 1);

	       int status;
	       s1310(sisl_sf1, sisl_sf2, intcv, epsge, maxstep, makecurv,
		     draw, &status);
#ifdef DEBUG
	       std::cout << "status: " << status << std::endl;
#endif
	       if (status >= 0 && intcv->pgeom != NULL)
		 {
		   shared_ptr<SplineCurve> pcurve1(SISLCurve2Go(intcv->ppar1));
		   shared_ptr<SplineCurve> pcurve2(SISLCurve2Go(intcv->ppar2));
		   shared_ptr<SplineCurve> space_curve1(SISLCurve2Go(intcv->pgeom));
		   shared_ptr<SplineCurve> space_curve2(space_curve1->clone());

#ifdef DEBUG
		   space_curve1->writeStandardHeader(of);
		   space_curve1->write(of);
#endif
		   // Make sure that the curves are k-regular
		   pcurve1->makeKnotStartRegular();
		   pcurve1->makeKnotEndRegular();
		   pcurve2->makeKnotStartRegular();
		   pcurve2->makeKnotEndRegular();
		   space_curve1->makeKnotStartRegular();
		   space_curve1->makeKnotEndRegular();
		   space_curve2->makeKnotStartRegular();
		   space_curve2->makeKnotEndRegular();

		   // We make sure the intersection cvs have the right direction (area below sfs
		   // to be trimmed away).
		   BoundedUtils::consistentIntersectionDir(*pcurve1, 
							   *space_curve1, 
							   *surf1,
							   *pcurve2, 
							   *space_curve2, *
							   surf2, epsge);
		   BoundedUtils::consistentIntersectionDir(*pcurve2, 
							   *space_curve2, 
							   *surf2,
							   *pcurve1, 
							   *space_curve1, 
							   *surf1, epsge);
		   int_segments1.push_back(shared_ptr<CurveOnSurface>
					   (new CurveOnSurface(surf1, pcurve1, 
							       space_curve1, 
							       false)));
		   int_segments2.push_back(shared_ptr<CurveOnSurface>
					   (new CurveOnSurface(surf2, pcurve2, 
							       space_curve2, 
							       false)));

		   // Do not reuse points
		   parvals.erase(parvals.begin()+kb, parvals.begin()+kb+4);
		   parvals.erase(parvals.begin()+ka, parvals.begin()+ka+4);
		 }
	       if (intcv != NULL)
		 {
		   intcv->epar1 = NULL;
		   intcv->epar2 = NULL;
		   freeIntcurve(intcv);
		 }
	       if (status >= 0)
		 {
		   break;
		 }
	     }
	   if (kb >= nmbvals)
	     ka += 4;
	 }
       if (sisl_sf1)
	 freeSurf(sisl_sf1);
       if (sisl_sf2)
	 freeSurf(sisl_sf2);
     }
       
   int stop_break = 1;
 }

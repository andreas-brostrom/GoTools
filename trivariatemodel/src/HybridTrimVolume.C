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

#include "GoTools/trivariatemodel/HybridTrimVolume.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeModelCreator.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"
#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/GapRemoval.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include <fstream>
#include <cstdlib>

#define DEBUG

using std::vector;
using std::set;
using std::make_pair;
using std::pair;

using namespace Go;

//==========================================================================
HybridTrimVolume::HybridTrimVolume(shared_ptr<SurfaceModel> model,
				   int material)
//==========================================================================
{
  model_ = model;
  material_ = material;

  // Unchanged backup model
  init_model_ = shared_ptr<SurfaceModel>(new SurfaceModel(*model_));
}

//==========================================================================
HybridTrimVolume::~HybridTrimVolume()
//==========================================================================
{

}

//==========================================================================
shared_ptr<VolumeModel> 
HybridTrimVolume::fetchVolModel(bool create_degen, bool refine_sharp)
//==========================================================================
{
  shared_ptr<VolumeModel> result;

  SurfaceModelUtils::limitUnderlyingSurfaces(model_);

  // Ensure that the initial input model is represented as a solid
  shared_ptr<Body> init_body;
  Body *model_body = init_model_->getBody();
  if ((!model_body) && init_model_->nmbBoundaries() == 0)
    init_body = shared_ptr<Body>(new Body(init_model_));

  // Identify rotational axis and classify all faces into three groups,
  // those that are rotational with respect, those that are rotational with
  // trimming curves destroying the rotational pattern and the rest
  vector<shared_ptr<ftSurface> > rotational_faces;
  vector<shared_ptr<ftSurface> > rotational_faces2;
  vector<shared_ptr<ftSurface> > other_faces;
  Point axis, centre, vec;
  double angle;
  bool found_axis = identifyRotationalAxis(centre, axis, vec, angle,
					   rotational_faces, 
					   rotational_faces2, other_faces);
#ifdef DEBUG
  std::ofstream ofrot("rot_faces1.g2");
  for (size_t ki=0; ki<rotational_faces.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = rotational_faces[ki]->surface();
      sf->writeStandardHeader(ofrot);
      sf->write(ofrot);
    }
  std::ofstream ofrot2("rot_faces2.g2");
  for (size_t ki=0; ki<rotational_faces2.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = rotational_faces2[ki]->surface();
      sf->writeStandardHeader(ofrot2);
      sf->write(ofrot2);
    }
  std::ofstream ofoth("other_faces.g2");
  for (size_t ki=0; ki<other_faces.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = other_faces[ki]->surface();
      sf->writeStandardHeader(ofoth);
      sf->write(ofoth);
    }
#endif

#ifdef DEBUG
  std::ofstream ofrotn1("rot_faces_ext1.g2");
  for (size_t ki=0; ki<rotational_faces.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = rotational_faces[ki]->surface();
      sf->writeStandardHeader(ofrotn1);
      sf->write(ofrotn1);
    }
  for (size_t ki=0; ki<rotational_faces2.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = rotational_faces2[ki]->surface();
      sf->writeStandardHeader(ofrotn1);
      sf->write(ofrotn1);
    }
#endif
  // Remove non-rotational inner trimming loops from rotational surfaces
  // and store for later use. Update topology accordingly
  removeInnerTrim(rotational_faces2, centre, axis, vec);

  // Make sure that the outer trimming curve is consistent with a
  // rotational model
  cleanOuterTrim(rotational_faces2, centre, axis, vec);

  // Create rotational model
  tpTolerances tptol = model_->getTolerances();
  shared_ptr<SurfaceModel> rotational_mod(new SurfaceModel(tptol.gap,
							   tptol.gap,
							   tptol.neighbour,
							   tptol.kink,
							   tptol.bend,
							   rotational_faces,
							   true));
  rotational_mod->append(rotational_faces2);
  int nmb_bd = rotational_mod->nmbBoundaries();
#ifdef DEBUG
  bool OK = rotational_mod->checkShellTopology();
  std::cout << "Rotational topology " << OK << std::endl;
  std::cout << "Number of boundaries " << nmb_bd << std::endl;
#endif

  shared_ptr<Body> rotational_body;
  if (nmb_bd == 0)
    {
      // Expected. Mark that the shell is a solid to simplify inside testing
      rotational_body = shared_ptr<Body>(new Body(rotational_mod));
    }

  // Sort non-rational faces into connected groups
  shared_ptr<SurfaceModel> other_mod(new SurfaceModel(tptol.gap,
						      tptol.gap,
						      tptol.neighbour,
						      tptol.kink,
						      tptol.bend,
						      other_faces));
  vector<shared_ptr<SurfaceModel> > other_sub =
    other_mod->getConnectedModels();

  vector<pair<shared_ptr<SurfaceModel>, bool> > classified_sub(other_sub.size());
  // Classify excluded sub models with regard to whether or not they
  // are internal to the rotational model
  for (size_t ki=0; ki<other_sub.size(); ++ki)
    {
      bool inside = checkSubModelInternal(other_sub[ki], rotational_mod);
      classified_sub[ki] = make_pair(other_sub[ki], inside);
    }

#ifdef DEBUG
  std::ofstream ofrotn2("rot_faces_ext2.g2");
  for (size_t ki=0; ki<rotational_mod->nmbEntities(); ++ki)
    {
      shared_ptr<ParamSurface> sf = rotational_mod->getSurface(ki);
      sf->writeStandardHeader(ofrotn2);
      sf->write(ofrotn2);
    }
#endif
  VolumeModelCreator::createRotationalModel(rotational_mod,
					    rotational_vol_);

  if (rotational_vol_.get())
    {
      // Check adjacency
      int nmb_vol_rot = rotational_vol_->nmbEntities();
      for (int ka=0; ka<nmb_vol_rot; ++ka)
	{
	  shared_ptr<ftVolume> tmp_vol = rotational_vol_->getBody(ka);
	  vector<ftVolume*> next_vol;
	  tmp_vol->getAdjacentBodies(next_vol);
	  int adj_break = 1;
	}

      // Expects only outer boundary
      vector<shared_ptr<ftSurface> > rot_bd_faces = 
	rotational_vol_->getBoundaryFaces();
      vector<shared_ptr<ftSurface> > rot_bd_faces2(rot_bd_faces.size()); 
      for (size_t ki=0; ki<rot_bd_faces.size(); ++ki)
	rot_bd_faces2[ki] = 
	  shared_ptr<ftSurface>(new ftSurface(rot_bd_faces[ki]->surface(), 
					       ki));
	rotational_shell_ = 
	   shared_ptr<SurfaceModel>(new SurfaceModel(tptol.gap,
						     tptol.gap,
						     tptol.neighbour,
						     tptol.kink,
						     tptol.bend,
						     rot_bd_faces2));

#ifdef DEBUG
      std::ofstream ofvol("rot_vol.g2");
      int nmb_vol = rotational_vol_->nmbEntities();
      for (int ka=0; ka<nmb_vol; ++ka)
	{
	  shared_ptr<ParamVolume> tmp = rotational_vol_->getVolume(ka);
	  tmp->writeStandardHeader(ofvol);
	  tmp->write(ofvol);
	}
      std::ofstream ofvolbd("rot_vol_faces.g2");
      for (size_t ki=0; ki<rot_bd_faces.size(); ++ki)
	{
	  shared_ptr<ParamSurface> tmp_sf = rot_bd_faces[ki]->surface();
	  tmp_sf->writeStandardHeader(ofvolbd);
	  tmp_sf->write(ofvolbd);
	}
#endif
    }

  // Fetch non-connected trimming loops from excluded sub models
#ifdef DEBUG
  std::ofstream ofedg("trim_edgs.g2");
#endif
  vector<shared_ptr<ParamSurface> > sfs_to_trim;
  vector<vector<shared_ptr<CurveOnSurface> > > trim_cvs;
  for (size_t ki=0; ki<other_sub.size(); ++ki)
    {
      vector<vector<shared_ptr<ftEdge> > > free_edgs;
      fetchFreeLoops(other_sub[ki], free_edgs, sfs_to_trim, trim_cvs);
      for (size_t kj=0; kj<free_edgs.size(); ++kj)
	for (size_t kr=0; kr<free_edgs[kj].size(); ++kr)
	  {
	    shared_ptr<ParamCurve> cv = free_edgs[kj][kr]->geomCurve();
#ifdef DEBUG
	    if (cv.get())
	      {
		shared_ptr<CurveOnSurface> sfcv = 
		  dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
		if (sfcv.get() && sfcv->spaceCurve().get())
		  {
		    sfcv->spaceCurve()->writeStandardHeader(ofedg);
		    sfcv->spaceCurve()->write(ofedg);
		  }
		else
		  {
		    cv->writeStandardHeader(ofedg);
		    cv->write(ofedg);
		  }
	      }
#endif
	  }
    }

  vector<vector<shared_ptr<ParamSurface> > > associated_surface;
  vector<vector<shared_ptr<ftSurface> > > added_face;
  vector<shared_ptr<ftSurface> > trim_faces;
  performTrimming(sfs_to_trim, trim_cvs, classified_sub, associated_surface,
		  added_face, trim_faces);

#ifdef DEBUG
  for (size_t ki=0; ki<classified_sub.size(); ++ki)
    {
      if (classified_sub[ki].second == true)
	continue;

      // int nmb_bd_sub = classified_sub[ki].first->nmbBoundaries();
      // std::cout << "Sub model " << ki << ", no. of boundaries: " << nmb_bd_sub << std::endl;
      std::ofstream ofout("out_sub.g2");
      int nmb = classified_sub[ki].first->nmbEntities();
      for (int ka=0; ka<nmb; ++ka)
	{
	  shared_ptr<ParamSurface> tmp_sf = 
	    classified_sub[ki].first->getSurface(ka);
	  tmp_sf->writeStandardHeader(ofout);
	  tmp_sf->write(ofout);
	}
      for (size_t kr=0; kr<added_face[ki].size(); ++kr)
	{
	  shared_ptr<ParamSurface> tmp_sf = 
	    added_face[ki][kr]->surface();
	  tmp_sf->writeStandardHeader(ofout);
	  tmp_sf->write(ofout);
	}
	  
      // if (nmb_bd_sub > 0)
      // 	{
      // 	  std::ofstream of_bd("submod_bd.g2");
      // 	  vector<shared_ptr<ftEdge> > bd_edgs = 
      // 	    classified_sub[ki].first->getBoundaryEdges();
      // 	  for (size_t kr=0; kr<bd_edgs.size(); ++kr)
      // 	    {
      // 	      shared_ptr<ParamCurve> tmp_cv = bd_edgs[kr]->geomCurve();
      // 	      shared_ptr<ParamCurve> tmp_cv2(tmp_cv->geometryCurve());
      // 	      tmp_cv2->writeStandardHeader(of_bd);
      // 	      tmp_cv2->write(of_bd);
      // 	    }
      // 	}
	int stop_break_out = 1;
	}
#endif

  vector<shared_ptr<ftSurface> > bd_faces = 
    rotational_vol_->getBoundaryFaces();
  vector<shared_ptr<ftSurface> > new_faces;
  vector<shared_ptr<ftSurface> > bd_removed;
  for (size_t ki=0; ki<trim_faces.size(); ++ki)
    {
      shared_ptr<BoundedSurface> trim_sf = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(trim_faces[ki]->surface());
      for (size_t kj=0; kj<bd_faces.size(); ++kj)
	{
	  shared_ptr<ParamSurface> bd_sf = bd_faces[kj]->surface();
	  shared_ptr<BoundedSurface> bd2_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(bd_sf);
	  if (bd2_sf.get())
	    bd_sf = bd2_sf->underlyingSurface();
	  if (trim_sf->underlyingSurface().get() == bd_sf.get())
	    {
	      // Replace block boundary surface. First find associated shell
	      Body* body = bd_faces[kj]->getBody();
	      shared_ptr<SurfaceModel> shell;
	      if (body)
		shell = body->getOuterShell();
	      if (shell.get())
		{
#ifdef DEBUG
		  std::ofstream shout1("shell_mod1.g2");
		  int nmb1 = shell->nmbEntities();
		  for (int ka1=0; ka1<nmb1; ++ka1)
		    {
		      shared_ptr<ParamSurface> tmp_sf = 
			shell->getSurface(ka1);
		      tmp_sf->writeStandardHeader(shout1);
		      tmp_sf->write(shout1);
		    }
#endif
		  size_t kr;
		  for (kr=0; kr<bd_removed.size(); ++kr)
		    if (bd_removed[kr].get() == bd_faces[kj].get())
		      break;
		  if (kr == bd_removed.size())
		    {
		      shell->removeFace(bd_faces[kj]);
		      bd_removed.push_back(bd_faces[kj]);
		    }
		  vector<shared_ptr<ftSurface> > new_bd_faces;
		  new_bd_faces.push_back(trim_faces[ki]);
		  shell->append(new_bd_faces);
		  new_faces.push_back(trim_faces[ki]);
#ifdef DEBUG
		  std::ofstream shout2("shell_mod2.g2");
		  int nmb2 = shell->nmbEntities();
		  for (int ka1=0; ka1<nmb2; ++ka1)
		    {
		      shared_ptr<ParamSurface> tmp_sf = 
			shell->getSurface(ka1);
		      tmp_sf->writeStandardHeader(shout2);
		      tmp_sf->write(shout2);
		    }
#endif
		  int stop_shell = 1;
		}
	    }
	}
    }

  // Add associated inner trimming sub models
  for (size_t kr=0; kr<associated_surface.size(); ++kr)
    {
      if (!classified_sub[kr].second)
	continue;   // External sub model

      vector<shared_ptr<ftSurface> > boundary_faces;
      for (size_t kj=0; kj<new_faces.size(); ++kj)
	{
	  shared_ptr<ParamSurface> bd_sf = new_faces[kj]->surface();
	  shared_ptr<BoundedSurface> bd_sf2 = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(bd_sf);
	  if (bd_sf2.get())
	    bd_sf = bd_sf2->underlyingSurface();

	  for (size_t kh=0; kh<associated_surface[kr].size(); ++kh)
	    {
	      if (associated_surface[kr][kh].get() == bd_sf.get())
		boundary_faces.push_back(new_faces[kj]);
	    }
	}
      insertTrimModel(classified_sub[kr].first, boundary_faces);
    }
	  
  // Add external sub models. 
  // Perform block structuring if possible
  for (size_t kr=0; kr<associated_surface.size(); ++kr)
    {
      if (classified_sub[kr].second)
	continue;   // Sub model representing trimming

      shared_ptr<VolumeModel> sub_volmod = 
	subModel2Trivariate(classified_sub[kr].first, added_face[kr]);

      rotational_vol_->append(sub_volmod, false);
#ifdef DEBUG
      std::ofstream ofvol("curr_vol.g2");
      int nmb_vols = rotational_vol_->nmbEntities();
      for (int ka=0; ka<nmb_vols; ++ka)
	{
	  shared_ptr<ParamVolume> curr_vol2 = rotational_vol_->getVolume(ka);
	  curr_vol2->writeStandardHeader(ofvol);
	  curr_vol2->write(ofvol);
	}
#endif
    }
  
  rotational_vol_->resetBoundarySfs();

  
  VolumeModelFileHandler filehandler;
  std::ofstream outfile("trimvolmod.g22");
  filehandler.writeStart(outfile);
  filehandler.writeHeader("Trimmed volume model", outfile);
  filehandler.writeVolumeModel(*rotational_vol_, outfile);
  filehandler.writeEnd(outfile);

  return rotational_vol_;
}

//==========================================================================
bool 
HybridTrimVolume::identifyRotationalAxis(Point& centre, Point& axis, 
					 Point& vec, double& angle,
					 vector<shared_ptr<ftSurface> >& rotational_faces, 
					 vector<shared_ptr<ftSurface> >& rotational_faces2, 
					 vector<shared_ptr<ftSurface> >& other_faces)
//==========================================================================
{
  rotational_faces.clear();
  other_faces.clear();

  double eps = model_->getTolerances().neighbour;
  double bend = model_->getTolerances().bend;
  double kink = model_->getTolerances().kink;

  int nmb = model_->nmbEntities();
  vector<Point> all_centre;
  vector<Point> all_axis;
  vector<vector<shared_ptr<ftSurface> > > rot_faces;
  vector<vector<int> > rot_type;
  vector<vector<pair<Point,double> > > slices;
  vector<vector<double> > radius;
  for (int ki=0; ki<nmb; ++ki)
    {
      // Fetch surface/underlying surface
      shared_ptr<ftSurface> face = model_->getFace(ki);
      shared_ptr<ParamSurface> surf = face->surface();
      Point curr_centre, curr_axis, curr_vec;
      double ang;
      double rad = -1;
      int rotational = surf->isAxisRotational(curr_centre, curr_axis,
					       curr_vec, ang);
      if (rotational == 0)
	{
	  shared_ptr<BoundedSurface> bd_surf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
	  if (bd_surf.get())
	    {
	      shared_ptr<ParamSurface> surf2 = bd_surf->underlyingSurface();

	      rotational = surf2->isAxisRotational(curr_centre, curr_axis,
						   curr_vec, ang);
	      if (rotational)
		rotational = 2;   // Trimming

	      // Must adjust curr_vec and ang according to the bounding
	      // domain of the trimmed surface
	      int stop_break = 1;
	    }
	}

      if (rotational)
	{
	  // Check if a new rotational axis is found
	  curr_axis.normalize();
	  size_t kj;
	  for (kj=0; kj<all_centre.size(); ++kj)
	    {
	      Point vec2 = curr_centre - all_centre[kj];
	      Point vec3 = vec2 - (vec2*curr_axis)*curr_axis;
	      double dd = vec3.length();
	      double ang2 = curr_axis.angle(all_axis[kj]);
	      if (dd < eps && M_PI-ang2 < bend)
		{
		  // The axes are oppositely oriented. Adjust the start vector
		  Array<double,3> tmp_vec(curr_vec[0], curr_vec[1], curr_vec[2]);
		  MatrixXD<double, 3> mat;
		  mat.setToRotation(ang, curr_axis[0], 
				    curr_axis[1], curr_axis[2]);  // Rotate the 
		  // start vector the angle ang around curr_axis
		  Array<double,3> tmp_vec2 = mat*tmp_vec;
		  curr_vec = Point(tmp_vec2[0], tmp_vec2[1], tmp_vec2[2]);
		}
		  
	      ang2 = std::min(ang2, M_PI-ang2);
	      ElementarySurface *elem = face->surface()->elementarySurface();
	      if (elem)
		{
		  RectDomain dom = face->surface()->containingDomain();
		  double par1[2], par2[2];
		  if (elem->isSwapped())
		    {
		      par1[0] = dom.umin();
		      par2[0] = dom.umax();
		      par1[1] = par2[1] = 0.5*(dom.vmin()+dom.vmax());
		    }
		  else
		    {
		      par1[0] = par2[0] = 0.5*(dom.umin()+dom.umax());
		      par1[1] = dom.vmin();
		      par2[1] = dom.vmax();
		    }
		  double curr_rad1 = elem->radius(par1[0], par1[1]);
		  double curr_rad2 = elem->radius(par2[0], par2[1]);
		  rad = 0.5*(curr_rad1+curr_rad2);
		}

	      if (dd < eps && ang2 < bend)
		{
		  rot_faces[kj].push_back(face);
		  rot_type[kj].push_back(rotational);
		  slices[kj].push_back(make_pair(curr_vec, ang));
		  radius[kj].push_back(rad);
		  break;
		}
	    }
	  if (kj == all_centre.size())
	    {
	      all_centre.push_back(curr_centre);
	      all_axis.push_back(curr_axis);
	      vector<shared_ptr<ftSurface> > curr_faces;
	      curr_faces.push_back(face);
	      rot_faces.push_back(curr_faces);
	      vector<int> curr_type;
	      curr_type.push_back(rotational);
	      rot_type.push_back(curr_type);
	      vector<pair<Point,double> > curr_slices;
	      curr_slices.push_back(make_pair(curr_vec, ang));
	      slices.push_back(curr_slices);
	      vector<double> curr_radius;
	      curr_radius.push_back(rad);
	      radius.push_back(curr_radius);
	    }
	}
      else
	other_faces.push_back(face);
    }

  // Check if there is one dominant rotational axis
  // // Too simple!
  // int max_nmb = 0, nmb_max = 0, idx_max = -1;
  // for (size_t kj=0; kj<rot_faces.size(); ++kj)
  //   {
  //     int nmb_faces = (int)rot_faces[kj].size();
  //     if (nmb_faces > max_nmb)
  // 	{
  // 	  max_nmb = nmb_faces;
  // 	  nmb_max = 1;
  // 	  idx_max = (int)kj;
  // 	}
  //     else if (nmb_faces == max_nmb)
  // 	nmb_max++;
  //   }
  
  // if (idx_max < 0 || nmb_max > 1)
  //   return false;

  int idx_max = -1;
  double max_rad = 0.0;
  for (size_t kj=0; kj<radius.size(); ++kj)
    {
      double curr_max = 0.0;
      for (size_t kr=0; kr<radius[kj].size(); ++kr)
	curr_max = std::max(curr_max, radius[kj][kr]);

      if (curr_max > max_rad)
	{
	  max_rad = curr_max;
	  idx_max = (int)kj;
	}
    }

  if (idx_max < 0)
    return false;

  // Sort faces into the ones agreeing with the found rotational axis
  // and the other
  for (size_t kj=0; kj<rot_faces.size(); ++kj)
    {
      if ((int)kj == idx_max)
	{
	  for (size_t kr=0; kr<rot_faces[kj].size(); ++kr)
	    {
	      if (rot_type[kj][kr] == 1)
		rotational_faces.push_back(rot_faces[kj][kr]);
	      else
		rotational_faces2.push_back(rot_faces[kj][kr]);
	    }
	}
      else
	other_faces.insert(other_faces.end(), rot_faces[kj].begin(),
			   rot_faces[kj].end());
    }
  centre = all_centre[idx_max];
  axis = all_axis[idx_max];

#ifdef DEBUG
  std::ofstream of("rot_faces.g2");
  for (size_t kj=0; kj<rotational_faces.size(); ++kj)
    {
      shared_ptr<ParamSurface> tmp_sf = rotational_faces[kj]->surface();
      tmp_sf->writeStandardHeader(of);
      tmp_sf->write(of);
    }
#endif
  // Compute rotational sector
  // Remove slice redundancies
  for (size_t kj=0; kj<slices[idx_max].size(); ++kj)
    {
      size_t kr;
      for (kr=kj+1; kr<slices[idx_max].size(); )
	{
	  if (fabs(radius[idx_max][kj]-radius[idx_max][kr]) > eps)
	    {
	      ++kr;
	      continue;  // Not the same surface
	    }

	  double ang2 = 
	    slices[idx_max][kj].first.angle(slices[idx_max][kr].first);
	  double ang3 = slices[idx_max][kj].second;
	  double ang4 = slices[idx_max][kr].second;
	  if ((ang2 <= kink && (fabs(ang3 - ang4) < kink || ang4 < ang3)) ||
	      ang3-ang2-ang4 > -kink)
	    {
	      slices[idx_max].erase(slices[idx_max].begin()+kr);
	      radius[idx_max].erase(radius[idx_max].begin()+kr);
	    }
	  else if ((ang2 <= kink && ang3 < ang4) || 
		   ang4-ang2-ang3 > -kink)
	    {
	      std::swap(slices[idx_max][kj], slices[idx_max][kr]);
	      std::swap(radius[idx_max][kj], radius[idx_max][kr]);
	      slices[idx_max].erase(slices[idx_max].begin()+kr);
	      radius[idx_max].erase(radius[idx_max].begin()+kr);
	      kr = kj+1;
	    }
	  else 
	    ++kr;
	}
    }

  // Join adjacent slices
  for (size_t kj=0; kj<slices[idx_max].size(); ++kj)
    {
      size_t kr;
      for (kr=kj+1; kr<slices[idx_max].size(); )
	{
	  if (fabs(radius[idx_max][kj]-radius[idx_max][kr]) > eps)
	    {
	      ++kr;
	      continue;  // Not the same surface
	    }

	  double ang2 = 
	    slices[idx_max][kj].first.angle(slices[idx_max][kr].first);
	  double ang3 = slices[idx_max][kj].second;
	  double ang4 = slices[idx_max][kr].second;
	  if (fabs(ang2-ang3) < kink)
	    {
	      slices[idx_max][kj].second += ang4;
	      slices[idx_max].erase(slices[idx_max].begin()+kr);
	      radius[idx_max].erase(radius[idx_max].begin()+kr);
	    }
	  else if (fabs(ang2-ang4) < kink)
	    {
	      std::swap(slices[idx_max][kj], slices[idx_max][kr]);
	      slices[idx_max][kj].second += ang3;
	      slices[idx_max].erase(slices[idx_max].begin()+kr);
	      std::swap(radius[idx_max][kj], radius[idx_max][kr]);
	      radius[idx_max].erase(radius[idx_max].begin()+kr);
	      kr = kj+1;
	    }
	  else
	    ++kr;
	}
    }
  
  // Find largest gap

  // Set rotational sector
  angle = 0.0;
  vec = slices[idx_max][0].first;
  double rot_ang = slices[idx_max][0].second;
  for (size_t kj=1; kj<slices[idx_max].size(); ++kj)
    {
      double ang2 = 
	slices[idx_max][0].first.angle(slices[idx_max][kj].first);
      angle += ang2;
      if (slices[idx_max][kj].second > rot_ang)
	{
	  rot_ang = slices[idx_max][kj].second;
	  vec = slices[idx_max][kj].first;
	}
    }
  angle += slices[idx_max][slices[idx_max].size()-1].second;
  angle = std::max(angle, rot_ang);
  angle = std::min(angle, 2.0*M_PI);

  // Make sure that the vector is orthogonal to the axis
  Point vec2 = vec % axis;
  vec = axis % vec2;

  return true;
}

//==========================================================================
void
HybridTrimVolume::cleanOuterTrim(vector<shared_ptr<ftSurface> >& faces, 
				 const Point& centre, const Point& axis, 
				 const Point& vec)
//==========================================================================
{
  for (size_t ki=0; ki<faces.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = faces[ki]->surface();
      shared_ptr<BoundedSurface> bdsf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
      if (!bdsf.get())
	{
	  return; // Surface has iso-parametric boundary curves
	}

      shared_ptr<CurveLoop> newloop;
      vector<vector<shared_ptr<ftEdge> > > removed_edgs;
      bool consistent = checkOuterLoop(faces[ki], centre, axis, vec,
				       newloop, removed_edgs);
      if (consistent)
	continue;

      // Replace loop and update topology
      //bool removed = model_->removeFace(faces[ki]);
      faces[ki]->isolateFace();
      shared_ptr<CurveLoop> bdloop = bdsf->loop(0);
      bdloop->swap(*newloop);
      shared_ptr<ftSurface> newface(new ftSurface(surf, 0));
      //model_->append(newface);
      faces[ki] = newface;

      // Store removed edges
      edge_seqs_.insert(edge_seqs_.end(), removed_edgs.begin(),
			removed_edgs.end());
    }
}
//==========================================================================
void
HybridTrimVolume::removeInnerTrim(vector<shared_ptr<ftSurface> >& faces, 
				   const Point& centre, const Point& axis, 
				   const Point& vec)
//==========================================================================
{
  for (size_t ki=0; ki<faces.size(); ++ki)
    {
      shared_ptr<BoundedSurface> bdsf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(faces[ki]->surface());
      if (!bdsf.get())
	{
	  return; // Surface shouldnt have any inner loops. Nothing to remove
	}

      // Identify non-rotational trimming loops
      int nmbloops = faces[ki]->nmbBoundaryLoops();
      for (int ka=1; ka<nmbloops; )
	{
	  // Only inner trimming loops
	  shared_ptr<Loop> loop = faces[ki]->getBoundaryLoop(ka);
	  bool consistent = checkRotationalLoop(loop, centre, axis, vec);
	  if (consistent)
	    {
	      ++ka;
	      continue;   // Consistent with rotational model
	    }

	  // Disconnect faces along trimming curve
	  disconnectAlongLoop(loop);

	  // Store trimming loop
	  removed_loops_.push_back(loop);
	  
	  // Remove trimming loop from surface and face
	  // First identify boundary loop
	  vector<CurveLoop> bdloops = bdsf->allBoundaryLoops();
	  shared_ptr<ParamCurve> cv1 = 
	    loop->getEdge(0)->geomEdge()->geomCurve();
	  size_t kj;
	  for (kj=1; kj<bdloops.size(); ++kj)
	    {
	      // Only inner loops
	      int n_cvs = bdloops[kj].size();
	      int kb;
	      for (kb=0; kb<n_cvs; ++kb)
		{
		  shared_ptr<ParamCurve> cv2 = bdloops[kj][kb];
		  if (cv1.get() == cv2.get())
		    break;
		}
	      if (kb < n_cvs)
		break;
	    }
	  if (kj < bdloops.size())
	    bdsf->removeBdLoop(kj);
	  
	  // Remove face loop
	  faces[ki]->removeLoop(ka);
	  nmbloops--;
	}
    }
}

//==========================================================================
bool
HybridTrimVolume::checkRotationalLoop(shared_ptr<Loop>& loop, 
				      const Point& centre, const Point& axis, 
				      const Point& vec)
//==========================================================================
{
  // Fetch all assosiated curves
  vector<shared_ptr<ParamCurve> > trimcvs;
  int size = (int)loop->size();
  int ki, kj;
  for ( ki=0; ki<size; ++ki)
    {
      ftEdge* edg = loop->getEdge(ki)->geomEdge();
      if (!edg)
	return false;  // Something strange is going on
      trimcvs.push_back(edg->geomCurve());
    }

  // Check if the loop is smooth
  double angtol = model_->getTolerances().kink;
  double eps = model_->getTolerances().gap;
  int nmb = (int)trimcvs.size();
  for (kj=0, kj=1; ki<(int)trimcvs.size(); ++ki, ++kj)
    {
      kj = kj % nmb;
      vector<Point> pt1(2);
      vector<Point> pt2(2);
      trimcvs[ki]->point(pt1, trimcvs[ki]->endparam(), 1);
      trimcvs[kj]->point(pt2, trimcvs[kj]->startparam(), 1);
      if (pt1[1].angle(pt2[1]) > angtol)
	{
	  return false;
	}
    }

  // Check for consistency with given rotational axis
  vector<Point> plane_norm(2);
  for (ki=0; ki<nmb; ++ki)
    {
      Point centre2, axis2, vec2, dir;
      double angle2;
      bool rot = trimcvs[ki]->isAxisRotational(centre2, axis2, vec2, angle2);

      if (!rot)
	return false;

      double tmp_ang = axis.angle(axis2);
      if (tmp_ang > angtol && fabs(M_PI-tmp_ang) > angtol)
	return false;   // Different axis
      
      Point tmp_vec = centre2 - centre;
      tmp_ang = axis.angle(tmp_vec);
      if (tmp_vec.length() > eps && tmp_ang > angtol &&
	  fabs(M_PI-tmp_ang) > angtol)
	return false;   // Centres not rotationally consistent
    }
  return true;
}

//==========================================================================
bool
HybridTrimVolume::checkOuterLoop(shared_ptr<ftSurface>& face,
				 const Point& centre, const Point& axis, 
				 const Point& vec, 
				 shared_ptr<CurveLoop>& newloop,
				 vector<vector<shared_ptr<ftEdge> > >& removed_edgs)
//==========================================================================
{
  double eps = model_->getTolerances().gap;
  double angtol = model_->getTolerances().kink;
  shared_ptr<BoundedSurface> bdsf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(face->surface());
  if (!bdsf.get())
    {
      return true; // Surface has iso-parametric boundary curves
    }
  shared_ptr<Loop> loop = face->getBoundaryLoop(0);
  size_t nmb = loop->size();
  shared_ptr<ParamSurface> undersf = bdsf->underlyingSurface();
  Point normal;
  bool planar = undersf->isPlanar(normal, eps);
  if (planar)
    {
      // Identify trimming curves not compliant with the rotational
      // axis
      vector<shared_ptr<ftEdge> > nonrot_edgs;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ftEdge> edge = 
	    dynamic_pointer_cast<ftEdge,ftEdgeBase>(loop->getEdge(ki));
	  if (!edge.get())
	    continue;
	  shared_ptr<ParamCurve> cv = edge->geomCurve();

	  Point centre2, axis2, vec2, dir2;
	  double angle;
	  bool rotational = cv->isAxisRotational(centre2, axis2, vec2, angle);
	  bool linear = cv->isLinear(dir2, eps);
	  if (linear)
	    {
	      Point pos = cv->point(cv->startparam());
	      Point dir = pos - centre;
	      Point dir3 = dir - (dir*axis)*dir;
	      if (dir.angle(dir3) > angtol)
		nonrot_edgs.push_back(edge);
	    }
	  else if (rotational)
	    {
	      // Check consistence
	      Point dir = centre - centre2;
	      double ang1 = axis.angle(dir);
	      double ang2 = axis.angle(axis2);
	      if ((ang1 > angtol && fabs(M_PI-ang1) > angtol) ||  
		  (ang2 > angtol && fabs(M_PI-ang2) > angtol))
		nonrot_edgs.push_back(edge);
	    }
	  else
	    nonrot_edgs.push_back(edge);
	}
      if (nonrot_edgs.size() > 0)
	{
	  // Define clean rotational outer loop
#ifdef DEBUG
	  std::ofstream of("nonrot_edgs.g2");
	  for (size_t ki=0; ki<nonrot_edgs.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = nonrot_edgs[ki]->geomCurve();
	      shared_ptr<CurveOnSurface> bdcv = 
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (bdcv.get())
		{
		  bdcv->spaceCurve()->writeStandardHeader(of);
		  bdcv->spaceCurve()->write(of);
		}
	      else
		{
		  cv->writeStandardHeader(of);
		  cv->write(of);
		}
	    }
#endif
	  int stop_break = 1;
	}
    }
  else
    {
      // Identify loop edges that are not consistent with iso-parametric
      // trimming curves
      RectDomain dom = bdsf->containingDomain();
      shared_ptr<ftEdge> last;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ftEdge> edge = 
	    dynamic_pointer_cast<ftEdge,ftEdgeBase>(loop->getEdge(ki));
	  if (!edge.get())
	    continue;
	  shared_ptr<ParamCurve> cv = edge->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  if (!sfcv.get())
	    continue;

	  int dir;
	  double val;
	  bool constant = sfcv->isConstantCurve(eps, dir, val);
	  if (constant)
	    {
	      // Check if the curve corresponds to the expected parameter domain
	      if (dir == 1)
		{
		  if (fabs(val-dom.umin()) > eps && fabs(dom.umax()-val) > eps)
		    constant = false;
		}
	      else if (dir == 2)
		{
		  if (fabs(val-dom.vmin()) > eps && fabs(dom.vmax()-val) > eps)
		    constant = false;
		}
	    }
	  if (!constant)
	    {
	      // Store edge
	      if (removed_edgs.size() == 0 || !last.get() ||
		  last->next() != edge.get())
		{
		  vector<shared_ptr<ftEdge> > curr_edgs;
		  curr_edgs.push_back(edge);
		  removed_edgs.push_back(curr_edgs);
		}
	      else
		{
		  removed_edgs[removed_edgs.size()-1].push_back(edge);
		}
	      last = edge;
	    }
	}
  
      if (removed_edgs.size() > 1 && last.get())
	{
	  if (last->next() == removed_edgs[0][0].get())
	    {
	      size_t ix = removed_edgs.size()-1;
	      removed_edgs[0].insert(removed_edgs[0].begin(),
				     removed_edgs[ix].begin(), 
				     removed_edgs[ix].end());
	      removed_edgs.pop_back();
	    }
	}

      if (removed_edgs.size() > 1)
	{
	  // Create containing boundary loop
	  shared_ptr<ParamSurface> surf = bdsf->underlyingSurface();
	  vector<shared_ptr<ParamSurface> > sub_sfs = 
	    surf->subSurfaces(dom.umin(), dom.vmin(), dom.umax(), dom.vmax());
	  int pardir[] = {2, 1, 2, 1};
	  double parval[] = {dom.vmin(), dom.umax(), dom.vmax(), dom.umin()};
	  double extent[] = {dom.umin(), dom.umax(), dom.vmin(), dom.vmax(),
			     dom.umin(), dom.umax(), dom.vmin(), dom.vmax()};
	  vector< shared_ptr< ParamCurve > >  vec;
	  for (int ka=0; ka<4; ++ka)
	    {
	      vector<shared_ptr<ParamCurve> > constcvs =
		sub_sfs[0]->constParamCurves(parval[ka], pardir[ka] == 2);
	      shared_ptr<ParamCurve> spacecv = constcvs[0];
	      Point p1, p2;
	      if (pardir[ka] == 1)
		{
		  p1 = Point(parval[ka],extent[2*ka]);
		  p2 = Point(parval[ka],extent[2*ka+1]);
		}
	      else
		{
		  p1 = Point(extent[2*ka],parval[ka]);
		  p2 = Point(extent[2*ka+1],parval[ka]);
		}
	      shared_ptr<ParamCurve> pcv(new SplineCurve(p1, spacecv->startparam(),
							 p2, spacecv->endparam()));
	      shared_ptr<ParamCurve> sfcv(new CurveOnSurface(surf, pcv, spacecv,
							     false, 3));
	      if (ka == 2 || ka == 3)
		sfcv->reverseParameterDirection();
	      vec.push_back(sfcv);
	    }
	  newloop = shared_ptr<CurveLoop>(new CurveLoop(vec, eps));
	}
    }
  return (removed_edgs.size() == 0);
}

//==========================================================================
void
HybridTrimVolume::disconnectAlongLoop(shared_ptr<Loop>& loop)
//==========================================================================
{
  size_t ki, kj;
  vector<shared_ptr<ftEdgeBase> > edgs = loop->getEdges();
  for (ki=0; ki<edgs.size(); ++ki)
    {
      if (edgs[ki]->twin())
	edgs[ki]->disconnectTwin();

      // Disconnect radial edge
      shared_ptr<EdgeVertex> radial_edge = 
	edgs[ki]->geomEdge()->getEdgeMultiplicityInstance();
      if (radial_edge.get())
	{
	  radial_edge->removeEdge(edgs[ki]->geomEdge());
	  edgs[ki]->geomEdge()->removeEdgeVertex();
	}
    }

  // Disconnect vertices
  vector<shared_ptr<Vertex> > vx = loop->getVertices();
  const ftFaceBase* face = loop->getFace();
  for (ki=0; ki<vx.size(); ++ki)
    {
      vector<ftEdge*> curr_edges = vx[ki]->getFaceEdges((ftSurface*)face);
      for (kj=0; kj<curr_edges.size(); ++kj)
	vx[ki]->removeEdge(curr_edges[kj]);

      shared_ptr<Vertex> new_vx = 
	shared_ptr<Vertex>(new Vertex(vx[ki]->getVertexPoint(), curr_edges));
      for (kj=0; kj<curr_edges.size(); ++kj)
	curr_edges[kj]->replaceVertex(vx[ki], new_vx);
    }
}

struct projectinfo
{
  shared_ptr<ParamCurve> crv_;
  double startpar_, endpar_;
  shared_ptr<ParamSurface> srf_;
  Point sfpar1_, sfpar2_;

  projectinfo(shared_ptr<ParamCurve> crv, double startpar, double endpar,
	      shared_ptr<ParamSurface> srf, Point sfpar1, Point sfpar2)
  {
    crv_ = crv;
    startpar_ = startpar;
    endpar_ = endpar;
    srf_ = srf;
    sfpar1_ = sfpar1;
    sfpar2_ = sfpar2;
  }
};
  
//==========================================================================
void
HybridTrimVolume::fetchFreeLoops(shared_ptr<SurfaceModel>& model,
				 vector<vector<shared_ptr<ftEdge> > >& edgs,
				 vector<shared_ptr<ParamSurface> >& sfs_to_trim,
				 vector<vector<shared_ptr<CurveOnSurface> > >& trim_cvs)
//==========================================================================
{
  int nmb_bd = model->nmbBoundaries();
  double eps1 = model_->getTolerances().neighbour;
  double eps2 = model_->getTolerances().gap;
  double a_tol = 1.0e-10;
  vector<projectinfo> project;
  for (int ki=0; ki<nmb_bd; ++ki)
    {
      vector<shared_ptr<ftEdge> > bd_edgs =
	model->getBoundaryEdges(ki);
      edgs.push_back(bd_edgs);

      for (size_t kj=0; kj<bd_edgs.size(); ++kj)
	{
	  // Project curve onto rotational model boundary 
	  // Check position of endpoints
	  double p1 = bd_edgs[kj]->tMin();
	  double p2 = bd_edgs[kj]->tMax();
	  Point pos1 = bd_edgs[kj]->point(p1);
	  Point pos2 = bd_edgs[kj]->point(p2);

	  Point close1, close2;
	  int idx1, idx2;
	  double dist1, dist2;
	  double par1[2], par2[2];
	  rotational_shell_->closestPoint(pos1, close1, idx1, par1, dist1);
	  rotational_shell_->closestPoint(pos2, close2, idx2, par2, dist2);
	  shared_ptr<ParamCurve> space_cv = bd_edgs[kj]->geomCurve();
	  shared_ptr<ParamCurve> space_cv2(space_cv->subCurve(bd_edgs[kj]->tMin(),
							      bd_edgs[kj]->tMax()));
	  std::ofstream of1("proj_pts.g2");
	  shared_ptr<CurveOnSurface> sfcv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(space_cv2);
	  if (sfcv.get())
	    space_cv2 = sfcv->spaceCurve();
#ifdef DEBUG
	  if (space_cv2)
	    {
	      space_cv2->writeStandardHeader(of1);
	      space_cv2->write(of1);
	    }
	  of1 << "400 1 0 4 255 0 0 255" << std::endl;
	  of1 << "2" << std::endl;
	  of1 << close1 << std::endl;
	  of1 << close2 << std::endl;
#endif
	  // Check if the curve intersects the surface(s) boundary loop(s)
	  vector<int> sfix;
	  vector<pair<double,double> > par_lim;
	  vector<pair<Point,Point> > sfpar_lim;
	  sfix.push_back(idx1);
	  par_lim.push_back(make_pair(p1,p2));
	  sfpar_lim.push_back(make_pair(Point(par1[0],par1[1]),
					Point(par2[0],par2[1])));
	  // if (idx2 != idx1)
	  //   {
	  //     sfix.push_back(idx2);
	  //     par_lim.push_back(make_pair(p1,p2));
	  //     sfpar_lim.push_back(make_pair(Point(par1[0],par1[1]),
	  // 				    Point(par2[0],par2[1])));
	  //   }

	  for (size_t kr=0; kr<sfix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> curr_sf = 
		rotational_shell_->getSurface(sfix[kr]);
	      CurveLoop bdloop = curr_sf->outerBoundaryLoop();
	      vector<double> cvint;
	      vector<pair<int,double> > loopint;
	      bdloop.intersect(space_cv2, eps1, loopint, cvint);
	      if (cvint.size() == 0)
		{
		  // Curve segment inside surface. Store projection
		  // info
		  projectinfo curr_project(space_cv2, par_lim[kr].first, 
					   par_lim[kr].second, curr_sf,
					   sfpar_lim[kr].first,
					   sfpar_lim[kr].second);
		  project.push_back(curr_project);
		}
	      else
		{
		  for (size_t kw=0; kw<cvint.size(); ++kw)
		    {
		      // Compute distance
		      shared_ptr<ParamCurve> tmpcv = bdloop[loopint[kw].first];
		      Point intpos1 = space_cv2->point(cvint[kw]);
		      Point intpos2 = tmpcv->point(loopint[kw].second);
		      double dist1 = intpos1.dist(intpos2);
		      // Post iteration
		      double par1, par2, dist2;
		      Point ptc1, ptc2;
		      ClosestPoint::closestPtCurves(space_cv2.get(), 
						    tmpcv.get(), 
						    space_cv2->startparam(), 
						    space_cv2->endparam(),
						    tmpcv->startparam(), 
						    tmpcv->endparam(),
						    cvint[kw], 
						    loopint[kw].second,
						    par1, par2, dist2, 
						    ptc1, ptc2);
		      if (dist2 < dist1)
			{
			  cvint[kw] = par1;
			  loopint[kw].second = par2;
			}
		    }

		  if (par_lim[kr].first < cvint[0]-eps2)
		    {
		      cvint.insert(cvint.begin(), par_lim[kr].first);
		      loopint.insert(loopint.begin(), 
				     make_pair(-1,0.0)); // Placeholder
		    }
		  if (par_lim[kr].second > cvint[cvint.size()-1]+eps2)
		    {
		      cvint.push_back(par_lim[kr].second);
		      loopint.push_back(make_pair(-1,0.0)); // Placeholder
		    }

		  for (size_t kh=1; kh<cvint.size(); ++kh)
		    {
		      // Check if segment belongs to current surface
		      double midpar = 0.5*(cvint[kh-1] + cvint[kh]);
		      Point mid = space_cv2->point(midpar);
		      double upar, vpar, dd;
		      double dd2 = std::numeric_limits<double>::max(); 
		      Point closemid, closemid2;
		      int idxmid;
		      double parmid[2];
		      curr_sf->closestPoint(mid, upar, vpar, closemid, dd,
					    eps2);
		      if (dd > eps1)
			{
			  // Check also with the remainder of the model
			  rotational_shell_->closestPoint(mid, closemid2, 
							  idxmid, parmid, dd2);
			}
		      Point pos3 = space_cv2->point(cvint[kh-1]);
		      Point pos4 = space_cv2->point(cvint[kh]);
		      double upar3, vpar3, upar4, vpar4, dd3, dd4;
		      Point close3, close4;
		      shared_ptr<ParamSurface> curr_sf2 = (dd < dd2) ? 
			curr_sf : rotational_shell_->getSurface(idxmid);
			  
		      curr_sf2->closestPoint(pos3, upar3, vpar3, close3,
					     dd3, eps2);
		      curr_sf2->closestPoint(pos4, upar4, vpar4, close4,
					     dd4, eps2);
		      if (dd < dd2+a_tol || idxmid == sfix[kr])
			{
			  projectinfo curr_project(space_cv2, cvint[kh-1],
						   cvint[kh], curr_sf2,
						   Point(upar3,vpar3),
						   Point(upar4,vpar4));
			  project.push_back(curr_project);
			  if (cvint[kh-1] > par_lim[kr].first)
			    {
			      par_lim[kr].first = cvint[kh-1];
			      sfpar_lim[kr].first = Point(upar3,vpar3);
			    }
			  if (cvint[kh] < par_lim[kr].second)
			    {
			      par_lim[kr].second = cvint[kh];
			      sfpar_lim[kr].second = Point(upar4,vpar4);
			    }
			}
		      else
			{
			  // Check if the curve piece exists already
			  size_t kw;
			  for (kw=0; kw<sfix.size(); ++kw)
			    {
			      if (sfix[kw] == idxmid)
				{
				  // Check consinstence of endpoints
				  if (fabs(cvint[kh-1]-par_lim[kw].first) < eps2 &&
				      fabs(cvint[kh]-par_lim[kw].second) < eps2)
				    break;
				}
			    }
			  if (kw == sfix.size())
			    {
			      sfix.push_back(idxmid);
			      par_lim.push_back(make_pair(cvint[kh-1],cvint[kh]));
			      sfpar_lim.push_back(make_pair(Point(upar3,vpar3),
							    Point(upar4,vpar4)));
			    }
			}
		    }
		}
	      int stop_break = 1;
	    }

	}
      int stop_break2 = 1;
    }

#ifdef DEBUG
  std::ofstream ofproj("proj_cvs.g2");
#endif  
  for (size_t kj=0; kj<project.size(); ++kj)
    {
      if (project[kj].endpar_ - project[kj].startpar_ <= a_tol)
	continue;

      shared_ptr<ParamCurve> subcv
	(project[kj].crv_->subCurve(project[kj].startpar_, project[kj].endpar_));
      shared_ptr<Point> startp(new Point(project[kj].sfpar1_[0], 
					 project[kj].sfpar1_[1]));
      shared_ptr<Point> endp(new Point(project[kj].sfpar2_[0], 
				       project[kj].sfpar2_[1]));
      shared_ptr<ParamCurve> curr_proj_cv(CurveCreators::projectSpaceCurve
					   (subcv,
					    project[kj].srf_,
					    startp, endp, eps1 /*eps2*/));
      if (!curr_proj_cv.get())
	{
	  // We might have hit a surface seem
	  shared_ptr<ParamSurface> under_sf = project[kj].srf_;
	  shared_ptr<SurfaceOnVolume> vol_sf =
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(under_sf);
	  if (vol_sf.get())
	    under_sf = vol_sf->spaceSurface();
	  RectDomain dom = under_sf->containingDomain();
	  bool closed_u, closed_v;
	  SurfaceTools::checkSurfaceClosed(*under_sf, closed_u, closed_v, eps1);

	  vector<shared_ptr<Point> > tmp_start, tmp_end;
	  if (closed_u && (*startp)[0]-dom.umin() < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point(dom.umax(),(*startp)[1]));
	      tmp_start.push_back(tmp_pt);
	      tmp_end.push_back(endp);
	    }
	  if (closed_u && dom.umax()-(*startp)[0] < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point(dom.umin(),(*startp)[1]));
	      tmp_start.push_back(tmp_pt);
	      tmp_end.push_back(endp);
	    }
	  if (closed_u && (*endp)[0]-dom.umin() < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point(dom.umax(),(*endp)[1]));
	      tmp_start.push_back(startp);
	      tmp_end.push_back(tmp_pt);
	    }
	  if (closed_u && dom.umax()-(*endp)[0] < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point(dom.umin(),(*endp)[1]));
	      tmp_start.push_back(startp);
	      tmp_end.push_back(tmp_pt);
	    }

	  if (closed_v && (*startp)[1]-dom.vmin() < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point((*startp)[0],dom.vmax()));
	      tmp_start.push_back(tmp_pt);
	      tmp_end.push_back(endp);
	    }
	  if (closed_v && dom.vmax()-(*startp)[1] < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point((*startp)[0],dom.vmin()));
	      tmp_start.push_back(tmp_pt);
	      tmp_end.push_back(endp);
	    }
	  if (closed_v && (*endp)[1]-dom.vmin() < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point((*endp)[0],dom.vmax()));
	      tmp_start.push_back(startp);
	      tmp_end.push_back(tmp_pt);
	    }
	  if (closed_v && dom.vmax()-(*endp)[1] < eps2)
	    {
	      shared_ptr<Point> tmp_pt(new Point((*endp)[0],dom.vmin()));
	      tmp_start.push_back(startp);
	      tmp_end.push_back(tmp_pt);
	    }
	      
	  for (size_t kw=0; kw<tmp_start.size(); ++kw)
	    {
	      curr_proj_cv = shared_ptr<SplineCurve>(CurveCreators::projectSpaceCurve
						     (subcv,
						      project[kj].srf_,
						      tmp_start[kw], 
						      tmp_end[kw], eps1));
	      if (curr_proj_cv.get())
		break;
	    }
	}
      if (curr_proj_cv.get())
	{
	  shared_ptr<ParamSurface> under_sf = project[kj].srf_;
	  shared_ptr<CurveOnSurface> curr_bdcv(new CurveOnSurface(under_sf,
								  curr_proj_cv,
								  true));
	  curr_bdcv->ensureSpaceCrvExistence(eps2);
#ifdef DEBUG
	  if (curr_bdcv->hasSpaceCurve())
	    {
	      curr_bdcv->spaceCurve()->writeStandardHeader(ofproj);
	      curr_bdcv->spaceCurve()->write(ofproj);
	    }
#endif
	  size_t kr;
	  for (kr=0; kr<sfs_to_trim.size(); ++kr)
	    if (sfs_to_trim[kr].get() == under_sf.get())
	      break;
	  if (kr < sfs_to_trim.size())
	    trim_cvs[kr].push_back(curr_bdcv);
	  else
	    {
	      vector<shared_ptr<CurveOnSurface> > curr_trim;
	      curr_trim.push_back(curr_bdcv);
	      sfs_to_trim.push_back(under_sf);
	      trim_cvs.push_back(curr_trim);
	    }
	}	
    }

  int stop_break3 = 1;
}

//==========================================================================
bool
HybridTrimVolume::checkSubModelInternal(shared_ptr<SurfaceModel>& sub_model,
					shared_ptr<SurfaceModel>& rot_shell)
//==========================================================================
{
  // Fetch internal point in sub_model
  shared_ptr<ParamSurface> surf = sub_model->getSurface(0);
  double upar, vpar;
  Point inner_pt = surf->getInternalPoint(upar, vpar);

  // Inside test
  double dist;
  return rot_shell->isInside(inner_pt, dist);
}

//==========================================================================
void
HybridTrimVolume::performTrimming(vector<shared_ptr<ParamSurface> >& sfs_to_trim,
				  vector<vector<shared_ptr<CurveOnSurface> > >& trim_cvs,
				  vector<pair<shared_ptr<SurfaceModel>, bool > >& sub_models,
				  vector<vector<shared_ptr<ParamSurface> > >& associated_surface,
				  vector<vector<shared_ptr<ftSurface> > >& added_face,
				  vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
  double eps = model_->getTolerances().gap;
  double eps2 = model_->getTolerances().neighbour;

  associated_surface.resize(sub_models.size());
  added_face.resize(sub_models.size());

  // Fetch bounding boxes for each sub model
  vector<BoundingBox> bbox(sub_models.size());
  for (size_t kj=0; kj<sub_models.size(); ++kj)
    bbox[kj] = sub_models[kj].first->boundingBox();

  for (size_t ki=0; ki<sfs_to_trim.size(); ++ki)
    {
      // Make trimmed surfaces
      vector<shared_ptr<BoundedSurface> > trim_sfs;
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs_to_trim[ki]);
      if (!bd_sf.get())
	{
	  vector<CurveLoop> loops = 
	    SurfaceTools::absolutelyAllBoundarySfLoops(sfs_to_trim[ki], eps);

	  bd_sf = shared_ptr<BoundedSurface>(new BoundedSurface(sfs_to_trim[ki],
								loops));
	}
      try {
	trim_sfs = 
	  BoundedUtils::splitWithTrimSegments(bd_sf, trim_cvs[ki], eps);
      }
      catch(...)
	{
#ifdef DEBUG
	  std::cout << "Trimmed surfaces missing" << std::endl;
#endif
	}

      for (size_t kj=0; kj<trim_sfs.size(); ++kj)
	{
	  // First check if the trimmed surface belongs to the 
	  // initial model
	  // Fetch internal point
	  double uinner, vinner, distinner;
	  Point inner = trim_sfs[kj]->getInternalPoint(uinner, vinner);
	  bool inside = init_model_->isInside(inner, distinner);
	  if (distinner < eps)
	    {
	      // Boundary point of initial model
	      // Probably the rotational surface
	      shared_ptr<ftSurface> face(new ftSurface(trim_sfs[kj], -1));
	      trim_faces.push_back(face);
	      continue;  
	    }

	  // Identify associated sub model
	  // Fetch point on the outer boundary of the trimming surface
	  CurveLoop bdloop = trim_sfs[kj]->outerBoundaryLoop();
	  if (bdloop.size() == 0)
	    continue;  // Doesn't make sense
	  shared_ptr<ParamCurve> cv = bdloop[0];
	  Point pos = cv->point(0.5*(cv->startparam()+cv->endparam()));

	  // Identify closest sub model boundary loop
	  double mindist = std::numeric_limits<double>::max(); 
	  int minix = -1;
	  size_t kr, kh;
	  int ka;
	  for (kh=0; kh<sub_models.size(); ++kh)
	    {
	      // Check box
	      double boxdist = bbox[kh].dist(pos);
	      if (boxdist > mindist)
		continue;   // Not this sub model

	      // Check edges
	      int nmb_bd = sub_models[kh].first->nmbBoundaries();
	      for (ka=0; ka<nmb_bd; ++ka)
		{
		  vector<shared_ptr<ftEdge> > bd_edgs =
		    sub_models[kh].first->getBoundaryEdges(ka);
		  for (kr=0; kr<bd_edgs.size(); ++kr)
		    {
		      double tpar, edg_dist;
		      Point edg_pt;
		      bd_edgs[kr]->closestPoint(pos, tpar, edg_pt, edg_dist);
		      if (edg_dist < mindist)
			{
			  mindist = edg_dist;
			  minix = (int)kh;
			}
		    }
		}
	    }
	  associated_surface[minix].push_back(sfs_to_trim[ki]);
	  if (sub_models[minix].second)
	    {
	      // Sub model represents subtraction from the rotational model
	      // Do not store trim surface
	      ;
	    }
	  else
	    {
	      // shared_ptr<BoundedSurface> trim2(trim_sfs[kj]->clone());
	      // trim2->swapParameterDirection();
	      // shared_ptr<SurfaceOnVolume> volsf =
	      // 	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(trim2->underlyingSurface());
	      // if (volsf.get())
	      // 	volsf->unsetParamSurf();  // Not associated the roational volume

	      // // Create face and associate it to the identified sub model
	      // shared_ptr<ftSurface> face2(new ftSurface(trim2, -1));
	      // //sub_models[minix].first->append(face2, true);
	      // added_face[minix].push_back(face2);

	      // Create and store trimming face for the rotational model
	      shared_ptr<ftSurface> face(new ftSurface(trim_sfs[kj], -1));
	      trim_faces.push_back(face);
	      shared_ptr<ftSurface> face2(new ftSurface(trim_sfs[kj], -1));
	      added_face[minix].push_back(face2);
	    }
	}
#ifdef DEBUG
      std::ofstream offace("trim_faces.g2");
      for (size_t kj=0; kj<trim_faces.size(); ++kj)
	{
	  shared_ptr<ParamSurface> tmp_sf = trim_faces[kj]->surface();
	  tmp_sf->writeStandardHeader(offace);
	  tmp_sf->write(offace);
	}
      std::ofstream offace2("added_faces.g2");
      for (size_t kj=0; kj<added_face.size(); ++kj)
	{
	  for (size_t kr=0; kr<added_face[kj].size(); ++kr)
	    {
	      shared_ptr<ParamSurface> tmp_sf = added_face[kj][kr]->surface();
	      tmp_sf->writeStandardHeader(offace2);
	      tmp_sf->write(offace2);
	    }
	}
#endif
      int stop_break = 1;
    }
}


//==========================================================================
void
HybridTrimVolume::insertTrimModel(shared_ptr<SurfaceModel>& trim_model,
				  vector<shared_ptr<ftSurface> >& bd_faces)
//==========================================================================
{
  if (bd_faces.size() == 0)
    return;

#ifdef DEBUG
  std::ofstream trim("trim_mod.g2");
  int nmb = trim_model->nmbEntities();
  for (int ka=0; ka<nmb; ++ka)
    {
      shared_ptr<ParamSurface> tmp_sf = 
	trim_model->getSurface(ka);
      tmp_sf->writeStandardHeader(trim);
      tmp_sf->write(trim);
    }
#endif

  // Fetch block associated to first boundary face
  shared_ptr<ParamSurface> bd_sf = bd_faces[0]->surface();
  shared_ptr<BoundedSurface> bd2_sf = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(bd_sf);
  if (bd2_sf.get())
    bd_sf = bd2_sf->underlyingSurface();

  vector<Body*> cand_bodies;
  Body* body = bd_faces[0]->getBody();
  if (body)
    cand_bodies.push_back(body);

  vector<shared_ptr<ftVolume> > trim_volumes;
  for (size_t ki=0; ki<cand_bodies.size(); ++ki)
    {
      vector<Body*> tmp_neighbours;
      cand_bodies[ki]->getAdjacentBodies(tmp_neighbours);

      shared_ptr<SurfaceModel> shell = cand_bodies[ki]->getOuterShell();
      if (shell.get())
	{
#ifdef DEBUG
	  std::ofstream shout("block_mod.g2");
	  std::ofstream ofbd0("bd_block.g2");
	  int nmb = shell->nmbEntities();
	  for (int ka=0; ka<nmb; ++ka)
	    {
	      shared_ptr<ParamSurface> tmp_sf = 
		shell->getSurface(ka);
	      tmp_sf->writeStandardHeader(shout);
	      tmp_sf->write(shout);
	      shared_ptr<ftSurface> tmp_face = shell->getFace(ka);
	      if (!tmp_face->twin())
		{
		  tmp_sf->writeStandardHeader(ofbd0);
		  tmp_sf->write(ofbd0);
		}
	    }
#endif
	  // // Create reduced surface model
	  // vector<shared_ptr<ftSurface> > faces2;
	  // int nmb_faces = shell->nmbEntities();
	  // tpTolerances tptol = shell->getTolerances();
	  // for (int ka=0; ka<nmb_faces; ++ka)
	  //   {
	  //     shared_ptr<ftSurface> shellface = shell->getFace(ka);
	  //     if (shellface.get() == bd_faces[0].get())
	  // 	continue;
	  //     faces2.push_back(shared_ptr<ftSurface>(new ftSurface(shellface->surface(), 
	  // 							   ka)));
	  //   }
	  // shared_ptr<SurfaceModel> model2(new SurfaceModel(tptol.gap, tptol.gap,
	  // 						   tptol.neighbour,
	  // 						   tptol.kink, tptol.bend,
	  // 						   faces2));
	  // vector<shared_ptr<SurfaceModel> > split_models =
	  //   model2->splitSurfaceModels(trim_model);

	  // Remove faces touching the trim model (split already performed)
	  vector<shared_ptr<ftSurface> > save_face;
	  int nmb_faces = shell->nmbEntities();
	  for (int ka=0; ka<nmb_faces; ++ka)
	    {
	      shared_ptr<ftSurface> shellface = shell->getFace(ka);
	      for (size_t kj=0; kj<bd_faces.size(); ++kj)
		if (shellface.get() == bd_faces[kj].get())
		  {
		    save_face.push_back(shellface);
		    break;
		  }
	    }
	  for (size_t kj=0; kj<save_face.size(); ++kj)
	    bool removed = shell->removeFace(save_face[kj]);
		  
	  vector<shared_ptr<SurfaceModel> > split_models =
	    shell->splitSurfaceModels(trim_model);

	  // Reenter removed faces
	  for (size_t kj=0; kj<save_face.size(); ++kj)
	    {
	      int nmb = shell->nmbEntities();
	      shared_ptr<ftSurface> face2(new ftSurface(save_face[kj]->surface(), 
							nmb));
	      shell->append(face2);
	    }
	      

#ifdef DEBUG
	  if (split_models[0].get())
	    {
	      std::ofstream split1("split_mod1.g2");
	      int nmb1 = split_models[0]->nmbEntities();
	      for (int ka=0; ka<nmb1; ++ka)
		{
		  shared_ptr<ParamSurface> tmp_sf = 
		    split_models[0]->getSurface(ka);
		  tmp_sf->writeStandardHeader(split1);
		  tmp_sf->write(split1);
		}
	    }
	  if (split_models[1].get())
	    {
	      std::ofstream split2("split_mod2.g2");
	      int nmb2 = split_models[1]->nmbEntities();
	      for (int ka=0; ka<nmb2; ++ka)
		{
		  shared_ptr<ParamSurface> tmp_sf = 
		    split_models[1]->getSurface(ka);
		  tmp_sf->writeStandardHeader(split2);
		  tmp_sf->write(split2);
		}
	    }
	  if (split_models[2].get())
	    {
	      std::ofstream split3("split_mod3.g2");
	      int nmb3 = split_models[2]->nmbEntities();
	      for (int ka=0; ka<nmb3; ++ka)
		{
		  shared_ptr<ParamSurface> tmp_sf = 
		    split_models[2]->getSurface(ka);
		  tmp_sf->writeStandardHeader(split3);
		  tmp_sf->write(split3);
		}
	    }
	  if (split_models[3].get())
	    {
	      std::ofstream split4("split_mod4.g2");
	      int nmb4 = split_models[3]->nmbEntities();
	      for (int ka=0; ka<nmb4; ++ka)
		{
		  shared_ptr<ParamSurface> tmp_sf = 
		    split_models[3]->getSurface(ka);
		  tmp_sf->writeStandardHeader(split4);
		  tmp_sf->write(split4);
		}
	    }
#endif
	  if (split_models[0].get() && split_models[1].get())
	    {
	      // A trimming is performed. Collect faces for trimmed block
	      if (save_face.size() > 0)
		{
		  split_models[0]->append(save_face);
		  // for (size_t kj=0; kj<save_face.size(); ++kj)
		  // 	{
		  // 	  int nmb = split_models[0]->nmbEntities();
		  // 	  shared_ptr<ftSurface> face2(new ftSurface(save_face[kj]->surface(), 
		  // 						    nmb));
		  // 	  split_models[0]->append(face2);
		  // 	}
		}
	      if (split_models[2]->nmbEntities() > 0)
		split_models[0]->append(split_models[2]);
#ifdef DEBUG
	      std::ofstream trim("trim_block.g2");
	      nmb = split_models[0]->nmbEntities();
	      for (int ka=0; ka<nmb; ++ka)
		{
		  shared_ptr<ParamSurface> tmp_sf = 
		    split_models[0]->getSurface(ka);
		  tmp_sf->writeStandardHeader(trim);
		  tmp_sf->write(trim);
		}
#endif
	      // Update twin pointers. To be done

	      // Fetch adjacent blocks to check if these needs trimming
	      vector<Body*> neighbours;
	      cand_bodies[ki]->getAdjacentBodies(neighbours);
	      for (size_t kj=0; kj<neighbours.size(); ++kj)
		{
		  size_t kr;
		  for (kr=0; kr<cand_bodies.size(); ++kr)
		    if (cand_bodies[kr] == neighbours[kj])
		      break;  // Check blocks only once
		  if (kr == cand_bodies.size())
		    cand_bodies.push_back(neighbours[kj]);
		}

	      // Create topological volume
	      shared_ptr<ftVolume> curr_vol = 
		rotational_vol_->fetchAsSharedPtr(cand_bodies[ki]);
	      if (curr_vol.get())
		{
		  shared_ptr<ParamVolume> geom_vol =
		    curr_vol->getVolume();
		  shared_ptr<ftVolume> trimmed_vol(new ftVolume(geom_vol,
								split_models[0]));
		  trim_volumes.push_back(trimmed_vol);

		  // Remove current volume block from rotational model
		  rotational_vol_->removeSolid(curr_vol, false);
		}
	    }
	  int stop_break = 1;
	}
    }

  if (trim_volumes.size() > 0)
    {
      // Add trimmed volumes to rotational model
      rotational_vol_->append(trim_volumes);
    }
#ifdef DEBUG
  rotational_vol_->checkModelTopology();
  std::ofstream oftrimall("trimvolmod.g2");
  int nmb_vol = rotational_vol_->nmbEntities();
  for (int ka=0; ka<nmb_vol; ++ka)
    {
      shared_ptr<ftVolume> ftvol = rotational_vol_->getBody(ka);
      int nmb_shell = ftvol->nmbOfShells();
      for (int kb=0; kb<nmb_shell; ++kb)
	{
	  shared_ptr<SurfaceModel> sfmod = ftvol->getShell(kb);
	  int nmb_face = sfmod->nmbEntities();
	  for (int kc=0; kc<nmb_face; ++kc)
	    {
	      shared_ptr<ParamSurface> rotsf = sfmod->getSurface(kc);
	      rotsf->writeStandardHeader(oftrimall);
	      rotsf->write(oftrimall);
	    }
	}
    }

  std::ofstream ofbd("bd_faces.g2");
  vector<shared_ptr<ftSurface> > bd_faces2 = 
    rotational_vol_->getBoundaryFaces();
  for (size_t ki=0; ki<bd_faces2.size(); ++ki)
    {
      shared_ptr<ParamSurface> rotsf = bd_faces2[ki]->surface();
      rotsf->writeStandardHeader(ofbd);
      rotsf->write(ofbd);
    }

  std::ofstream ofbd2("bd_blocks2.g2");
  for (size_t ki=0; ki<trim_volumes.size(); ++ki)
    {
      shared_ptr<SurfaceModel> sfmod2 = trim_volumes[ki]->getShell(0);
      int nmb_face2 = sfmod2->nmbEntities();
      for (int kc=0; kc<nmb_face2; ++kc)
	    {
	      shared_ptr<ftSurface> rotf = sfmod2->getFace(kc);
	      if (!rotf->twin())
		{
		  shared_ptr<ParamSurface> rotsf = sfmod2->getSurface(kc);
		  rotsf->writeStandardHeader(ofbd2);
		  rotsf->write(ofbd2);
		}
	    }
	}
      
#endif
  int stop_finish = 1;
}

//==========================================================================
shared_ptr<VolumeModel>
HybridTrimVolume::subModel2Trivariate(shared_ptr<SurfaceModel>& sub_model,
				      vector<shared_ptr<ftSurface> >& added_face)
//==========================================================================
{
  tpTolerances tptol = sub_model->getTolerances();
  int degree = 3;
  if (added_face.size() == 1)
    {
      shared_ptr<ParamSurface> addsf(added_face[0]->surface()->clone());
      addsf->swapParameterDirection();
      shared_ptr<BoundedSurface> bdsf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(addsf);
      shared_ptr<ParamSurface> under = (bdsf.get()) ? 
	bdsf->underlyingSurface() : addsf;
      shared_ptr<SurfaceOnVolume> volsf = 
	dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(under);
      if (volsf.get())
	volsf->unsetParamSurf();
      shared_ptr<ftSurface> newface(new ftSurface(addsf, 
						  sub_model->nmbEntities()));
      
      sub_model->append(newface);
    }
  else if (added_face.size() > 1)
    {
#ifdef DEBUG
      std::ofstream of1("before_merge.g2");
      for (size_t ki=0; ki<added_face.size(); ++ki)
	{
	  shared_ptr<ParamSurface> tmpsf = added_face[ki]->surface();
	  tmpsf->writeStandardHeader(of1);
	  tmpsf->write(of1);
	}
#endif
      // Try to merge trimming faces to create a better starting
      // point for block structurering
      shared_ptr<SurfaceModel> addmod(new SurfaceModel(tptol.gap, tptol.gap,
						       tptol.neighbour,
						       tptol.kink, tptol.bend,
						       added_face));
      RegularizeFaceSet regularize(addmod, false, 0);
      regularize.removeExtraDiv(true);

#ifdef DEBUG
	std::ofstream of2("after_merge.g2");
	int tmpnmb = addmod->nmbEntities();
	for (int ka=0; ka<tmpnmb; ++ka)
	  {
	  shared_ptr<ParamSurface> tmpsf = addmod->getSurface(ka);
	  tmpsf->writeStandardHeader(of2);
	  tmpsf->write(of2);
	}
#endif
	int nmbface = addmod->nmbEntities();
	for (int ka=0; ka<nmbface; ++ka)
	  {
	    shared_ptr<ParamSurface> tmpsf = addmod->getSurface(ka);
	    shared_ptr<ParamSurface> addsf(tmpsf->clone());
	    addsf->swapParameterDirection();
	    shared_ptr<BoundedSurface> bdsf = 
	      dynamic_pointer_cast<BoundedSurface, ParamSurface>(addsf);
	    shared_ptr<ParamSurface> under = (bdsf.get()) ? 
	      bdsf->underlyingSurface() : addsf;
	    shared_ptr<SurfaceOnVolume> volsf = 
	      dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(under);
	    if (volsf.get())
	      volsf->unsetParamSurf();
	    shared_ptr<ftSurface> newface(new ftSurface(addsf, 
							sub_model->nmbEntities()));
      
	    sub_model->append(newface);
	  }
	int stop_break = 1;
    }

  // Check if the model can be represented by one block
  shared_ptr<ftVolume> ftvol(new ftVolume(sub_model));
  bool reg = ftvol->isRegularized();
  shared_ptr<VolumeModel> volmod;
  vector<shared_ptr<ftVolume> > reg_vols;
  if (!reg)
    {
      vector<SurfaceModel*> modified_adjacent;
      try {
	reg_vols = 
	  ftvol->replaceWithRegVolumes(degree, modified_adjacent,
				       false, 1, false, true);
      }
      catch (...)
	{
	  ;
	}
    }
  if (reg_vols.size() == 0)
    reg_vols.push_back(ftvol);

  volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, 
						   tptol.gap,
						   tptol.neighbour,
						   tptol.kink, 
						   tptol.bend));

#ifdef DEBUG
  std::ofstream of("submod_vol.g2");
#endif
  int nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
      bool bd_trim = curr_vol->isBoundaryTrimmed();
      bool iso_trim = curr_vol->isIsoTrimmed();
      bool reg = curr_vol->isRegularized(true);

      if (reg)
	{
	  curr_vol->untrimRegular(degree, true);
#ifdef DEBUG
	  shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	  curr_vol2->writeStandardHeader(of);
	  curr_vol2->write(of);
#endif
	}
    }

  volmod->makeCommonSplineSpaces();
  volmod->averageCorrespondingCoefs();

  if (added_face.size() > 1)
    {
      // Reinsert added faces in the outer shell
    }

  return volmod;
}

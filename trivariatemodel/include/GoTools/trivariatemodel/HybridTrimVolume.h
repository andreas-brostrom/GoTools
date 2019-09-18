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

#ifndef _HYBRIDTRIMVOLUME_H
#define _HYBRIDTRIMVOLUME_H

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class ftSurface;
  class Loop;
  class ftVolume;
  class ParamVolume;
  class ParamSurface;
  class VolumeModel;

  enum sf_type
  {
    UNKNOWN = -1,
    FREEFORM = 0,
    PLANAR = 1,
    ROTATIONAL = 2
  };

  /// Create a multi-block volume model with trimmed blocks from a boundary 
  /// represented solid respecting characteristica of the input shape

  class HybridTrimVolume
  {
  public:
    /// Constructor
    /// model : The trimming shell / outer shell of the brep solid
    /// material: Eventual material specification associated with the brep solid
    ///           that should be maintained in the trivariate model
    HybridTrimVolume(shared_ptr<SurfaceModel> model, int material=-1);

    /// Destructor
    ~HybridTrimVolume();

    shared_ptr<VolumeModel> 
      fetchVolModel(bool create_degen = true, bool refine_sharp = false);

  private:
    shared_ptr<SurfaceModel> model_;
    int material_;
    shared_ptr<SurfaceModel> init_model_;
    BoundingBox bigbox_;
    shared_ptr<VolumeModel> rotational_vol_;
    shared_ptr<SurfaceModel> rotational_shell_;

    std::vector<shared_ptr<Loop> > removed_loops_;
    std::vector<std::vector<shared_ptr<ftEdge> > > edge_seqs_;

    bool 
      identifyRotationalAxis(Point& centre, Point& axis, 
			     Point& vec, double& angle,
			     std::vector<shared_ptr<ftSurface> >& rotational_faces, 
			     std::vector<shared_ptr<ftSurface> >& rotational_faces2, 
			     std::vector<shared_ptr<ftSurface> >& other_faces);

    void removeInnerTrim(std::vector<shared_ptr<ftSurface> >& faces, 
			 const Point& centre, const Point& axis, 
			 const Point& vec);

    void cleanOuterTrim(std::vector<shared_ptr<ftSurface> >& faces, 
			 const Point& centre, const Point& axis, 
			 const Point& vec);

    bool checkRotationalLoop(shared_ptr<Loop>& loop, 
			     const Point& centre, const Point& axis, 
			     const Point& vec);

    bool checkOuterLoop(shared_ptr<ftSurface>& face,
			const Point& centre, const Point& axis, 
			const Point& vec, 
			shared_ptr<CurveLoop>& newloop,
			std::vector<std::vector<shared_ptr<ftEdge> > >& removed_edgs);

    void disconnectAlongLoop(shared_ptr<Loop>& loop);

    void fetchFreeLoops(shared_ptr<SurfaceModel>& model,
			std::vector<std::vector<shared_ptr<ftEdge> > >& edgs,
			std::vector<shared_ptr<ParamSurface> >& sfs_to_trim,
			std::vector<std::vector<shared_ptr<CurveOnSurface> > >& trim_cvs);

  bool
    checkSubModelInternal(shared_ptr<SurfaceModel>& sub_model,
			  shared_ptr<SurfaceModel>& rot_shell);

  void
    performTrimming(std::vector<shared_ptr<ParamSurface> >& sfs_to_trim,
		    std::vector<std::vector<shared_ptr<CurveOnSurface> > >& trim_cvs,
		    std::vector<std::pair<shared_ptr<SurfaceModel>, bool > >& sub_models,
		    std::vector<std::vector<shared_ptr<ParamSurface> > >& associated_surface,
		    std::vector<std::vector<shared_ptr<ftSurface> > >& added_face,
		    std::vector<shared_ptr<ftSurface> >& trim_faces);

  void
    insertTrimModel(shared_ptr<SurfaceModel>& trim_model,
		    std::vector<shared_ptr<ftSurface> >& bd_faces);

  shared_ptr<VolumeModel>
    subModel2Trivariate(shared_ptr<SurfaceModel>& sub_model,
			std::vector<shared_ptr<ftSurface> >& added_face);

  };

} // namespace Go


#endif // _HYBRIDTRIMVOLUME_H

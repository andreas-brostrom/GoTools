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

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/FaceUtilities.h"
#include "GoTools/geometry/PointCloud.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 5) {
    std::cout << "Input parameters : Input file, output file, output file 2, density"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream fileout(argv[2]);
  std::ofstream fileout2(argv[3]);
  
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;
  double density = atof(argv[4]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);
  if (!sfmodel)
    exit(-1);

  std::vector<SamplePointData2> samples;
  sfmodel->fetchSamplePoints2(density, samples);

  vector<double> points;
  points.reserve(samples.size());
  for (size_t ki=0; ki<samples.size(); ++ki)
    {
      points.insert(points.end(), samples[ki].pos_.begin(), samples[ki].pos_.end());

      vector<Point> der(3);
      samples[ki].face_->surface()->point(der, samples[ki].face_par_[0],
					  samples[ki].face_par_[1], 1);
      Point max_curv = samples[ki].kvec1_[0]*der[1]+samples[ki].kvec1_[1]*der[2];
      (void)max_curv.normalize_checked();
      Point min_curv = samples[ki].kvec2_[0]*der[1]+samples[ki].kvec2_[1]*der[2];
      (void)min_curv.normalize_checked();
      fileout2 << samples[ki].pos_ << " " << samples[ki].norm_ << " ";
      fileout2 << max_curv << " " << samples[ki].k1_ << " ";
      fileout2 << min_curv << " " << samples[ki].k2_;
      fileout2 << " " << samples[ki].type_ << std::endl;
    }
  PointCloud3D pt_cloud(&points[0], (int)points.size()/3);
  pt_cloud.writeStandardHeader(fileout);
  pt_cloud.write(fileout);
  


  int break_point;
  break_point = 1;
}


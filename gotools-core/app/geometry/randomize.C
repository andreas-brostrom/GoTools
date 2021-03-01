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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/BoundingBox.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
  if (argc != 4)
    {
    std::cout << "Input parameters : Input file, output file, factor"  << std::endl;
    exit(-1);
  }
  
  // Read the point cloud from file
  ifstream input(argv[1]);
  ofstream output(argv[2]);
  double fac = atof(argv[3]);
  
  ObjectHeader header;
  PointCloud3D cloud;
  input >> header >> cloud;

  vector<Point> dir(10);
  dir[0] = Point(1.0, 0.0, 0.0);
  dir[1] = Point(0.0, 1.0, 0.0);
  dir[2] = Point(0.0, 0.0, 1.0);
  dir[3] = Point(1.0, 1.0, 0.0);
  dir[4] = Point(1.0, 0.0, 1.0);
  dir[5] = Point(0.0, 1.0, 1.0);
  dir[6] = Point(1.0, 1.0, 1.0);
  dir[7] = Point(1.0, 0.5, 0.5);
  dir[8] = Point(0.5, 1.0, 0.5);
  dir[9] = Point(0.5, 0.5, 1.0);
  
  BoundingBox box = cloud.boundingBox();
  double mlen = fac*(box.low().dist(box.high()));
  int r1;

  int num = cloud.numPoints();
  double *points = cloud.rawData();
  int dim = 3;
  int ki;
  double *curr;
  vector<double> rand_pts;
  for (ki=0, curr=points; ki<num; ++ki, curr+=dim)
    {
      r1 = std::rand();
      double f1 = ((double)r1/RAND_MAX)*mlen;
      int sgn = (r1 % 2 == 0) ? 1 : -1;
      int div = (int)(((double)r1/RAND_MAX + 0.00001)*10.0);
      div = std::min(div, 9);
      Point pos(curr, curr+dim);
      Point pos2 = pos +sgn*f1*dir[div];
      rand_pts.insert(rand_pts.end(), pos2.begin(), pos2.end());
    }

  PointCloud3D cloud2(&rand_pts[0], num);
  cloud2.writeStandardHeader(output);
  cloud2.write(output);
}

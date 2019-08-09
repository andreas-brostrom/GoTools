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

#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include <fstream>
#include <algorithm>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{
  if (argc != 5)
    {
      std::cout << "Input parameters: infile, level value, tolerance, outfile" << std:: endl;
      exit(1);
    }

    // Read the surfaces from file
    std::ifstream input(argv[1]);
    if (input.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    
    double val = atof(argv[2]); 

    double tol = atof(argv[3]);

    std::ofstream out(argv[4]);

   GoTools::init();

   shared_ptr<ObjectHeader> header(new ObjectHeader());

    header->read(input);
    shared_ptr<GeomObject> geom_obj(Factory::createObject(header->classType()));
    shared_ptr<SplineSurface> surf =
      dynamic_pointer_cast<SplineSurface,GeomObject>(geom_obj);
    if (!surf.get())
      {
	std::cout << "Object one is not a spline surface" << std::endl;
	exit(1);
      }
    surf->read(input);
    SISLSurf* ss_sisl = GoSurf2SISL(*surf, false);

    // Finding intersections, using SISL function sh1761.
    int jpt = 0; // number of single intersection point
    int jcrv = 0; // number of intersection curves
    int jsurf = 0; 
    double* gpar = SISL_NULL; // parameter values of the single intersection points
    double* spar = SISL_NULL; // dummy array
    int *pretop=SISL_NULL;
    SISLIntcurve** wcurve; // array containing descriptions of the intersection curves
    SISLIntsurf** wsurf; 

    auto qp  = newPoint(&val, 1, 1);
    auto qo1 = newObject(SISLSURFACE);
    qo1->s1 = ss_sisl;
    qo1->o1 = qo1;
    auto qo2 = newObject(SISLPOINT);
    qo2->p1 = qp;
  
    SISLIntdat* qintdat = SISL_NULL; // intersection result

    double epsge = 1e-6;
    int kstat = 0;
  
    // find intersections
    sh1761(qo1, qo2, epsge, &qintdat, &kstat);
    if (kstat < 0) 
      {
	std::cout << "Level value intersection failed" << std::endl;
	return 1;
      }

    // express intersections on output format (SISL)
    const int kdeg = 1;
    if (qintdat)
      hp_s1880 (qo1, qo2, kdeg, 2, 0, qintdat, &jpt, &gpar, &spar, &pretop,
		&jcrv, &wcurve, &jsurf, &wsurf, &kstat);
  
    if (gpar)                         free(gpar);
    if (spar)                         free(spar);
    if (pretop)                       free(pretop);
    if (qintdat)                      freeIntdat(qintdat);
    if ((bool)wcurve && (jcrv > 0)) freeIntcrvlist(wcurve, jcrv);
    for (int i = 0; i < jsurf; ++i)   freeIntsurf(wsurf[i]);
    if ((bool)wsurf && (jsurf > 0))          free(wsurf);
    
}

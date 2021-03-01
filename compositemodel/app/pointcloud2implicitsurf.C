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

#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "sisl.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
  if (argc != 5)
    {
    std::cout << "Input parameters : Input file, output file, output file viz, degree"  << std::endl;
    exit(-1);
  }
  
    // Read the point cloud from file
    ifstream input(argv[1]);
    ofstream output(argv[2]);
    ofstream outviz(argv[3]);
    int degree = atoi(argv[4]);
    ObjectHeader header;
    PointCloud3D cloud;
    input >> header >> cloud;

    // Implicitize
    ImplicitizePointCloudAlgo implicitize(cloud, degree);
    implicitize.perform();

    // Get result
    BernsteinTetrahedralPoly implicit;
    BaryCoordSystem3D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Check accuracy
    double avdist = 0.0;
    double maxdist = 0.0;
    int dim = 3;
    int numpt = cloud.numPoints();
    for (int ki=0; ki<numpt; ++ki)
      {
	Vector3D curr = cloud.point(ki);
	Vector4D bary = bc.cartToBary(curr);
	double dist = implicit(bary);
	maxdist = std::max(maxdist, fabs(dist));
	avdist += dist;
      }
    avdist /= (double)numpt;
    std::cout << "Maximum distance: " << maxdist << std::endl;
    std::cout << "Average distance: " << avdist << std::endl;

    // Write out implicit function
    std::cout << "Sigma_min: " << sigma_min << std::endl;
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to output" << endl;

    // Differentiate
    Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
    Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
    Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
    Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
    BernsteinTetrahedralPoly deriv1, deriv2, deriv3, deriv4;
    BernsteinTetrahedralPoly deriv2_1, deriv2_2, deriv2_3, deriv2_4;
    implicit.deriv(1, bdir1, deriv1);
    implicit.deriv(1, bdir2, deriv2);
    implicit.deriv(1, bdir3, deriv3);
    implicit.deriv(1, bdir4, deriv4);
    implicit.deriv(2, bdir1, deriv2_1);
    implicit.deriv(2, bdir2, deriv2_2);
    implicit.deriv(2, bdir3, deriv2_3);
    implicit.deriv(2, bdir4, deriv2_4);
    
    // Fetch points on the implicit surface, fetch also points where the gradient vanishes
    BoundingBox ptbox = cloud.boundingBox();
    Point low = ptbox.low();
    Point high = ptbox.high();

    Point dir(3);
    int nmb_sample;
    std::cout << "Box min: " << low << ", box max: " << high << std::endl;
    std::cout << "Give view direction: " << std::endl;
    std::cin >> dir;
    std::cout << "Give number of points: " << std::endl;
    std::cin >> nmb_sample;

    // Find extension of view area
    std::ofstream tmp("tmplines.g2");
    dir.normalize();
    Point bmid = 0.5*(low + high);

    tmp << "400 1 0 4 255 0 0 255" << std::endl;
    tmp << "1" << std::endl;
    tmp << bmid << std::endl;

    // Define bounding box as a surface model
    double gap = 1.0e-6;
    Point xdir(1.0, 0.0, 0.0);
    Point ydir(0.0, 1.0, 0.0);
    Point zdir(0.0, 0.0, 1.0);
    CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
    shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low, xdir, ydir, high[0]-low[0],
							  high[1]-low[1], high[2]-low[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-0.5*len*dir2, 0.0, bmid+0.5*len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-0.5*len*dir3, 0.0, bmid+0.5*len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;

    cv1->writeStandardHeader(tmp);
    cv1->write(tmp);
    cv2->writeStandardHeader(tmp);
    cv2->write(tmp);
    ssf->writeStandardHeader(tmp);
    ssf->write(tmp);

    int ik = degree + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    vector<double> points;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);
	    tmp << "400 1 0 4 0 255 0 255 " << std::endl;
	    tmp << "1" << std::endl;
	    tmp << sfpos << std::endl;
	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;
	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

	    Vector4D bary1 = bc.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc.baryToCart(bary1);
	    Vector3D tp2 = bc.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = implicit.pickLine(bary1, bary2);

	    // Compute zeroes of bernstein polynomial
	    // First make sisl curve
	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;

	    // Intersect
	    double eps = 1.0e-6;
	    int kstat = 0;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);

	    // Compute cartesian points and curves associated with intersections
	    for (kr=0; kr<kpt; ++kr)
	      {
		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		Vector3D pos = bc.baryToCart(barypt);
		points.insert(points.end(), pos.begin(), pos.end());

		// Check
		Vector4D bary = bc.cartToBary(pos);
		double dist = implicit(bary);
		if (dist > 1.0e-4)
		  std::cout << "dist: " << dist << ", ki= " << ki << ", kj= " << kj << std::endl;
		double dist2 = implicit(barypt);
		if (dist2 > 1.0e-4)
		  std::cout << "dist: " << dist2 << ", ki= " << ki << ", kj= " << kj << std::endl;
	      }
	    for (kr=0; kr<kcrv; ++kr)
	      {
		int ipt = intcv[kr]->ipoint;
		double par1 = intcv[kr]->epar1[0];
		double par2 = intcv[kr]->epar1[ipt-1];
		Point pp1 = (1.0-par1)*cart1 + par1*cart2;
		Point pp2 = (1.0-par2)*cart1 + par2*cart2;
		linesegs.insert(linesegs.end(), pp1.begin(), pp1.end());
		linesegs.insert(linesegs.end(), pp2.begin(), pp2.end());
		  
		// Check
		Point pmid = 0.5*pp1 + 0.5*pp2;
		Vector4D bary = bc.cartToBary(Vector3D(pmid[0], pmid[1], pmid[2]));
		double dist = implicit(bary);
		if (dist > 1.0e-4)
		  std::cout << "line dist: " << dist << ", ki= " << ki << ", kj= " << kj << std::endl;
	      }
	    
	    SISLCurve *qc2 = 0;
	    s1720(qc, 1, &qc2, &kstat);
		
	    int kcrv2=0, kpt2=0;
	    double *epar2 = 0;
	    SISLIntcurve **intcv2 = 0;
	    if (qc2)
	      s1871(qc2, &zero, 1, eps, &kpt2, &epar2, &kcrv2, &intcv2, &kstat);
	    double mindist = 1.0e8;
	    for (kr=0; kr<kpt2; ++kr)
	      {
		Vector4D barypt = (1.0 - epar2[kr])*bary1 + epar2[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		Vector3D pos = bc.baryToCart(barypt);
		der.insert(der.end(), pos.begin(), pos.end());
		double dist = implicit(barypt);
		mindist = std::min(dist, mindist);
		// if (fabs(dist) < std::max(maxdist,1.0e-7))
		//   der2.insert(der2.end(), pos.begin(), pos.end());
	      }
	    for (kr=0; kr<kpt2; ++kr)
	      {
		Vector4D barypt = (1.0 - epar2[kr])*bary1 + epar2[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		Vector3D pos = bc.baryToCart(barypt);
		double dist = implicit(barypt);
		mindist = std::min(dist, mindist);
		if (fabs(dist) < std::max(sigma_min,1.0e-8))//std::min(1.5*mindist,std::max(0.1*maxdist,1.0e-7)))
		  {
		    std::ofstream f1("iter.g2");
		    std::ofstream f2("iter.txt");
		    for (int kb=0; kb<50; ++kb)
		      {
			f1 << "400 1 0 4 255 0 0 255" << std::endl;
			f1 << "1" << std::endl;
			f1 << pos << std::endl;
			
			double d1 = deriv1(barypt);
			double d2 = deriv2(barypt);
			double d3 = deriv3(barypt);
			double d4 = deriv4(barypt);
			Vector4D dv(d1,d2,d3,d4);
			f2 << barypt << std::endl;
			f2 << dv << std::endl;
			f2 << dist << "  " << dv.length() << std::endl;
			double fac = 3000;
			// if (dist > 0.0)
			//   fac *= -1;
			barypt = barypt - fac*dist*dv;
			pos = bc.baryToCart(barypt);
			dist = implicit(barypt);
		      }
		    der2.insert(der2.end(), pos.begin(), pos.end());
		  }
	      }
	    for (kr=0; kr<kcrv2; ++kr)
	      {
		int ipt = intcv2[kr]->ipoint;
		double par1 = intcv2[kr]->epar1[0];
		double par2 = intcv2[kr]->epar1[ipt-1];
		Point pp1 = (1.0-par1)*cart1 + par1*cart2;
		Point pp2 = (1.0-par2)*cart1 + par2*cart2;
		lineder.insert(lineder.end(), pp1.begin(), pp1.end());
		lineder.insert(lineder.end(), pp2.begin(), pp2.end());
	      }
	    
	    if (qc) freeCurve(qc);
	    if (qc2) freeCurve(qc2);
	    if (intcv) freeIntcrvlist(intcv, kcrv);
	    if (intcv2) freeIntcrvlist(intcv2, kcrv2);
	    if (epar) free(epar);
	    if (epar2) free(epar2);
	  }
      }

    // Output
    if (points.size() > 0)
      {
	PointCloud3D ptcloud(&points[0], points.size()/3);
	outviz << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(outviz);
      }
    if (linesegs.size() > 0)
      {
	LineCloud lines(&linesegs[0], linesegs.size()/6);
	outviz << "410 1 0 4 255 0 0 255 " << std::endl;
	lines.write(outviz);
      }

    if (der.size() > 0)
      {
	PointCloud3D ptcloud(&der[0], der.size()/3);
	outviz << "400 1 0 4 0 255 0 255" << std::endl;
	ptcloud.write(outviz);
      }
    if (lineder.size() > 0)
      {
	LineCloud lines(&lineder[0], lineder.size()/6);
	outviz << "410 1 0 4 0 255 0 255 " << std::endl;
	lines.write(outviz);
      }
    if (der2.size() > 0)
      {
	PointCloud3D ptcloud(&der2[0], der2.size()/3);
	outviz << "400 1 0 4 100 155 0 255" << std::endl;
	ptcloud.write(outviz);
      }

    LineCloud tmp2(&tmpline[0], tmpline.size()/6);
    tmp2.writeStandardHeader(tmp);
    tmp2.write(tmp);
    
    return 0;
}

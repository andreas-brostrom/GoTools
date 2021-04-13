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


#include "GoTools/lrsplines3D/LRFeature3DUtils.h"

#include <fstream>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

//#define DEBUG


//===========================================================================
void  LRFeature3DUtils::writeCellInfo(const LRSplineVolume& vol, 
				      double tol, int ncell1, int ncell2, 
				      int ncell3, std::ostream &out)
//===========================================================================
{
  double u1 = vol.startparam_u();
  double u2 = vol.endparam_u();
  double v1 = vol.startparam_v();
  double v2 = vol.endparam_v();
  double w1 = vol.startparam_w();
  double w2 = vol.endparam_w();
  double udel = (u2 - u1)/(double)ncell1;
  double vdel = (v2 - v1)/(double)ncell2;
  double wdel = (w2 - w1)/(double)ncell3;
  int ix = (vol.dimension() == 3) ? 5 : 3;
  int del = 4 + vol.dimension();  // Parameter triple, position and
  // distance between volume and point

  int nc2 = ncell1*ncell2*ncell3;
  vector<int> nmb_pts(3*nc2, 0);
  vector<double> cellinfo(18*nc2, 0.0);
  int *npt = &nmb_pts[0];
  int *nout_over = npt+nc2;
  int *nout_under = nout_over+nc2;
  double *avdist = &cellinfo[0];
  double *max_over = avdist+nc2;
  double *max_under = max_over+nc2;
  double *minheight = max_under+nc2;
  double *maxheight = minheight+nc2;
  double *nel = maxheight+nc2;
  double *stdd = nel+nc2;
  double *slope = stdd+nc2;
  double *avheight = slope+nc2;
  double *stddheight = avheight+nc2;
  double *av2 = stddheight+nc2;
  double *av_over = av2+nc2;
  double *av_under = av_over+nc2;
  double *av_height_vol = av_under+nc2;
  double *min_height_vol = av_height_vol+nc2;
  double *max_height_vol = min_height_vol+nc2;
  double *laplace = max_height_vol+nc2;

  vector<int> nhit(nc2, 0);
  
  std::cout <<"Starting feature output. Non-parallel" << std::endl;
  // std::ofstream of("height-slope.txt");
  
  // Adjust default
  for (int kr=0; kr<nc2; ++kr)
    {
      max_over[kr] = std::numeric_limits<double>::lowest();
      max_under[kr] = std::numeric_limits<double>::max();
      maxheight[kr] = std::numeric_limits<double>::lowest();
      minheight[kr] = std::numeric_limits<double>::max();
      max_height_vol[kr] = std::numeric_limits<double>::lowest();
      min_height_vol[kr] = std::numeric_limits<double>::max();
    }

  // Traverse elements and relate point and element information to
  // cell
  int nelem = 0;
  std::cout << "Number of elements: " << vol.numElements() << std::endl;
  for (LRSplineVolume::ElementMap::const_iterator it=vol.elementsBegin();
       it != vol.elementsEnd(); ++it, ++nelem)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      double wel1 = it->second->wmin();
      double wel2 = it->second->wmax();
      int i1 = (int)((uel1 - u1)/udel);
      int i2 = (int)((uel2 - u1)/udel);
      i1 = std::max(0, std::min(ncell1-1, i1));
      i2 = std::max(0, std::min(ncell1-1, i2));
      int j1 = (int)((vel1 - v1)/vdel);
      int j2 = (int)((vel2 - v1)/vdel);
      j1 = std::max(0, std::min(ncell2-1, j1));
      j2 = std::max(0, std::min(ncell2-1, j2));
      int k1 = (int)((wel1 - w1)/wdel);
      int k2 = (int)((wel2 - w1)/wdel);
      k1 = std::max(0, std::min(ncell3-1, k1));
      k2 = std::max(0, std::min(ncell3-1, k2));

      // Compute element fraction in cell
      for (int ki=i1; ki<=i2; ++ki)
	{
	  double ufrac = (std::min(uel2, u1+(ki+1)*udel) -
			  std::max(uel1, u1+ki*udel))/(uel2 - uel1);
	  for (int kj=j1; kj<=j2; ++kj)
	    {
	      double vfrac = (std::min(vel2, v1+(kj+1)*vdel) -
			      std::max(vel1, v1+kj*vdel))/(vel2 - vel1);
	      for (int kr=k1; kr<=k2; ++kr)
		{
		  double wfrac = (std::min(wel2, w1+(kr+1)*wdel) -
				  std::max(wel1, w1+kr*wdel))/(wel2 - wel1);
		  nel[(kr*ncell2+kj)*ncell1+ki] += (ufrac*vfrac*wfrac);
		}
	    }
	}

     // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell1-1, i1));
	  j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell2-1, j1));
	  k1 = (int)((points[kr*del+2] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  int ixcell = (k1*ncell2+j1)*ncell1 + i1;
	  av2[ixcell] += points[kr*del+ix+1];
	  avdist[ixcell] += fabs(points[kr*del+ix+1]);
	  avheight[ixcell] += points[kr*del+ix];
	  max_over[ixcell] = std::max(max_over[ixcell], points[kr*del+ix+1]);
	  max_under[ixcell] = std::min(max_under[ixcell], points[kr*del+ix+1]);
	  maxheight[ixcell] = std::max(maxheight[ixcell], points[kr*del+ix]);
	  minheight[ixcell] = std::min(minheight[ixcell], points[kr*del+ix]);
	  npt[ixcell]++;
	  if (points[kr*del+ix+1] > tol)
	    nout_over[ixcell]++;
	  if (points[kr*del+ix+1] > 0.0)
	    av_over[ixcell] += points[kr*del+ix+1];
	  if (points[kr*del+ix+1] < -tol)
	  if (points[kr*del+ix+1] < -tol)
	    nout_under[ixcell]++;
	  if (points[kr*del+ix+1] < 0.0)
	    av_under[ixcell] -= points[kr*del+ix+1];
	}

      // Define parameter values for computation of slope, laplacian and volume
      // height information in element
      int i3 = (int)((uel1 - u1)/udel);
      int i4 = (int)((uel2 - u1)/udel);
      int j3 = (int)((vel1 - v1)/vdel);
      int j4 = (int)((vel2 - v1)/vdel);
      int k3 = (int)((wel1 - w1)/wdel);
      int k4 = (int)((wel2 - w1)/wdel);
      int nsample = 2;
      int ns3 = nsample*nsample*nsample;
      double udel2 = udel/(double)nsample;
      double vdel2 = vdel/(double)nsample;
      double wdel2 = wdel/(double)nsample;
      vector<double> paru, parv, parw;
      paru.reserve(nsample*(i4-i3+1));
      parv.reserve(nsample*(j4-j3+1));
      parw.reserve(nsample*(k4-k3+1));
      double upar2 = u1 + i3*udel + 0.5*udel2;
      while (upar2-udel2 > uel1)
      	upar2 -= udel2;
      while (upar2 <= uel1)
      	upar2 += udel2;
      for (/*upar2=std::max(upar2,uel1)*/; upar2 <= uel2; upar2+=udel2)
	paru.push_back(upar2);
      double vpar2 = v1 + j3*vdel + 0.5*vdel2;
      while (vpar2-vdel2 > vel1)
      	vpar2 -= vdel2;
      while (vpar2 <= vel1)
      	vpar2 += vdel2;
      for (/*vpar2=std::max(vpar2,vel1)*/; vpar2 <= vel2; vpar2+=vdel2)
	parv.push_back(vpar2);
      double wpar2 = w1 + k3*wdel + 0.5*wdel2;
      while (wpar2-wdel2 > wel1)
     	wpar2 -= wdel2;
      while (wpar2 <= wel1)
     	wpar2 += wdel2;
       for (/*wpar2=std::max(wpar2,wel1)*/; wpar2 <= wel2; wpar2+=wdel2)
	parw.push_back(wpar2);

      // Compute point and derivative information
      vector<double> der;
      vol.elementGridEvaluate(it->second.get(), &paru[0], paru.size(),
			      &parv[0], parv.size(),
      			      &parw[0], parw.size(), 2, der);
      // vol.elementGridEvaluate(it->second.get(), paru, parv,
      // 			      parw, 1, der);
      //std::cout << "Grid evaluate end, nelem=" << nelem << std::endl;

      // Compute slope, laplacian and height and distribute to associated
      // cell
      size_t ka, kb, kc, kd;
      for (kc=0, kd=0; kc<parw.size(); ++kc)
	{
	  k1 = (int)((parw[kc] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  for (kb=0; kb<parv.size(); ++kb)
	    {
	      j1 = (int)((parv[kb] - v1)/vdel);
	      j1 = std::max(0, std::min(ncell2-1, j1));
	      for (ka=0; ka<paru.size(); ++ka, kd+=10)
		{
		  double volheight = der[kd];
		  // vector<Point> tmppt(4);
		  // vol.point(tmppt, paru[ka], parv[kb], parw[kc], 1);
		  //volheight = tmppt[0][0];
		  double slope2 = sqrt(der[kd+1]*der[kd+1]+der[kd+2]*der[kd+2]+
				       der[kd+3]*der[kd+3]);
		   double laplace2 = der[kd+4]*der[kd+4]+der[kd+7]*der[kd+7]+
		     der[kd+9]*der[kd+9];
		  // double slope3 = sqrt(tmppt[1][0]*tmppt[1][0]+tmppt[2][0]*tmppt[2][0]+
		  // 		       tmppt[3][0]*tmppt[3][0]);
		  // of << volheight << " " << tmppt[0][0] << " " << slope2 << " " << slope3 << std::endl;
		  i1 = (int)((paru[ka] - u1)/udel);
		  i1 = std::max(0, std::min(ncell1-1, i1));
		  int ixcell = (k1*ncell2+j1)*ncell1 + i1;
		  
		  slope[ixcell] += slope2/(double)ns3;
		  laplace[ixcell] += laplace2/(double)ns3;
		  av_height_vol[ixcell] += volheight/(double)ns3;
		  min_height_vol[ixcell] =
		    std::min(min_height_vol[ixcell], volheight);
		  max_height_vol[ixcell] =
		    std::max(max_height_vol[ixcell], volheight);
		  nhit[ixcell]++;
		}
	    }
	}
    }

  // std::cout << "hit: ";
  // for (size_t ki2=0; ki2<nhit.size(); ++ki2)
  //   std::cout << nhit[ki2] << " ";
  // std::cout << std::endl;
  
  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  av2[kr] /= (double)npt[kr];
	  avdist[kr] /= (double)npt[kr];
	  avheight[kr] /= (double)npt[kr];
	  av_over[kr] /= (double)npt[kr];
	  av_under[kr] /= (double)npt[kr];
	}
    }

  // Standard deviation
  for (LRSplineVolume::ElementMap::const_iterator it=vol.elementsBegin();
       it != vol.elementsEnd(); ++it)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      double wel1 = it->second->wmin();
      double wel2 = it->second->wmax();

      // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  int i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell1-1, i1));
	  int j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell2-1, j1));
	  int k1 = (int)((points[kr*del+2] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  int ixcell = (k1*ncell2+j1)*ncell1 + i1;
	  double tmp = points[kr*del+ix+1] - av2[ixcell];
	  stdd[ixcell] += tmp*tmp;
	  double tmp2 = points[kr*del+ix] - avheight[ixcell];
	  stddheight[ixcell] += tmp2*tmp2;
	}
    }

  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  stdd[kr] /= (double)npt[kr];
	  stddheight[kr] /= (double)npt[kr];
	}
    }

  double upar, vpar, wpar;

  // Normalize to harmonize the different entries
  double minval = 0.0;
  double maxval = 10.0;
  vector<vector<double> > outval(18);
  for (int kh=0; kh<18; ++kh)
    outval[kh].resize(nc2);
  double currmin = std::numeric_limits<double>::max();
  double currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(slope[kr], currmin);
      currmax = std::max(slope[kr], currmax);
      outval[0][kr] = slope[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[0][kr] = (currmax < 1.0e-9) ? 0 : outval[0][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_height_vol[kr], currmin);
      currmax = std::max(av_height_vol[kr], currmax);
      outval[1][kr] = av_height_vol[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[1][kr] = (currmax < 1.0e-9) ? 0 : outval[1][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double diff = max_height_vol[kr] - min_height_vol[kr];
      currmin = std::min(diff, currmin);
      currmax = std::max(diff, currmax);
      outval[2][kr] = diff;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[2][kr] = (currmax < 1.0e-9 || nhit[kr] == 0) ? 0 : outval[2][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avdist[kr], currmin);
      currmax = std::max(avdist[kr], currmax);
      outval[3][kr] = avdist[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[3][kr] = (currmax < 1.0e-9) ? 0 : outval[3][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 :
	std::max(fabs(max_over[kr]), fabs(max_under[kr]));
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[4][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[4][kr] = (currmax < 1.0e-9) ? 0 : outval[4][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avheight[kr], currmin);
      currmax = std::max(avheight[kr], currmax);
      outval[5][kr] = avheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    {
      if (currmin < 0)
	outval[5][kr] -= currmin;
      outval[5][kr] =(currmax < 1.0e-9) ? 0 :  outval[5][kr]*maxval/currmax;
    }

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : maxheight[kr] - minheight[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[6][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[6][kr] =(currmax < 1.0e-9) ? 0 :  outval[6][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stdd[kr], currmin);
      currmax = std::max(stdd[kr], currmax);
      outval[7][kr] = stdd[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[7][kr] = (currmax < 1.0e-9) ? 0 : outval[7][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stddheight[kr], currmin);
      currmax = std::max(stddheight[kr], currmax);
      outval[8][kr] = stddheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[8][kr] = (currmax < 1.0e-9) ? 0 : outval[8][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double div = std::max(fabs(max_under[kr]), fabs(max_over[kr]));
      double tmp = (npt[kr] == 0 || div < 1.0e-9) ? 0.0 : 
	avdist[kr]/div;
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[9][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[9][kr] = (currmax < 1.0e-9) ? 0 : outval[9][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : 
	std::max(max_over[kr], 0.0) - std::min(max_under[kr], 0.0);
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[10][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[10][kr] = (currmax < 1.0e-9) ? 0 : outval[10][kr]*maxval/currmax;
  
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_over[kr], currmin);
      currmax = std::max(av_over[kr], currmax);
      outval[11][kr] = av_over[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[11][kr] = (currmax < 1.0e-9) ? 0 : outval[11][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_under[kr], currmin);
      currmax = std::max(av_under[kr], currmax);
      outval[12][kr] = av_under[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[12][kr] = (currmax < 1.0e-9) ? 0 : outval[12][kr]*maxval/currmax;
      
   currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_under[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[13][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[13][kr] = (currmax < 1.0e-9) ? 0 : outval[13][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_over[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[14][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[14][kr] = (currmax < 1.0e-9) ? 0 : outval[14][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(nel[kr], currmin);
      currmax = std::max(nel[kr], currmax);
      outval[15][kr] = nel[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[15][kr] = (currmax < 1.0e-9) ? 0 : outval[15][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(laplace[kr], currmin);
      currmax = std::max(laplace[kr], currmax);
      outval[16][kr] = laplace[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[16][kr] = (currmax < 1.0e-9) ? 0 : outval[16][kr]*maxval/currmax;

 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min((double)npt[kr], currmin);
      currmax = std::max((double)npt[kr], currmax);
      outval[17][kr] = (double)npt[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[17][kr] = (currmax < 1.0e-9) ? 0 : outval[17][kr]*maxval/currmax;

 
  // Write to file
  out << ncell1 << "  " << ncell2 << "  " << ncell3 << " " << "17" << std::endl;
  for (int kr=0; kr<nc2; ++kr)
    {
      for (int kh=0; kh<18; ++kh)
	out << (float)outval[kh][kr] << " ";
      out << std::endl;
    }

  std::cout <<"Finished feature output" << std::endl;
}

//===========================================================================
void  LRFeature3DUtils::writeCellInfo_omp(const LRSplineVolume& vol, 
					  double tol, int ncell1, int ncell2, 
					  int ncell3, std::ostream &out)
//===========================================================================
{
  double u1 = vol.startparam_u();
  double u2 = vol.endparam_u();
  double v1 = vol.startparam_v();
  double v2 = vol.endparam_v();
  double w1 = vol.startparam_w();
  double w2 = vol.endparam_w();
  double udel = (u2 - u1)/(double)ncell1;
  double vdel = (v2 - v1)/(double)ncell2;
  double wdel = (w2 - w1)/(double)ncell3;
  int ix = (vol.dimension() == 3) ? 5 : 3;
  int del = 4 + vol.dimension();  // Parameter triple, position and
  // distance between volume and point

  printf("(%f,%f) x (%f,%f) x (%f,%f)\n",u1,u2,v1,v2,w1,w2);
  int nc2 = ncell1*ncell2*ncell3;
  vector<int> nmb_pts(3*nc2, 0);
  vector<double> cellinfo(18*nc2, 0.0);
  int *npt = &nmb_pts[0];
  int *nout_over = npt+nc2;
  int *nout_under = nout_over+nc2;
  double *avdist = &cellinfo[0];
  double *max_over = avdist+nc2;
  double *max_under = max_over+nc2;
  double *minheight = max_under+nc2;
  double *maxheight = minheight+nc2;
  double *nel = maxheight+nc2;
  double *stdd = nel+nc2;
  double *slope = stdd+nc2;
  double *avheight = slope+nc2;
  double *stddheight = avheight+nc2;
  double *av2 = stddheight+nc2;
  double *av_over = av2+nc2;
  double *av_under = av_over+nc2;
  double *av_height_vol = av_under+nc2;
  double *min_height_vol = av_height_vol+nc2;
  double *max_height_vol = min_height_vol+nc2;
  double *laplace = max_height_vol+nc2;

  std::cout <<"Starting feature output. Parallel" << std::endl;
  
  // Adjust default
  for (int kr=0; kr<nc2; ++kr)
    {
      max_over[kr] = std::numeric_limits<double>::lowest();
      max_under[kr] = std::numeric_limits<double>::max();
      maxheight[kr] = std::numeric_limits<double>::lowest();
      minheight[kr] = std::numeric_limits<double>::max();
      max_height_vol[kr] = std::numeric_limits<double>::lowest();
      min_height_vol[kr] = std::numeric_limits<double>::max();
    }

  // Traverse elements and relate point and element information to
  // cell
  int nelem = 0;
  std::cout << "Number of elements: " << vol.numElements() << std::endl;
  for (LRSplineVolume::ElementMap::const_iterator it=vol.elementsBegin();
       it != vol.elementsEnd(); ++it, ++nelem)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      double wel1 = it->second->wmin();
      double wel2 = it->second->wmax();
      int i1 = (int)((uel1 - u1)/udel);
      int i2 = (int)((uel2 - u1)/udel);
      i1 = std::max(0, std::min(ncell1-1, i1));
      i2 = std::max(0, std::min(ncell1-1, i2));
      int j1 = (int)((vel1 - v1)/vdel);
      int j2 = (int)((vel2 - v1)/vdel);
      j1 = std::max(0, std::min(ncell2-1, j1));
      j2 = std::max(0, std::min(ncell2-1, j2));
      int k1 = (int)((wel1 - w1)/wdel);
      int k2 = (int)((wel2 - w1)/wdel);
      k1 = std::max(0, std::min(ncell3-1, k1));
      k2 = std::max(0, std::min(ncell3-1, k2));

      // Compute element fraction in cell
      for (int ki=i1; ki<=i2; ++ki)
	{
	  double ufrac = (std::min(uel2, u1+(ki+1)*udel) -
			  std::max(uel1, u1+ki*udel))/(uel2 - uel1);
	  for (int kj=j1; kj<=j2; ++kj)
	    {
	      double vfrac = (std::min(vel2, v1+(kj+1)*vdel) -
			      std::max(vel1, v1+kj*vdel))/(vel2 - vel1);
	      for (int kr=k1; kr<=k2; ++kr)
		{
		  double wfrac = (std::min(wel2, w1+(kr+1)*wdel) -
				  std::max(wel1, w1+kr*wdel))/(wel2 - wel1);
		  nel[(kr*ncell2+kj)*ncell1+ki] += (ufrac*vfrac*wfrac);
		}
	    }
	}

     // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell1-1, i1));
	  j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell2-1, j1));
	  k1 = (int)((points[kr*del+2] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  int ixcell = (k1*ncell2+j1)*ncell1 + i1;
	  av2[ixcell] += points[kr*del+ix+1];
	  avdist[ixcell] += fabs(points[kr*del+ix+1]);
	  avheight[ixcell] += points[kr*del+ix];
	  max_over[ixcell] = std::max(max_over[ixcell], points[kr*del+ix+1]);
	  max_under[ixcell] = std::min(max_under[ixcell], points[kr*del+ix+1]);
	  maxheight[ixcell] = std::max(maxheight[ixcell], points[kr*del+ix]);
	  minheight[ixcell] = std::min(minheight[ixcell], points[kr*del+ix]);
	  npt[ixcell]++;
	  if (points[kr*del+ix+1] > tol)
	    nout_over[ixcell]++;
	  if (points[kr*del+ix+1] > 0.0)
	    av_over[ixcell] += points[kr*del+ix+1];
	  if (points[kr*del+ix+1] < -tol)
	  if (points[kr*del+ix+1] < -tol)
	    nout_under[ixcell]++;
	  if (points[kr*del+ix+1] < 0.0)
	    av_under[ixcell] -= points[kr*del+ix+1];
	}
    }

 
  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  av2[kr] /= (double)npt[kr];
	  avdist[kr] /= (double)npt[kr];
	  avheight[kr] /= (double)npt[kr];
	  av_over[kr] /= (double)npt[kr];
	  av_under[kr] /= (double)npt[kr];
	}
    }
  
  vector<LRSplineVolume::ElementMap::const_iterator> el1_vec;
  int num_elem = vol.numElements();
  el1_vec.reserve(num_elem);
  int max_num_bsplines = 0;
  for (LRSplineVolume::ElementMap::const_iterator iter = vol.elementsBegin();
       iter != vol.elementsEnd(); ++iter)
  {
      el1_vec.push_back(iter);
  }
  LRSplineVolume::ElementMap::const_iterator el1;
  int nsample = 2;
  int ns3 = nsample*nsample*nsample;
  double udel2 = udel/(double)nsample;
  double vdel2 = vdel/(double)nsample;
  double wdel2 = wdel/(double)nsample;

  vector<double> elem_slope(num_elem*nc2, 0.0);
  vector<double> elem_laplace(num_elem*nc2, 0.0);
  vector<double> elem_avheight(num_elem*nc2, 0.0);
  vector<double> elem_minheight(num_elem*nc2, std::numeric_limits<double>::max());
  vector<double> elem_maxheight(num_elem*nc2, std::numeric_limits<double>::lowest());
  vector<double> paru, parv, parw;
  int a1, a2, a3;
  printf("Starting parallel processing %d \n", num_elem*nc2);
#pragma omp parallel default(none) private(nelem, el1, paru, parv, parw, a1, a2, a3) shared(vol, el1_vec, num_elem, u1, u2, v1, v2, w1, w2, udel, vdel, wdel, udel2, vdel2, wdel2, nsample, ncell1, ncell2, ncell3, nc2, ns3, elem_slope, elem_laplace, elem_avheight, elem_minheight, elem_maxheight)
  {
    int i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4;
    double uel1, uel2, vel1, vel2, wel1, wel2;
    double upar2, vpar2, wpar2;
    vector<double> der;
    size_t ka, kb, kc, kd;
    double volheight, slope2, laplace2;
    int ixcell;
    #pragma omp for schedule(auto)//guided)//static,8)//runtime)//dynamic,4)
    //#pragma omp for schedule(static,2)
    for (nelem=0; nelem<num_elem;++nelem)
    {
      el1 = el1_vec[nelem];
      
      // Identify overlapping cells
      uel1 = el1->second->umin();
      uel2 = el1->second->umax();
      vel1 = el1->second->vmin();
      vel2 = el1->second->vmax();
      wel1 = el1->second->wmin();
      wel2 = el1->second->wmax();

      // Define parameter values for computation of slope, laplacian and volume
      // height information in element
      i3 = (int)((uel1 - u1)/udel);
      i3 = std::max(0, std::min(ncell1-1, i3));
      i4 = (int)((uel2 - u1)/udel);
      i4 = std::max(0, std::min(ncell1-1, i4));
      j3 = (int)((vel1 - v1)/vdel);
      j3 = std::max(0, std::min(ncell2-1, j3));
      j4 = (int)((vel2 - v1)/vdel);
      j4 = std::max(0, std::min(ncell2-1, j4));
      k3 = (int)((wel1 - w1)/wdel);
      k3 = std::max(0, std::min(ncell3-1, k3));
      k4 = (int)((wel2 - w1)/wdel);
      k4 = std::max(0, std::min(ncell3-1, k4));
      paru.resize(nsample*(i4-i3+1));
      parv.resize(nsample*(j4-j3+1));
      parw.resize(nsample*(k4-k3+1));
      upar2 = u1 + i3*udel;
      while (upar2-udel2 >= uel1)
	upar2 -= udel2;
      for (a1=0, upar2=std::max(upar2,uel1); upar2 <= uel2; upar2+=udel2, ++a1)
	paru[a1] = upar2;
      if (a1 > (int)paru.size())
	printf("Mismatch error length u %d %d \n",a1,paru.size());
      vpar2 = v1 + j3*vdel;
      while (vpar2-vdel2 >= vel1)
	vpar2 -= vdel2;
      for (a2=0, vpar2=std::max(vpar2,vel1); vpar2 <= vel2; vpar2+=vdel2, ++a2)
	parv[a2] = vpar2;
      if (a2 > (int)parv.size())
	printf("Mismatch error length v %d %d \n",a2,parv.size());
      wpar2 = w1 + k3*wdel;
      while (wpar2-wdel2 >= wel1)
	wpar2 -= wdel2;
      for (a3=0, wpar2=std::max(wpar2,wel1); wpar2 <= wel2; wpar2+=wdel2, ++a3)
	parw[a3] = wpar2;
      if (a3 > (int)parw.size())
	printf("Mismatch error length w %d %d \n",a3,parw.size());

      //printf("Evaluate: %i\n", nelem);

      // Compute po and derivative information
      vol.elementGridEvaluate(el1->second.get(), &paru[0], a1,
			      &parv[0], a2, &parw[0], a3, 2, der);

      // vol.elementGridEvaluate(el1->second.get(), paru, parv,
      // 			      parw, 1, der);
      //std::cout << "Grid evaluate end, nelem=" << nelem << std::endl;

      // Compute slope, laplacian and height and distribute to associated
      // cell
     for (kc=0, kd=0; kc<parw.size(); ++kc)
	{
	  k1 = (int)((parw[kc] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  for (kb=0; kb<parv.size(); ++kb)
	    {
	      j1 = (int)((parv[kb] - v1)/vdel);
	      j1 = std::max(0, std::min(ncell2-1, j1));
	      for (ka=0; ka<paru.size(); ++ka, kd+=10)
		{
		  volheight = der[kd];
		  slope2 = sqrt(der[kd+1]*der[kd+1]+der[kd+2]*der[kd+2]+
				der[kd+3]*der[kd+3]);
		   laplace2 = der[kd+4]*der[kd+4]+der[kd+7]*der[kd+7]+
		     der[kd+9]*der[kd+9];
		   i1 = (int)((paru[ka] - u1)/udel);
		   i1 = std::max(0, std::min(ncell1-1, i1));
		   ixcell = (k1*ncell2+j1)*ncell1 + i1;
		  
		   elem_slope[nelem*nc2+ixcell] += slope2/(double)ns3;
		   elem_laplace[nelem*nc2+ixcell] += laplace2/(double)ns3;
		   elem_avheight[nelem*nc2+ixcell] += volheight/(double)ns3;
		   elem_minheight[nelem*nc2+ixcell] =
		     std::min(elem_minheight[nelem*nc2+ixcell], volheight);
		   elem_maxheight[nelem*nc2+ixcell] =
		     std::max(elem_maxheight[nelem*nc2+ixcell], volheight);
		}
	    }
	}
    }
  }

 printf("Finished parallel \n");
 
// Sequential 
   for (int kr=0; kr<nc2; ++kr)
    {
      for (int ka=0; ka<num_elem; ++ka)
	{
	  slope[kr] += elem_slope[ka*nc2 + kr];
	  laplace[kr] += elem_laplace[ka*nc2 + kr];
	  av_height_vol[kr] += elem_avheight[ka*nc2 + kr];
	  min_height_vol[kr] = std::min(min_height_vol[kr], elem_minheight[ka*nc2+kr]);
	  max_height_vol[kr] = std::max(max_height_vol[kr], elem_maxheight[ka*nc2+kr]);
	}
    }
   
  // Standard deviation
  for (LRSplineVolume::ElementMap::const_iterator it=vol.elementsBegin();
       it != vol.elementsEnd(); ++it)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      double wel1 = it->second->wmin();
      double wel2 = it->second->wmax();

      // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  int i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell1-1, i1));
	  int j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell2-1, j1));
	  int k1 = (int)((points[kr*del+2] - w1)/wdel);
	  k1 = std::max(0, std::min(ncell3-1, k1));
	  int ixcell = (k1*ncell2+j1)*ncell1 + i1;
	  double tmp = points[kr*del+ix+1] - av2[ixcell];
	  stdd[ixcell] += tmp*tmp;
	  double tmp2 = points[kr*del+ix] - avheight[ixcell];
	  stddheight[ixcell] += tmp2*tmp2;
	}
    }

  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  stdd[kr] /= (double)npt[kr];
	  stddheight[kr] /= (double)npt[kr];
	}
    }

  double upar, vpar, wpar;

  // Normalize to harmonize the different entries
  double minval = 0.0;
  double maxval = 10.0;
  vector<vector<double> > outval(17);
  for (int kh=0; kh<17; ++kh)
    outval[kh].resize(nc2);
  double currmin = std::numeric_limits<double>::max();
  double currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(slope[kr], currmin);
      currmax = std::max(slope[kr], currmax);
      outval[0][kr] = slope[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[0][kr] = (currmax < 1.0e-9) ? 0 : outval[0][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_height_vol[kr], currmin);
      currmax = std::max(av_height_vol[kr], currmax);
      outval[1][kr] = av_height_vol[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[1][kr] = (currmax < 1.0e-9) ? 0 : outval[1][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double diff = max_height_vol[kr] - min_height_vol[kr];
      currmin = std::min(diff, currmin);
      currmax = std::max(diff, currmax);
      outval[2][kr] = diff;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[2][kr] = (currmax < 1.0e-9) ? 0 : outval[2][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avdist[kr], currmin);
      currmax = std::max(avdist[kr], currmax);
      outval[3][kr] = avdist[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[3][kr] = (currmax < 1.0e-9) ? 0 : outval[3][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 :
	std::max(fabs(max_over[kr]), fabs(max_under[kr]));
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[4][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[4][kr] = (currmax < 1.0e-9) ? 0 : outval[4][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avheight[kr], currmin);
      currmax = std::max(avheight[kr], currmax);
      outval[5][kr] = avheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    {
      if (currmin < 0)
	outval[5][kr] -= currmin;
      outval[5][kr] =(currmax < 1.0e-9) ? 0 :  outval[5][kr]*maxval/currmax;
    }

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : maxheight[kr] - minheight[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[6][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[6][kr] =(currmax < 1.0e-9) ? 0 :  outval[6][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stdd[kr], currmin);
      currmax = std::max(stdd[kr], currmax);
      outval[7][kr] = stdd[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[7][kr] = (currmax < 1.0e-9) ? 0 : outval[7][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stddheight[kr], currmin);
      currmax = std::max(stddheight[kr], currmax);
      outval[8][kr] = stddheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[8][kr] = (currmax < 1.0e-9) ? 0 : outval[8][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double div = std::max(fabs(max_under[kr]), fabs(max_over[kr]));
      double tmp = (npt[kr] == 0 || div < 1.0e-9) ? 0.0 : 
	avdist[kr]/div;
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[9][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[9][kr] = (currmax < 1.0e-9) ? 0 : outval[9][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : 
	std::max(max_over[kr], 0.0) - std::min(max_under[kr], 0.0);
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[10][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[10][kr] = (currmax < 1.0e-9) ? 0 : outval[10][kr]*maxval/currmax;
  
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_over[kr], currmin);
      currmax = std::max(av_over[kr], currmax);
      outval[11][kr] = av_over[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[11][kr] = (currmax < 1.0e-9) ? 0 : outval[11][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_under[kr], currmin);
      currmax = std::max(av_under[kr], currmax);
      outval[12][kr] = av_under[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[12][kr] = (currmax < 1.0e-9) ? 0 : outval[12][kr]*maxval/currmax;
      
   currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_under[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[13][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[13][kr] = (currmax < 1.0e-9) ? 0 : outval[13][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_over[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[14][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[14][kr] = (currmax < 1.0e-9) ? 0 : outval[14][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(nel[kr], currmin);
      currmax = std::max(nel[kr], currmax);
      outval[15][kr] = nel[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[15][kr] = (currmax < 1.0e-9) ? 0 : outval[15][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(laplace[kr], currmin);
      currmax = std::max(laplace[kr], currmax);
      outval[16][kr] = laplace[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[16][kr] = (currmax < 1.0e-9) ? 0 : outval[16][kr]*maxval/currmax;

 
  // Write to file
  out << ncell1 << "  " << ncell2 << "  " << ncell3 << " " << "17" << std::endl;
  for (int kr=0; kr<nc2; ++kr)
    {
      for (int kh=0; kh<17; ++kh)
	out << outval[kh][kr] << " ";
      out << std::endl;
    }

  std::cout <<"Finished feature output" << std::endl;
}


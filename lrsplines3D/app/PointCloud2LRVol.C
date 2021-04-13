#include <iostream>
#include <fstream>
#include <string.h>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines3D/LRVolApprox.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"

using namespace std;
using namespace Go;

void print_help_text()
{
  std::cout << "Purpose: Approximate a trivariate point cloud by an LR B-spline volume. \n";
  std::cout << "Mandatory parameters: input point cloud, output volume (.g2), tolerance, number of iterations. \n";
  std::cout << "An adaptive approximation procedure is applied which for the";
  std::cout << " specified number of iterations: \n";
  std::cout << " - Approximates the points with a volume in the current spline space \n";
  std::cout << " - Computes the approximation accuracy \n";
  std::cout << " - Refines the volume in areas where the tolerance is not met \n";
  std::cout << "The approximation is completed when the tolerance is met or the";
  std::cout << " specified number of iterations is exceeded.";
  std::cout << "The number of iterations is recommended to lie in the interval [4:7]. \n";
  std::cout << "The first line in the point file reports on the number of points.\n";
   std::cout << "The points follow and is to be given as x, y, z and w.\n";
  std::cout << "Optional input parameters: \n";
  std::cout << "-verbose <0/1>: Detailed output, default 0 \n";
  std::cout << "-par <dim>: Dimension of geometry space, used for parameterized data (default 1) \n";
  std::cout << "-magnitude  <0/1>: Approximate length of point vector \n";
  std::cout << "-dist <filename (.txt)> : Write distance field to file (x, y, z, distance) \n";
  std::cout << "-info <filename> : Write accuracy information to file \n";
  std::cout << "-initmba <0/1>: 0 = initiate with least squares method \n";
  std::cout << "                1 = apply only multilevel B-spline approximation (MBA) (Default) \n";
  std::cout << "-degree <2/3> : degree of polynomial segments, default = 2\n";
  std::cout << "-minsize <length>: Minimum element size (all directions) \n";
  std::cout << "-tolfile: File specifying domains with specific tolerances, global tolerance apply outside domains. PointCloud2LR -tolfile for file format \n";
  std::cout << "-toldoc: Documentation on file format for tolerance domains. \n";
  std::cout << "-outfrac <percentage>: Local measure for when the fraction of points outside the tolerance should lead to volume splitting \n";
  std::cout << "-feature <ncell1> <ncell2> <ncell3>: Specify 3D grid for feature output \n";
  std::cout << "-featurelevels <number of levels> <level 1> ... <level n> \n";
  std::cout << "-bezout <0/1>: Bezier file output for modhot vizualization, default 0 \n";
  std::cout << "-surfin <filename.g2>: LR spline volume to update \n";
  std::cout << "-spacein <filename.g2>: LR spline volume from which to fetch initial spline space \n";
  std::cout << "surfin and spacein can not be applied simultanously \n";
  std::cout << "-combineinout <0/1>: Combine resulting volume and volume given with -spacein or -surfin in one higher-dimensional volume. Requires number of iterations = 0 \n";
  std::cout << "-h or --help : Write this text\n";
}

void print_tol_file_format()
{
  std::cout << "File specifying domains/boxes with different tolerances. \n";
  std::cout << "If not otherwise stated, the global tolerance will applies outside the boxes \n";
  std::cout << "Numbers are given as floats and separated by space \n";
  std::cout << "Line1: Number of boxes (integer), whether or not the tolerance is specified outside the boxes (0/1) \n";
  std::cout << "Following lines: \n";
  std::cout << "xmin xmax ymin ymax zmin zmax tolerance \n";
  std::cout << "Positive numbers for tolerance means absolute value, \n";
  std::cout << "negative numbers mean multiplication factor to standard deviation of points \n";
  std::cout << "Last line (if given): Tolerance, positive or negative float as for previous lines. Value overrules global tolerance. \n";
  std::cout << "Ensure non-overlapping boxes. No test applied. \n";
}

int fetchIntParameter(int argc, char *argv[], int ki, int& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atoi(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchDoubleParameter(int argc, char *argv[], int ki, double& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atof(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchCharParameter(int argc, char *argv[], int ki, char*& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = argv[ki+1];
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int main (int argc, char *argv[]) {

  char *pointfile = 0;     // Input point file
  char *volfile = 0;       // Volume output file
  char *infofile = 0;      // Accuracy information output file
  char *tolfile = 0;       // File specifying varying tolerances
  char *field_out = 0;     // Distance field output file
  char *invol = 0;         // 

  double epsge;
  int levels;
  int initMBA = 1;
  int del = 4;
  int degree = 2;
  int verbose = 0;
  double minsize = -1.0;
  double outfrac = 0.0;
  int ncell1=0, ncell2=0, ncell3=0;
  bool features = false;
  vector<int> feature_levels;
  int dim = 1;
  int magnitude = 0;
  int bezout = 0;
  int initvol = 0;
  int combine = 0;

  int ki, kj;
  vector<bool> par_read(argc-1, false);

  // Read optional parameters
  int nmb_par = argc-1;
  for (ki=1; ki<argc; ++ki)
    {
      string arg(argv[ki]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
      else if (arg == "-toldoc")
	{
	  print_tol_file_format();
	  exit(0);
	}
      else if (arg == "-dist")
	{
	  int stat = fetchCharParameter(argc, argv, ki, field_out, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-info")
	{
	  int stat = fetchCharParameter(argc, argv, ki, infofile, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-initmba")
	{
	  int stat = fetchIntParameter(argc, argv, ki, initMBA, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-verbose")
	{
	  int stat = fetchIntParameter(argc, argv, ki, verbose, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-par")
	{
	  int stat = fetchIntParameter(argc, argv, ki, dim, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-magnitude")
	{
	  int stat = fetchIntParameter(argc, argv, ki, magnitude, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-degree")
	{
	  int stat = fetchIntParameter(argc, argv, ki, degree, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-minsize")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, minsize, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
       else if (arg == "-tolfile")
	{
	  int stat = fetchCharParameter(argc, argv, ki, tolfile, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-outfrac")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, outfrac, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  outfrac = std::max(0.0, outfrac);
	  outfrac /= 100.0;
	}
      else if (arg == "-feature")
	{
	  if (ki == argc-1)
	    {
	      std::cout << "ERROR: Missing input" << std::endl;
	      print_help_text();
	      return 1;
	    }
	  features = true;
	  ncell1 = atoi(argv[ki+1]);
	  ncell2 = atoi(argv[ki+2]);
	  ncell3 = atoi(argv[ki+3]);
	  par_read[ki-1] = par_read[ki] = par_read[ki+1] = par_read[ki+2] = true;
	  nmb_par -= 4;
	}
     else if (arg == "-featurelevels")
	{
	  if (ki == argc-1)
	    {
	      std::cout << "ERROR: Missing input" << std::endl;
	      print_help_text();
	      return 1;
	    }
	  int fsize = atoi(argv[ki+1]);
	  par_read[ki-1] = par_read[ki] = true;
	  feature_levels.resize(fsize);
	  for (int ka=0; ka<fsize; ++ka)
	    {
	      feature_levels[ka] = atoi(argv[ki+ka+2]);
	      par_read[ki+ka+1] = true;
	    }
	  nmb_par -= (fsize+2);
	}
     else if (arg == "-bezout")
	{
	  int stat = fetchIntParameter(argc, argv, ki, bezout, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
     else if (arg == "-surfin")
       {
	  int stat = fetchCharParameter(argc, argv, ki, invol, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  initvol = 1;
       }
     else if (arg == "-spacein")
       {
	  int stat = fetchCharParameter(argc, argv, ki, invol, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  initvol = 2;
       }
     else if (arg == "-combineinout")
       {
	  int stat = fetchIntParameter(argc, argv, ki, combine, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
       }
    }

  // Read remaining parameters
  if (nmb_par != 4)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

   for (ki=1; ki<argc; ++ki)
    {
      if (par_read[ki-1])
	continue;
      if (nmb_par == 4)
	{
	  pointfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 3)
	{
	  volfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  epsge = atof(argv[ki]);
	  nmb_par--;
	}
      else
	{
	  levels = atoi(argv[ki]);
	}
    }

   if ((!invol) || levels > 0)
     combine = 0;

   // Read data points
   ifstream ifs(pointfile);
  int nmb_pts;
  ifs >> nmb_pts;

  int dim2 = (magnitude) ? 1 : dim;
  vector<double> pc4d;

  double domain[6];
  domain[0] = domain[2] = domain[4] = 1.0e8;
  domain[1] = domain[3] = domain[5] = -1.0e8;
  vector<double> mba_level(dim2,0.0); //0.5*(minval+maxval); //0.0;
  vector<double> minval(dim2, std::numeric_limits<double>::max());
  vector<double> maxval(dim2, std::numeric_limits<double>::lowest());
  vector<double> tmpq(dim);
  for (int ix=0; ix!=nmb_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      
      domain[0] = std::min(domain[0], p0);
      domain[1] = std::max(domain[1], p0);
      domain[2] = std::min(domain[2], p1);
      domain[3] = std::max(domain[3], p1);
      domain[4] = std::min(domain[4], p2);
      domain[5] = std::max(domain[5], p2);

      for (int ka=0; ka<dim; ++ka)
	{
	  ifs >> q0;
	  tmpq[ka] = q0;
	}

      if (magnitude)
	{
	  Point tmp(tmpq.begin(), tmpq.end());
	  //double tmp = Utils::sum_squared(tmpq.begin(), tmpq.end());
	  tmpq[0] = tmp.length(); //sqrt(tmp);
	}

      for (int ka=0; ka<dim2; ++ka)
	{
	  pc4d.push_back(tmpq[ka]);

	  mba_level[ka] += tmpq[ka]/(double)nmb_pts;
	  minval[ka] = std::min(minval[ka], tmpq[ka]);
	  maxval[ka] = std::max(maxval[ka], tmpq[ka]);
	}
    }

 
  bool use_stdd = false;
  vector<LRVolApprox::TolBox> tolerances;
  if (tolfile != 0)
    {
      std::ifstream tolin(tolfile);
      int nmb_box, last;
      double umin, umax, vmin, vmax, wmin, wmax, tol;
      tolin >> nmb_box >> last;
      tolerances.resize(nmb_box);
      for (int ka=0; ka<nmb_box; ++ka)
	{
	  tolin >> umin >> umax >> vmin >> vmax >> wmin >> wmax;
	  tolin >> tol;
	  std::cout << "Tol: " << tol << std::endl;
	  tolerances[ka].setVal(std::max(umin,domain[0]), std::min(umax,domain[1]),
				std::max(vmin,domain[2]), std::min(vmax,domain[3]),
				std::max(wmin,domain[4]), std::min(wmax,domain[5]), tol);
	  std::cout << "Box tolerance: " << tolerances[ka].tol << std::endl;
	  if (tol < 0.0)
	    use_stdd = true;
	}
      if (last > 0)
	{
	  tolin >> epsge;
	  if (epsge < 0.0)
	    use_stdd = true;
	}
    }
  if (epsge < 0.0)
    use_stdd = true;

  double mba_level2 = 0.0; //(dim > 1) ? 0.0 : mba_level[0];
  
  if (use_stdd)
    {
      double stdd = 0.0;
      for (ki=0; ki<nmb_pts; ++ki)
	{
	  double q1 = pc4d[del*ki+del-1];
	  stdd += (pow(mba_level2-q1,2)/(double)nmb_pts);
	}
      stdd = sqrt(stdd);
      for (size_t kj=0; kj<tolerances.size(); ++kj)
	{
	  if (tolerances[kj].tol < 0.0)
	    tolerances[kj].setTol(fabs(tolerances[kj].tol)*stdd);
	  std::cout << "Box tolerance2: " << tolerances[kj].tol << std::endl;
	}
      if (epsge < 0.0)
	epsge = fabs(epsge)*stdd;

      std::cout << "Standard deviation: " << stdd << std::endl;
     }
     
  std::cout << "Domain: [" << domain[0] << "," << domain[1] << "]x[" << domain[2];
  std::cout << "," << domain[3] << "]x[" << domain[4] << "," << domain[5] << "]" << std::endl;
  for (int ka=0; ka<dim2; ++ka)
    std::cout << "Range(" << ka << "): [" << minval[ka] << "," << maxval[ka] << "]" << std::endl;
  shared_ptr<LRVolApprox> vol_approx;
  shared_ptr<LRSplineVolume> volin;
  if (invol && initvol > 0)
    {
      std::ifstream ifs(invol);
      GoTools::init();
      ObjectHeader oh;
      oh.read(ifs);
      volin = shared_ptr<LRSplineVolume>(new LRSplineVolume());
      volin->read(ifs);

      shared_ptr<LRSplineVolume> volspace = volin;
      if (initvol == 2)
	{
	  // Set coefficients to zero and weights to one. Update dimension if
	  // necessary
	  if (combine)
	    {
	      // Must keep coefficient information in input volume
	      volspace = shared_ptr<LRSplineVolume>(volin->clone());
	    }
	  Point coef(dim2);
	  coef.setValue(0.0);
	  for (LRSplineVolume::BSplineMap::const_iterator it1 = volspace->basisFunctionsBegin();
	       it1 != volspace->basisFunctionsEnd(); ++it1)
	    {
	      double gamma = it1->second->gamma();
	      volspace->setCoefAndDim(coef, gamma, it1->second.get());
	    }
	}
      vol_approx = shared_ptr<LRVolApprox>(new LRVolApprox(volspace, pc4d, epsge,
							   false, false, true));
    }
  else
    {
      int ncoef = 6; //6; //8; //6
      int order = degree+1; //5; //3
      int nm = ncoef*ncoef*ncoef;
      double dom = (domain[1]-domain[0])*(domain[3]-domain[2])*(domain[5]-domain[4]);
      double c1 = std::pow((double)nm/dom, 1.0/3.0);
      int nc[3];
      for (kj=0; kj<3; ++kj)
	{
	  double tmp = c1*(domain[2*kj+1]-domain[2*kj]);
	  nc[kj] = (int)tmp;
	  //if (tmp - (double)nc[kj] > 0.5)
	  ++nc[kj];
	  nc[kj] = std::max(nc[kj], order);
	}
      std::cout << "Number of coefficients: " << nc[0] << ", " << nc[1] << ", " << nc[2] << std::endl;
      vol_approx = shared_ptr<LRVolApprox>(new LRVolApprox(nc[0], order, nc[1], order, nc[2], order,
							   pc4d, dim2, domain, epsge, mba_level2));
      // vol_approx = shared_ptr<LRVolApprox>(new LRVolApprox(ncoef, order, ncoef, order,
      // 							   ncoef, order, pc4d, dim2, domain,
      // 							   epsge, mba_level2));
      vol_approx->setInitMBA(initMBA);
    }

  if (tolerances.size() > 0)
    vol_approx->setVarTolBox(tolerances);
  if (minsize > 0.0)
    vol_approx->setMinimumElementSize(minsize, minsize, minsize);
  if (outfrac > 0.0)
    vol_approx->setOutFraction(outfrac);
  if (verbose)
    vol_approx->setVerbose(true);
  else
    vol_approx->setVerbose(false);

  // Feature output
  if (features)
    {
      vol_approx->setFeatureOut(ncell1, ncell2, ncell3);
      vol_approx->setFeatureLevel(feature_levels);
    }
  
  double max, average, av_all;
  double maxout, avout;
  int num_out;
  cout << "Starting approximation..." << endl;

  shared_ptr<LRSplineVolume> result = vol_approx->getApproxVol(max,av_all,average,num_out,levels);

  vol_approx->fetchOutsideTolInfo(maxout, avout);

  if (combine)
    {
      // Combine initial volume with result in one volume
      Point coef(volin->dimension()+result->dimension());
      LRSplineVolume::BSplineMap::const_iterator it1 = volin->basisFunctionsBegin();
      LRSplineVolume::BSplineMap::const_iterator it2 = result->basisFunctionsBegin();
      for (; it1 != volin->basisFunctionsEnd(); ++it1, ++it2)
	{
	  Point cf1 = it1->second->Coef();
	  double gamma = it1->second->gamma();
	  Point cf2 = it2->second->Coef();
	  int ki, kr;
	  for (ki=0, kr=0; ki<volin->dimension(); ++ki, ++kr)
	    coef[kr] = cf1[ki];
	  for (ki=0; ki<result->dimension(); ++ki, ++kr)
	    coef[kr] = cf2[ki];
	  result->setCoefAndDim(coef, gamma, it2->second.get());
	}
    }
  
  if (infofile)
    {
      std::ofstream infoout(infofile);   // Accuracy information output stream

      infoout << "Total number of points: " << nmb_pts << std::endl;
      infoout << "Number of elements: " << result->numElements() << std::endl;
      infoout << "Maximum distance: " << max << std::endl;
      infoout << "Average distance: " << av_all << std::endl;
      infoout << "Average distance for points outside of the tolerance: " << average << std::endl;
      infoout << "Number of points outside the tolerance: " << num_out << std::endl;
      infoout << "Maximum distance exceeding tolerance (dist-tol): " << maxout << std::endl;
      infoout << "Average distance exceeding tolerance (dist-tol): " << avout << std::endl;
    } 
  else
    {
      std::cout << "Total number of points: " << nmb_pts << std::endl;
      std::cout << "Number of elements: " << result->numElements() << std::endl;
      std::cout << "Maximum distance: " << max << std::endl;
      std::cout << "Average distance: " << av_all << std::endl;
      std::cout << "Average distance for points outside of the tolerance: " << average << std::endl;
      std::cout << "Number of points outside the tolerance: " << num_out << std::endl;
      std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout << std::endl;
      std::cout << "Average distance exceeding tolerance (dist-tol): " << avout << std::endl;
     } 
  ofstream ofs(volfile);
  result->writeStandardHeader(ofs);
  result->write(ofs);

  if (bezout)
    {
      LRSpline3DBezierCoefs bez(*result);
      
      bez.getBezierCoefs();

      char bezfile[200];
      int len = (int)strlen(volfile);
      strcpy(bezfile, volfile);
      int len2 = (int)strlen(bezfile);
      strncat(bezfile,".bb", 3);
      bez.writeToFile(bezfile);
    }

  if (field_out)
    {
      // Fetch data points with distance information
      vector<double> pnts_dist;
      pnts_dist.reserve(5*nmb_pts);
      try {
	LRSplineVolume::ElementMap::const_iterator elem = result->elementsBegin();
	LRSplineVolume::ElementMap::const_iterator last = result->elementsEnd();
	for (; elem != last; ++elem)
	  {
	    if (!elem->second->hasDataPoints())
	      continue;
	    vector<double>& points = elem->second->getDataPoints();
	    pnts_dist.insert(pnts_dist.end(), points.begin(), points.end());
	  }
      }
      catch (...)
	{
	  std::cout << "ERROR: Extraction of distance information failed" << std::endl;
	  return 1;
	}

      // Write to file
      std::ofstream field_info(field_out);
      (void)field_info.precision(15);
      for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	{
	  for (ki=0; ki<4; ++ki)
	    field_info << pnts_dist[kj+ki] << " ";
	  field_info << std::endl;
	}
    }

}

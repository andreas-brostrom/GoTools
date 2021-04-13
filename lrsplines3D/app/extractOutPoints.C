#include <iostream>
#include <fstream>
#include <string.h>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/Point.h"

using namespace std;
using namespace Go;

int main (int argc, char *argv[]) {

  if (argc != 6) {
    cout << "usage: ./extractOutsidePoints <input 4d pt cloud(.txt)> <dim> <input LR field> <level value> <output txt file>" << endl;
    return -1;
  }

  ifstream ifs1(argv[1]);
  int dim = atoi(argv[2]);
  ifstream ifs2(argv[3]);
  double level = atof(argv[4]);
  char* outfile(argv[5]);

  GoTools::init();
  ObjectHeader oh;
  oh.read(ifs2);
  shared_ptr<LRSplineVolume> volin(new LRSplineVolume());
  volin->read(ifs2);

  Array<double,6> param = volin->parameterSpan();
  
  int num_pts;
  ifs1 >> num_pts;
  vector<double> outpts;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2;
      vector<double> qval(dim);
      ifs1 >> p0 >> p1 >> p2;
      for (int ka=0; ka<dim; ++ka)
	{
	  ifs1 >> qval[ka];
	}

      if (p0 < param[0] || p0 > param[1] || p1 < param[2] ||
	  p1 > param[3] || p2 < param[4] || p2 > param[5])
	continue;   // Outside
      
      Point pt;
      volin->point(pt, p0, p1, p2);
      if (pt.length() < level)
	{
	  outpts.push_back(p0);
	  outpts.push_back(p1);
	  outpts.push_back(p2);
	  for (int ka=0; ka<dim; ++ka)
	    outpts.push_back(qval[ka]);
	}
    }

  std::ofstream ofs(outfile);
  std::streamsize prev = ofs.precision(15);
  int nmbout = (int)outpts.size()/(3+dim);
  ofs <<  nmbout << std::endl;
  for (int ka=0; ka<(int)outpts.size(); ka+=(3+dim))
    {
      for (int kb=0; kb<3+dim; ++kb)
	ofs << outpts[ka+kb] << " ";
      ofs << std::endl;
    }


  
}

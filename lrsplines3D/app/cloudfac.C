#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;

int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main (int argc, char *argv[]) {

  if (argc != 9) {
    cout << "usage: ./clodfac <input 4d pt cloud(.g2)> <point dimension> <input density lrspline(.g2)> <density level> <input velocity lrspline> <output distance out> <output distance in> <no levels>" << endl;
    return -1;
  }

  ifstream ifs1(argv[1]);
  int ptdim = atoi(argv[2]);
  ifstream ifs2(argv[3]);
  double dlev = atof(argv[4]);
  ifstream ifs3(argv[5]);
  ofstream ofs1(argv[6]);
  ofstream ofs2(argv[7]);

  int level = atoi(argv[8]);


  int num_pts;
  ifs1 >> num_pts;

  GoTools::init();

  vector<double> pc4din, pc4dout;

  ObjectHeader oh;
  oh.read(ifs2);

  shared_ptr<LRSplineVolume> vol1(new LRSplineVolume());
  vol1->read(ifs2);
  int dim1 = vol1->dimension();

  oh.read(ifs3);

  shared_ptr<LRSplineVolume> vol2(new LRSplineVolume());
  vol2->read(ifs3);
  int dim2 = vol2->dimension();

  vector<double> minh(ptdim, std::numeric_limits<double>::max());
  vector<double> maxh(ptdim, std::numeric_limits<double>::lowest());
  double maxmagn = 0.0;
  double mindin = 1.0e8, maxdin = -1.0e8;
  double mindout = 1.0e8, maxdout = -1.0e8;
  double avdin = 0.0, avdout = 0.0;
  int numin = 0, numout = 0;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs1 >> p0 >> p1 >> p2;
      Point pd;
      vol1->point(pd, p0, p1, p2);
      if (pd[0] < dlev)
	{
	  pc4dout.push_back(p0);
	  pc4dout.push_back(p1);
	  pc4dout.push_back(p2);
	}
      else
	{
	  pc4din.push_back(p0);
	  pc4din.push_back(p1);
	  pc4din.push_back(p2);
	}

      Point ptval(ptdim);
      for (int ka=0; ka<ptdim; ++ka)
	{
	  ifs1 >> q0;
	  if (pd[0] < dlev)
	    pc4dout.push_back(q0);
	  else
	    pc4din.push_back(q0);
	  ptval[ka] = q0;
	  minh[ka] = std::min(minh[ka], q0);
	  maxh[ka] = std::max(maxh[ka], q0);
	  ptval[ka] = q0;
	}
      double magn = ptval.length();
      maxmagn = std::max(maxmagn, magn);
      
      Point pos;
      vol2->point(pos, p0, p1, p2);
      double dist = 0.0;
      if (ptdim == 1 && dim2 == 1)
	dist = q0 - pos[0];
      else if (ptdim == dim2)
	dist = pos.dist(ptval);
      else if (dim2 == 1)
	{
	  double len = ptval.length();
	  dist = fabs(pos[0] - len);
	}
      
       if (pd[0] < dlev)
	{
	  mindout = std::min(mindout, dist);
	  maxdout = std::max(maxdout, dist);
	  avdout += fabs(dist);
	  numout++;
	}
       else
 	{
	  mindin = std::min(mindin, dist);
	  maxdin = std::max(maxdin, dist);
	  avdin += fabs(dist);
	  numin++;
	}
   }
  avdout /= (double)numout;
  avdin /= (double)numin;
  
  std::cout << "Min value: ";
  for (int ka=0; ka<ptdim; ++ka)
    std::cout << minh[ka] << "  ";
  std::cout << std::endl <<"Max value: ";
  for (int ka=0; ka<ptdim; ++ka)
    std::cout << maxh[ka] << "  ";
  std::cout <<  std::endl;
  std::cout << "Max magnitude: " << maxmagn << std::endl;
  std::cout << "Max d out: " << maxdout << ", mind out: " << mindout << ", num out: " << numout << ", average: " << avdout << std::endl;
  std::cout << "Max d in: " << maxdin << ", mind in: " << mindin << ", num in: " << numin << ", average: " << avdin << std::endl;

  int ki;

  vector<double> limits_d(2*level+1);
  vector<vector<double> > level_din(2*level+2);
  vector<vector<double> > level_dout(2*level+2);

  double mind = std::min(mindout, mindin);
  double maxd = std::max(maxdout, maxdin);
  double mleveld = std::max(fabs(mind), fabs(maxd));
  double deld = mleveld/(double)level;
  limits_d[level] = 0;
  for (ki=1; ki<=level; ++ki)
    {
      limits_d[level-ki] = -ki*deld;
      limits_d[level+ki] = ki*deld;
    }

  int kj;
  double *curr = &pc4din[0];
  for (int ix=0; ix!=numin; ++ix, curr+=3+ptdim)
    {
      Point pos;
      vol2->point(pos, curr[0], curr[1], curr[2]);
      double dist = 0.0;
      if (ptdim == 1 && dim2 == 1)
	dist = curr[3] - pos[0];
      else if (ptdim == dim2)
	{
	  Point ptval(curr[3], curr[4], curr[5]);
	  dist = pos.dist(ptval);
	}
      else if (dim2 == 1)
	{
	  Point ptval(curr[3], curr[4], curr[5]);
	  double len = ptval.length();
	  dist = fabs(pos[0] - len);
	}

     // Find classifications
      for (kj=0; kj<limits_d.size(); ++kj)
	if (dist < limits_d[kj])
	  {
	    level_din[kj].push_back(curr[0]);
	    level_din[kj].push_back(curr[1]);
	    level_din[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits_d.size())
	{
	  level_din[kj].push_back(curr[0]);
	  level_din[kj].push_back(curr[1]);
	  level_din[kj].push_back(curr[2]);
	}
    }

  curr = &pc4dout[0];
  for (int ix=0; ix!=numout; ++ix, curr+=3+ptdim)
    {
      Point pos;
      vol2->point(pos, curr[0], curr[1], curr[2]);
      double dist = 0.0;
      if (ptdim == 1 && dim2 == 1)
	dist = curr[3] - pos[0];
      else if (ptdim == dim2)
	{
	  Point ptval(curr[3], curr[4], curr[5]);
	  dist = pos.dist(ptval);
	}
      else if (dim2 == 1)
	{
	  Point ptval(curr[3], curr[4], curr[5]);
	  double len = ptval.length();
	  dist = fabs(pos[0] - len);
	}

     // Find classifications
      for (kj=0; kj<limits_d.size(); ++kj)
	if (dist < limits_d[kj])
	  {
	    level_dout[kj].push_back(curr[0]);
	    level_dout[kj].push_back(curr[1]);
	    level_dout[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits_d.size())
	{
	  level_dout[kj].push_back(curr[0]);
	  level_dout[kj].push_back(curr[1]);
	  level_dout[kj].push_back(curr[2]);
	}
    }

  // Write to files
  for (ki=0; ki<level_dout.size(); ++ki)
    {
      std::cout << "Level out [";
      if (ki == 0)
	std::cout << "-inf";
      else
	std::cout << limits_d[ki-1];
      if (ki < level_dout.size()-1)
	std::cout << ", " << limits_d[ki] << "]: ";
      else
	std::cout << ", inf]: ";
      std::cout << level_dout[ki].size()/3 << std::endl;
      if (level_dout[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_dout[ki].begin(), level_dout[ki].size()/3);

      double cc[3];
      if (ki <= level)
	{
	  cc[0] = ((level-ki)*colors[0][0] + ki*colors[1][0])/level;
	  cc[1] = ((level-ki)*colors[0][1] + ki*colors[1][1])/level;
	  cc[2] = ((level-ki)*colors[0][2] + ki*colors[1][2])/level;
	}
      else
	{
	  cc[0] = ((ki-level-1)*colors[2][0] + 
		   (2*level-ki+1)*colors[1][0])/level;
	  cc[1] = ((ki-level-1)*colors[2][1] + 
		   (2*level-ki+1)*colors[1][1])/level;
	  cc[2] = ((ki-level-1)*colors[2][2] + 
		   (2*level-ki+1)*colors[1][2])/level;
	}

      ofs1 << "400 1 0 4 " << cc[0] << " " << cc[1];
      ofs1 << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(ofs1);
    }

  for (ki=0; ki<level_din.size(); ++ki)
    {
      std::cout << "Level in [";
      if (ki == 0)
	std::cout << "-inf";
      else
	std::cout << limits_d[ki-1];
      if (ki < level_din.size()-1)
	std::cout << ", " << limits_d[ki] << "]: ";
      else
	std::cout << ", inf]: ";
      std::cout << level_din[ki].size()/3 << std::endl;
      if (level_din[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_din[ki].begin(), level_din[ki].size()/3);

      double cc[3];
      if (ki <= level)
	{
	  cc[0] = ((level-ki)*colors[0][0] + ki*colors[1][0])/level;
	  cc[1] = ((level-ki)*colors[0][1] + ki*colors[1][1])/level;
	  cc[2] = ((level-ki)*colors[0][2] + ki*colors[1][2])/level;
	}
      else
	{
	  cc[0] = ((ki-level-1)*colors[2][0] + 
		   (2*level-ki+1)*colors[1][0])/level;
	  cc[1] = ((ki-level-1)*colors[2][1] + 
		   (2*level-ki+1)*colors[1][1])/level;
	  cc[2] = ((ki-level-1)*colors[2][2] + 
		   (2*level-ki+1)*colors[1][2])/level;
	}

      ofs2 << "400 1 0 4 " << cc[0] << " " << cc[1];
      ofs2 << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(ofs2);
    }

}

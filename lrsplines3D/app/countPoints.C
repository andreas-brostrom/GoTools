#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;


int main (int argc, char *argv[]) {

  if (argc != 5) {
    cout << "usage: ./countPoints <input 4d pt cloud(.txt)> <dim> <output txt file> <average number of points in cell>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  int dim = atoi(argv[2]);
  ofstream ofs(argv[3]);
  int av_no = atoi(argv[4]);

  // std::cout << argv[1] << std::endl;
  // std::cout << dim << std::endl;
  // std::cout << argv[3] << std::endl;
  // std::cout << av_no << std::endl;
  int num_pts;
  ifs >> num_pts;

  vector<double> pc;
  double domain[6];
  domain[0] = domain[2] = domain[4] = 1.0e8;
  domain[1] = domain[3] = domain[5] = -1.0e8;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2;
      pc.push_back(p0);
      pc.push_back(p1);
      pc.push_back(p2);
      domain[0] = std::min(domain[0], p0);
      domain[1] = std::max(domain[1], p0);
      domain[2] = std::min(domain[2], p1);
      domain[3] = std::max(domain[3], p1);
      domain[4] = std::min(domain[4], p2);
      domain[5] = std::max(domain[5], p2);

      for (int ka=0; ka<dim; ++ka)
	{
	  ifs >> q0;
	}
    }
  
  std::cout << "Domain: [" << domain[0] << "," << domain[1] << "]x[" << domain[2];
  std::cout << "," << domain[3] << "]x[" << domain[4] << "," << domain[5] << "]" << std::endl;
  double dom = (domain[1]-domain[0])*(domain[3]-domain[2])*(domain[5]-domain[4]);
  double nm = (double)num_pts/(double)av_no;
  double c1 = std::pow(nm/dom, 1.0/3.0);
  double tmp = c1*(domain[1]-domain[0]);
  int nx = (int)tmp;
  nx = std::max(nx, 10);
  tmp = c1*(domain[3]-domain[2]);
  int ny = (int)tmp;
  ny = std::max(ny, 10);
  tmp = c1*(domain[5]-domain[4]);
  int nz = (int)tmp;
  nz = std::max(nz, 10);

  vector<double> lev_x(nx+1);
  vector<double> lev_y(ny+1);
  vector<double> lev_z(nz+1);
  double delx = (domain[1] - domain[0])/(double)(nx);
  double dely = (domain[3] - domain[2])/(double)(ny);
  double delz = (domain[5] - domain[4])/(double)(nz);
  lev_x[0] = domain[0];
  for (int ix=1; ix<=nx; ++ix)
    lev_x[ix] = lev_x[ix-1] + delx;
  lev_y[0] = domain[2];
  for (int iy=1; iy<=ny; ++iy)
    lev_y[iy] = lev_y[iy-1] + dely;
  lev_z[0] = domain[4];
  for (int iz=1; iz<=nz; ++iz)
    lev_z[iz] = lev_z[iz-1] + delz;
  
  vector<int> nmbpt(nx*ny*nz, 0);

  for (size_t ki=0; ki<pc.size(); ki+=3)
    {
      int ix;
      for (ix=0; ix<nx-1; ++ix)
	if (pc[ki] >= lev_x[ix] && pc[ki] < lev_x[ix+1])
	  break;
      ix = std::min(ix,nx-1);
      int iy;
      for (iy=0; iy<ny-1; ++iy)
	if (pc[ki+1] >= lev_y[iy] && pc[ki+1] < lev_y[iy+1])
	  break;
      iy = std::min(iy,ny-1);
       int iz;
      for (iz=0; iz<nz-1; ++iz)
	if (pc[ki+2] >= lev_z[iz] && pc[ki+2] < lev_z[iz+1])
	  break;
      iz = std::min(iz,nz-1);

      nmbpt[(iz*ny+iy)*nx+ix]++;
    }

  double av_no2 = 0.0;
  int max_no = 0;
  std::streamsize prev = ofs.precision(15);
  ofs << nx*ny*nz + 8 << std::endl;
  for (int kk=0; kk<nz; ++kk)
    {
      double midz = 0.5*(lev_z[kk]+lev_z[kk+1]);
      for (int kj=0; kj<ny; ++kj)
	{
	  double midy = 0.5*(lev_y[kj]+lev_y[kj+1]);
	  for (int ki=0; ki<nx; ++ki)
	    {
	      ofs << 0.5*(lev_x[ki]+lev_x[ki+1]) << " " << midy << " " << midz;
	      ofs << " " << nmbpt[(kk*ny+kj)*nx+ki] << std::endl;

	      av_no2 += (double)nmbpt[(kk*ny+kj)*nx+ki];
	      max_no = std::max(max_no, nmbpt[(kk*ny+kj)*nx+ki]);
	    }
	}
    }

  int currno = nmbpt[0] == 0 ? -1 : nmbpt[0];
  ofs << domain[0] << " " << domain[2] << " " << domain[4] << " " << currno << std::endl;
  currno = nmbpt[nx-1] == 0 ? -1 : nmbpt[nx-1];
  ofs << domain[1] << " " << domain[2] << " " << domain[4] << " " << currno << std::endl;
  currno = nmbpt[(ny-1)*nx] == 0 ? -1 : nmbpt[(ny-1)*nx];
  ofs << domain[0] << " " << domain[3] << " " << domain[4] << " " << currno << std::endl;
  currno = nmbpt[ny*nx-1] == 0 ? -1 : nmbpt[ny*nx-1];
  ofs << domain[1] << " " << domain[3] << " " << domain[4] << " " << currno << std::endl;
  currno = nmbpt[(nz-1)*nx*ny] == 0 ? -1 : nmbpt[(nz-1)*nx*ny];
  ofs << domain[0] << " " << domain[2] << " " << domain[5] << " " << currno << std::endl;
  currno = nmbpt[(nz-1)*nx*ny+nx-1] == 0 ? -1 : nmbpt[(nz-1)*nx*ny+nx-1];
  ofs << domain[1] << " " << domain[2] << " " << domain[5] << " " << currno << std::endl;
  currno = nmbpt[(nz-1)*nx*ny+(ny-1)*nx] == 0 ? -1 : nmbpt[(nz-1)*nx*ny+(ny-1)*nx];
  ofs << domain[0] << " " << domain[3] << " " << domain[5] << " " << currno << std::endl;
  currno = nmbpt[nx*ny*nz-1] == 0 ? -1 : nmbpt[nx*ny*nz-1];
  ofs << domain[1] << " " << domain[3] << " " << domain[5] << " " << currno << std::endl;
  
  av_no2 /= (double)(nx*ny*nz);
  std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
  std::cout << "Average no of pts: " << av_no2 << std::endl;
  std::cout << "Maximum no of pts: " << max_no << std::endl;
}

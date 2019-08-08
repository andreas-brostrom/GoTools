#include <fstream>
#include <iostream>
#include <chrono>

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"


using namespace std;
using namespace Go;

int main(int varnum, char* vararg[])
{

  if (varnum != 3) 
    {
      std::cout << "<LR B-spline surface in> <LR B-spline surface out>" << std::endl;
      exit(-1);
    }
  
  // Establishing test LR spline surface
  ifstream is(vararg[1]);
  ObjectHeader header;
  header.read(is);
  LRSplineSurface lrsurf(is);
  is.close();

  ofstream os(vararg[2]);

  RectDomain dom = lrsurf.parameterDomain();
  std::cout << "Domain: (" << dom.umin() << ", " << dom.umax() << ") x (";
  std::cout << dom.vmin() << ", " << dom.vmax() << ")" << std::endl;
  std::cout << "Reduced parameter domain (umin, umax, vmin, vmax):" << std::endl;
  double u1, u2, v1, v2;
  std::cin >> u1 >> u2 >> v1 >> v2;
  double fuzzy = 0.01;
  shared_ptr<LRSplineSurface> sub(lrsurf.subSurface(u1, v1, u2, v2, fuzzy));

  sub->writeStandardHeader(os);
  sub->write(os);
}

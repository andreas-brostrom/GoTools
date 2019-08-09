#include "GoTools/lrsplines2D/LRSplinePlotUtils2.h"
#include "GoTools/utils/checks.h"
#include <map>
#include <array>
#include <functional>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace Go;

extern map<array<int, 4>, int> boxdraw_map;

namespace {

struct Segment {int ix; int start; int end;};
typedef vector<Segment> segvec;

// given the gridpoint (x, y), return the multiplicities (0 or higher) of 
// meshrectangles meeting at this gridpoint.  The four components are, in order:
// multiplicity at right ([x, x+1]), multiplicity upwards ([y, y+1]), multiplicity
// at left, and multiplicity downwards.
array<int, 4> gridpoint_multiplicities(const Mesh2D& m, int x, int y) ;  

// return true if (x, y) is on the border of m in direction 'dir', where
// dir \in [0, 3] (representing right, up, left, down).
bool on_border(const Mesh2D& m, int x, int y, int dir);

// returns a vector containing the segments in 'input' that are terminated on 
// both sides by segments in 'terminators'.  The segments returned are removed
// from 'input'.
segvec terminated_by(segvec& input, const segvec& terminators);

void plot_segments_on_mesh(const Mesh2D& m, const segvec& xf_segs, const segvec& yf_segs);
  
// Convert all 'longest consecutive segments' in mesh m to the format 
// represented by the struct Segment.  (What is returned from m is only index
// pairs representing starts and ends, but does not give the index of the line the
// segment is localized on.  This information is added by the below function.
segvec all_segments(const Mesh2D& m, Direction2D d);

// Helper function for the 'supports_at_corner' function.
// compute the knots of the deduced knotvector perpendicular to the fixed direction,
// and with ortogonal segments of at least length 'orto_len'
vector<int> derived_knotvec(const Mesh2D&m,
			    Direction2D dfixed, int ix_fixed, 
			    int beg, int max_run, 
			    int orto_len, int numknots);

// Helper function for the 'supports_at_corner' function.  Computes the knotvecs
// associated with ortogonal directions up to max_orto, and returns the result.
// If 'result' is the resulting vector, then 'result[i]' is the knotvector associated
// with an ortogonal direction of length 'i'.  Naturally, result[0] is then undefined,
// but is present in the vector to make the indexing of the other elements more natural.
vector<vector<int> > compute_supports(const Mesh2D& m, 
				      Direction2D dfixed, int ix_fixed, int max_orto, int ix_start_run,
				      int max_run, int numknots);

// Helper function for the 'supports_at_corner' function.
// Identify all nonzero spans in the knotvec 'kvec' that cover 'knots + 1' knots.
// The results are returned as a vector of integer pairs.  The first integer in a returned pair
// represents the length of the span, whereas the second represents the start knot index of the span.
vector<pair<int, int> > spans(const vector<int>& kvec, int knots);

// Compute all supports in mesh 'm' that have their lower leftmost corner in (xpos, ypos), 
// and with number of knots equal to 'xknots' in the x-direction and 'yknots' in the y-direction.
vector<pair<vector<int>, vector<int> > > 
supports_at_corner(const Mesh2D& m, int xpos, int ypos, int xknots, int yknots);

// arguments of functor (x, y, dir), where dir \in [0, 3] (representing right, up, left, down).
void plot_mesh_generic(const Mesh2D& m, const function<bool(int, int, int)>&f);



}; // end anonymous namespace


namespace Go
{

// =============================================================================
void generate_mesh_texture(const Mesh2D& m, int u_pixel_res, int v_pixel_res, const char* filename)
// =============================================================================
{
  assert(u_pixel_res > 0 && v_pixel_res > 0);
  const int BL = 0;   // black - represents multiplicities > 1
  const int GR = 210; // grey  - represents multiplicity = 1
  const int WH = 255; // white - represents background
  const double v_span = m.maxParam(YFIXED) - m.minParam(YFIXED);
  const double u_span = m.maxParam(XFIXED) - m.minParam(XFIXED);
  const double half_pixel_width_v = v_span / (v_pixel_res - 1) / 2;
  const double half_pixel_width_u = u_span / (u_pixel_res - 1) / 2;
    
  vector<int> pixels; // we store the pixel values in a vector before writing them to file, due to differing
                      // direction of y-axis in images and in conventional parameter planes  
  pixels.reserve(u_pixel_res * v_pixel_res);
  for (int v_ix = 0, next_vknot = 0; v_ix != v_pixel_res; ++v_ix) {
    const double vix_mid = m.minParam(YFIXED) + (v_ix / double(v_pixel_res - 1)) * v_span;
    const bool on_meshline_v = (abs(vix_mid - m.kval(YFIXED, next_vknot)) <= half_pixel_width_v);
    for (int u_ix = 0, next_uknot = 0; u_ix != u_pixel_res; ++u_ix) {
      const double uix_mid = m.minParam(XFIXED) + (u_ix / double(u_pixel_res - 1)) * u_span;
      const bool on_meshline_u = (abs(uix_mid - m.kval(XFIXED, next_uknot)) <= half_pixel_width_u);
      const int pixel_value =
	on_meshline_u ? ((next_vknot == 0)                                         ? GR :
			 (m.nu(XFIXED, next_uknot, next_vknot-1, next_vknot) == 1) ? GR :
			 (m.nu(XFIXED, next_uknot, next_vknot-1, next_vknot) >  1) ? BL : WH) :
	on_meshline_v ? ((next_uknot == 0)                                         ? GR :
			 (m.nu(YFIXED, next_vknot, next_uknot-1, next_uknot) == 1) ? GR :
			 (m.nu(YFIXED, next_vknot, next_uknot-1, next_uknot) >  1) ? BL : WH) : 
	                WH;
      pixels.push_back(pixel_value);

      if (on_meshline_u && next_uknot < m.numDistinctKnots(XFIXED) - 1) ++next_uknot;
    }
    if (on_meshline_v && next_vknot < m.numDistinctKnots(YFIXED) - 1)  ++next_vknot;
  }

  // writing pixels to file
  ofstream os(filename);
  if (!os) throw runtime_error("generate_mesh_texture() : unable to open file.");
  else os << "P2 " << u_pixel_res << " " << v_pixel_res << " 255\n";
  for (int v_ix = v_pixel_res - 1; v_ix >= 0; --v_ix)
    for (int u_ix = 0; u_ix != v_pixel_res; ++u_ix) 
      os << pixels[v_ix * u_pixel_res + u_ix] << " ";
  os.close();
}

// =============================================================================
  void plot_mesh(const Mesh2D& m, int thick_threshold)
// =============================================================================
{
  plot_mesh_generic(m, 
		    [&](int x, int y, int dir) 
		        {return gridpoint_multiplicities(m, x, y)[dir] >= thick_threshold;});
}

// =============================================================================
void plot_all_supports(const Mesh2D& m, int xdeg, int ydeg)
// =============================================================================
{
  for (int y_ix = 0; y_ix != m.numDistinctKnots(YFIXED); ++y_ix)
    for (int x_ix = 0; x_ix != m.numDistinctKnots(XFIXED); ++x_ix)
      plot_supports_at_corner(m, x_ix, y_ix, xdeg, ydeg);
}

// =============================================================================
  void plot_rect_domain(const Mesh2D& m, int xmin, int ymin, int xmax, int ymax)
// =============================================================================
{
  plot_mesh_generic(m, 
		    [=](int x, int y, int dir)
		    {
		      if (on_border(m, x, y, dir)) return false;
		      return 
			(dir == 0 && (y == ymin || y == ymax) && nondecreasing(xmin, x, xmax-1)) ? true :
			(dir == 1 && (x == xmin || x == xmax) && nondecreasing(ymin, y, ymax-1)) ? true :
			(dir == 2 && (y == ymin || y == ymax) && nondecreasing(xmin+1, x, xmax)) ? true :
			(dir == 3 && (x == xmin || x == xmax) && nondecreasing(ymin+1, y, ymax)) ? true :
			false;
		    });
}

// // =============================================================================
// void plot_bspline_function(const Mesh2D& m, const BSplineFunction& b)
// // =============================================================================
// {
//   plot_support(m, 
// 	       &(b.kvec(XFIXED)[0]), 
// 	       &(b.kvec(YFIXED))[0], 
// 	       b.kvec(XFIXED).size(), 
// 	       b.kvec(YFIXED).size());
// }

// =============================================================================
void plot_support(const Mesh2D& m, const int* const kvec1, const int* const kvec2, int len1, int len2)
// =============================================================================
{
    wcout << "\nX-knotvec (indices, not values): "; 
    for (auto k = kvec1; k != kvec1 + len1; ++k) wcout << *k << " ";
    wcout << "\nY-knotvec (indices, not values): " ;
    for (auto k = kvec2; k != kvec2 + len2; ++k) wcout << *k << " ";

    plot_mesh_generic(m,
		      [&](int x, int y, int dir)
		      {
			if (on_border(m, x, y, dir)) return false;
			//auto& xknots = s.first;
			//auto& yknots = s.second;

			if ((dir%2 == 0) && 
			    (find(kvec2, kvec2 + len2, y) != kvec2 + len2))
			  return nondecreasing(kvec1[0], 
					       (dir == 0) ? x : x-1,
					       kvec1[len1 - 1] - 1);
			if ((dir%2 == 1) && (find (kvec1, kvec1 + len1, x) != kvec1 + len1))
			  return nondecreasing(kvec2[0],
					       (dir == 1) ? y : y-1,
					       kvec2[len2 - 1] - 1);
			// if we got this far, we are not located on a meshrectangle
			// that participates in the support
			return false;
		      });

}

// =============================================================================
void plot_supports_at_corner(const Mesh2D& m, int xpos, int ypos, int xdeg, int ydeg)
// =============================================================================
{
  auto supps = supports_at_corner(m, xpos, ypos, xdeg+2, ydeg+2);
  for (auto s : supps) {
    plot_support(m, &(s.first[0]), &(s.second[0]), s.first.size(), s.second.size());
  }
}


// =============================================================================
void plot_largest_tensorgrid(const Mesh2D& m)
// =============================================================================
{
  auto vx = orig_tensorgrid_knotpositions(m, XFIXED);
  auto vy = orig_tensorgrid_knotpositions(m, YFIXED);
  
  plot_mesh_generic(m, 
  		    [&m, &vx, &vy] (int x, int y, int dir)
		    { if (on_border(m, x, y, dir)) return false;
		      else return (dir == 1 || dir == 3) ?
			     find(vx.begin(), vx.end(), x) != vx.end() :
			     find(vy.begin(), vy.end(), y) != vy.end();});
}

// =============================================================================
void plot_history(const Mesh2D& m)
// =============================================================================
{
  const int xfnum = m.numDistinctKnots(XFIXED) - 1;
  const int yfnum = m.numDistinctKnots(YFIXED) - 1;
  segvec xf_all = all_segments(m, XFIXED);
  segvec yf_all = all_segments(m, YFIXED);
  segvec xf_active, yf_active;
  xf_active.push_back({0, 0, yfnum});
  xf_active.push_back({xfnum, 0, yfnum});
  yf_active.push_back({0, 0, xfnum});
  yf_active.push_back({yfnum, 0, xfnum});
  segvec xf_new, yf_new;

  for (int step = 0; true; step++) {
    xf_new = terminated_by(xf_all, yf_active);
    yf_new = terminated_by(yf_all, xf_active);
     if (xf_new.empty() && yf_new.empty()) return; // no further development possible

    wcout << "Now plotting history step " << step;
    wcout << ", with " << xf_new.size() + yf_new.size() << " new segments.";
    plot_segments_on_mesh(m, xf_new, yf_new);

    copy(xf_new.begin(), xf_new.end(), back_inserter(xf_active));
    copy(yf_new.begin(), yf_new.end(), back_inserter(yf_active));
  } 
}

// =============================================================================
vector<int> orig_tensorgrid_knotpositions(const Mesh2D& m, Direction2D d)
// =============================================================================
{
  vector<int> result;
  for (int i = 0; i != m.numDistinctKnots(d); ++i) {
    if (m.nu(d, i, 0, m.numDistinctKnots(flip(d)) - 1) > 0)
      result.push_back(i);
  }
  return result;
}



}; // end namespace LR_test


namespace {

// =============================================================================
array<int, 4> gridpoint_multiplicities(const Mesh2D& m, int x, int y)
// =============================================================================
{
  return array<int, 4> { m.nu(YFIXED, y, x, (x+1 < m.numDistinctKnots(XFIXED)) ? x+1 : x),
                         m.nu(XFIXED, x, y, (y+1 < m.numDistinctKnots(YFIXED)) ? y+1 : y),
                         m.nu(YFIXED, y, (x>0) ? x-1 : x, x),
                         m.nu(XFIXED, x, (y>0) ? y-1 : y, y)
                       };
}

// =============================================================================
void plot_segments_on_mesh(const Mesh2D& m, const segvec& xf_segs, const segvec& yf_segs)
// =============================================================================
{
  auto found = [](const segvec& sv, int ix, int pos)->bool 
  {
    for (auto seg : sv) 
      if (seg.ix == ix && nondecreasing(seg.start, pos, seg.end - 1)) return true;
    return false;
  };

  plot_mesh_generic(m,
		    [&] (int x, int y, int dir)
		    { if (on_border(m, x, y, dir)) return false;
		      else return (dir == 0) ? found(yf_segs, y, x) :
			          (dir == 1) ? found(xf_segs, x, y) :
			          (dir == 2) ? found(yf_segs, y, x-1) : // OK since not on border
			                       found(xf_segs, x, y-1);  // OK since not on border
		    });
}

// =============================================================================
void plot_mesh_generic(const Mesh2D& m, const std::function<bool(int, int, int)>&f)
// =============================================================================
{
  map<array<int, 4>, int>& bmap = boxdraw_map;
  int prev_right=0;
  wcout << endl;

  for (int y = m.numDistinctKnots(YFIXED) - 1; y >= 0; --y) {
    for (int x = 0; x != m.numDistinctKnots(XFIXED); ++x) {
      auto tmp = gridpoint_multiplicities(m, x, y);
      for (int i = 0; i != 4; ++i) { tmp[i] = f(x, y, i) ? 2 : (tmp[i] > 0) ? 1 : 0;}

      // drawing intermediate dash to get a better x/y aspect ratio
      if (x > 0) wcout << (wchar_t)bmap[array<int, 4>{tmp[2], 0, prev_right, 0}];

      // draw the corner point
      wcout << (wchar_t) bmap[tmp];
      prev_right = tmp[0];
    }
    wcout << '\n';
  }
}

// =============================================================================
bool on_border(const Mesh2D& m, int x, int y, int dir)
// =============================================================================
{
  return 
    (dir == 0) ? (x == m.numDistinctKnots(XFIXED) - 1) :
    (dir == 1) ? (y == m.numDistinctKnots(YFIXED) - 1) :
    (dir == 2) ? (x == 0) :
                 (y == 0); // (dir == 3) here
}

// =============================================================================
segvec terminated_by(segvec& input, const segvec& trm) // 'trm' refers to the 'terminators'
// =============================================================================
{
  auto terminates = [&](const Segment& t, const Segment& s, bool at_start) { 
    return t.ix == (at_start ? s.start : s.end) && nondecreasing(t.start, s.ix, t.end);
  };

  auto is_terminated = [&](Segment s) { 
    bool terminated_at_start = false;
    bool terminated_at_end = false;
    for (auto t : trm) {
      if (terminates(t, s, true))       terminated_at_start = true;
      else if (terminates(t, s, false)) terminated_at_end   = true;
      if (terminated_at_start && terminated_at_end) return true;
    }
    return false;
  };

  segvec result;
  segvec input_reduced;
  for (Segment s : input) 
    if (is_terminated(s)) result.push_back(s);
    else input_reduced.push_back(s);

  input.swap(input_reduced);
  return result;

}

// =============================================================================
segvec all_segments(const Mesh2D& m, Direction2D d)
// =============================================================================
{
  segvec result;
  vector<pair<int, int> > tmp;
  for (int ix = 0; ix != m.numDistinctKnots(d); ++ix) {
    tmp = m.segments(d, ix, 1);
    for (auto s : tmp) result.push_back({ix, s.first, s.second});
  }
  return result;
}

// =============================================================================
vector<int> derived_knotvec(const Mesh2D&m,
			      Direction2D dfixed, int ix_fixed, 
			      int beg, int max_run, 
			      int orto_len, int numknots)
// =============================================================================
{
  vector<int> result;
  int start_mul = 0; // will be set (and fixed) at first iteration
  for (int pos = beg; pos <= beg + max_run && result.size() < numknots + start_mul - 1; ++pos) {
    const int mult = m.nu(flip(dfixed), pos, ix_fixed, ix_fixed + orto_len);
    if (pos == beg) {
      start_mul = mult;
      if (start_mul == 0) return {}; // knot multipicity at start is zero.  No support possible.
    } 
    result.insert(result.end(), mult, pos);
  }
  // If resulting knotvector is insufficient for establishing a support, we return the empty vector.
  return (result.size() >= numknots + start_mul - 1) ? result : vector<int>();
}

// =============================================================================
vector<vector<int> > compute_supports(const Mesh2D& m, 
				      Direction2D dfixed, int ix_fixed, int max_orto, int ix_start_run,
				      int max_run, int numknots)
// =============================================================================
{
  vector<vector<int> > result(1); // first index is a dummy - only used to make other indexing easier 
  for (int cur_orto = 1; cur_orto <= max_orto; ++cur_orto) {
    result.push_back(derived_knotvec(m, dfixed, ix_fixed, ix_start_run, max_run, cur_orto, numknots));
  }
  return result;
}

// =============================================================================
vector<pair<int, int> > spans(const vector<int>& kvec, int knots)
// =============================================================================
{
  if (kvec.empty()) return {};
  vector<pair<int,int> > result;
  for (auto it = kvec.begin(), it2 = it + knots - 1; *it == kvec[0] && it2 < kvec.end(); ++it, ++it2)
    if (*it < *it2) result.push_back({{*it2 - *it}, int(it - kvec.begin())});
  return result;
}

// =============================================================================
vector<pair<vector<int>, vector<int> > > 
supports_at_corner(const Mesh2D& m, int xpos, int ypos, int xknots, int yknots)
// =============================================================================
{
  vector<pair<vector<int>, vector<int> > > result;
  const int x_extent = m.extent(YFIXED, ypos, xpos, 1); // number of consecutive non-zero intervals
  const int y_extent = m.extent(XFIXED, xpos, ypos, 1);
  
  vector<vector<int> > xvecs = compute_supports(m, YFIXED, ypos, y_extent, xpos, x_extent, xknots);
  vector<vector<int> > yvecs = compute_supports(m, XFIXED, xpos, x_extent, ypos, y_extent, yknots);
  
  for (int y_orto = 1; y_orto <= y_extent; ++y_orto) {
    for (auto x_op : spans(xvecs[y_orto], xknots)) {// all minimal x-support spans for the y-orto
      const int x_orto = x_op.first;
      auto y_op_vec = spans(yvecs[x_orto], yknots);
      auto y_op = find_if(y_op_vec.begin(), y_op_vec.end(), 
			  [y_orto](pair<int,int>& x) {return x.first == y_orto;});
      if (y_op != y_op_vec.end()) {
	auto xstart = &(xvecs[y_orto][x_op.second]);
	auto ystart = &(yvecs[x_orto][y_op->second]);
	result.push_back({{xstart, xstart + xknots}, {ystart, ystart + yknots}});
      }
    }
  }
  return result;
}


}; //end anonymous namespace 

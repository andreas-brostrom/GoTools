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

#ifndef _ENRICHEDTRIANG_H
#define _ENRICHEDTRIANG_H

#include "GoTools/utils/Point.h"
#include <vector>

namespace Go
{


/// \brief Class that represents a triangulation enriched with a scalar
/// value

  class Edge;
  class Triang;

  class Node
  {
    friend class Edge;
    friend class Triang;
    friend class EnrichedTriang;

  private:
    Point pos_;    // A trangulation node
    bool val_set_;
    double val_;   // Value in node
    std::vector<Triang*> triang_;  // Surrounding triangles
    std::vector<Edge*> edge_;    // Surrounding edges

    Node()
    {
      val_set_ = false;
    }

    Node(double *pos, int dim)
    {
      pos_ = Point(pos, pos+dim, true);
      val_set_ = false;
    }

    void addTriangle(Triang *triangle);

    std::vector<shared_ptr<Edge> > defineEdges();

    double getVal();

    void computeVal();

    Point getPos()
    {
      return pos_;
    }

    std::pair<Point,double> getPosAndVal();
  };

  class Edge
  {
    friend class Node;
    friend class Triang;
    friend class EnrichedTriang;

  private:
    Node* node_[2];

    Edge(Node *n1, Node *n2)
      {
	node_[0] = n1;
	node_[1] = n2;
      }

    double getVal(double fac)
    {
      return fac*node_[0]->getVal() + (1.0-fac)*node_[1]->getVal();
    }

    std::pair<Point,double> getPosAndVal(double fac);
  };

  class Triang
  {
    friend class Node;
    friend class Edge;
    friend class EnrichedTriang;

  private:
    Node* node_[3];
    double val_;

    Triang(Node *n1, Node *n2, Node *n3, double val);

    double getVal()
    {
      return val_;
    }

    std::pair<Point,double> getMid();

    std::vector<std::pair<Point,double> > getPosAndVal(int level);
  };

class EnrichedTriang
{
 public:
  /// Constructor given triangle and value information
  EnrichedTriang(int dimension, std::vector<double>& vertices, 
		 std::vector<int>& triangles, std::vector<double>& values);

  /// Destructor
  ~EnrichedTriang();

  /// Fetch position and value in nodes
  void evalNode(std::vector<std::pair<Point, double> >& nodeval);

  /// Fetch position and values in edges
  void evalEdge(int level, std::vector<std::pair<Point, double> >& edgeval);

  /// Evaluate position and values in triangles equally distributed in
  /// the inner, level gives number of rows : 1, 2 or 3
  void evalTriang(int level, std::vector<std::pair<Point, double> >& triangval);

  /// Evaluate position and values in the triangulation
  /// Level: 1 gives one point per triangle and nodes
  ///        2 gives 7 points per triangle, 3 per edge and nodes
  ///        3 gives 19 points per triangle, 5 per edge and nodes
  void evaluate(int level, std::vector<std::pair<Point, double> >& result);
  
  

 private:
  std::vector<shared_ptr<Node> > node_;  // Triangle nodes
  std::vector<shared_ptr<Edge> > edge_;  // Edges between triangle nodes
  std::vector<shared_ptr<Triang> > triang_; // Triangles


};

} // namespace Go



#endif  //  _ENRICHEDTRIANG_H

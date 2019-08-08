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

#include "GoTools/geometry/EnrichedTriang.h"

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;

//===========================================================================
EnrichedTriang::EnrichedTriang(int dimension, vector<double>& vertices, 
			       vector<int>& triangles, 
			       vector<double>& values)
//===========================================================================
{
  int nmb_nodes = (int)vertices.size()/dimension;

  // Prepare storage
  node_.reserve(nmb_nodes);
  triang_.reserve(triangles.size());

  // Store nodes
  for (int ki=0; ki<nmb_nodes; ++ki)
    node_.push_back(shared_ptr<Node>(new Node(&vertices[ki*dimension], 
  					      dimension)));

  // Store triangles
  int nmb_triang = (int)triangles.size()/3;
  for (int ki=0; ki<nmb_triang; ++ki)
    {
      triang_.push_back(shared_ptr<Triang>(new Triang(node_[triangles[3*ki]].get(),
						      node_[triangles[3*ki+1]].get(),
						      node_[triangles[3*ki+2]].get(),							    
						      values[ki])));
    }

  // Compute value information in nodes and make sure that all edges are
  // defined
  for (int ki=0; ki<nmb_nodes; ++ki)
    {
      node_[ki]->computeVal();
      vector<shared_ptr<Edge> > new_edges = node_[ki]->defineEdges();
      if (new_edges.size() > 0)
	edge_.insert(edge_.end(), new_edges.begin(), new_edges.end());
    }
}

//===========================================================================
EnrichedTriang::~EnrichedTriang()
//===========================================================================
{
}

//===========================================================================
void EnrichedTriang::evalNode(vector<pair<Point, double> >& nodeval)
//===========================================================================
{
  for (size_t ki=0; ki<node_.size(); ++ki)
    nodeval.push_back(node_[ki]->getPosAndVal());
}


//===========================================================================
void EnrichedTriang::evalEdge(int level, 
			      vector<pair<Point, double> >& edgeval)
//===========================================================================
{
  int nmb = (level == 1) ? 1 : 3;  // To be extended
  double del = 1.0/(double)(nmb+1);

  for (size_t kj=0; kj<edge_.size(); ++kj)
    {
      for (int ki=1; ki<=nmb; ++ki)
	edgeval.push_back(edge_[kj]->getPosAndVal(ki*del));
    }
}

//===========================================================================
void EnrichedTriang::evalTriang(int level, 
				vector<pair<Point, double> >& triangval)
//===========================================================================
{
  if (level <= 1)
    {
      for (size_t ki=0; ki<triang_.size(); ++ki)
	triangval.push_back(triang_[ki]->getMid());
    }
  else 
    {
      for (size_t ki=0; ki<triang_.size(); ++ki)
	{
	  vector<pair<Point, double> > currval;
	  currval = triang_[ki]->getPosAndVal(level);
	  triangval.insert(triangval.end(), currval.begin(), currval.end());
	}
    }
      
}

//===========================================================================
void EnrichedTriang::evaluate(int level, vector<pair<Point, double> >& result)
//===========================================================================
{
  vector<pair<Point, double> > nodeval;
  // evalNode(nodeval);
  // result.insert(result.end(), nodeval.begin(), nodeval.end());

  if (level > 1)
    {
      vector<pair<Point, double> > edgeval;
      evalEdge(level-1, edgeval);
      result.insert(result.end(), edgeval.begin(), edgeval.end());
    }

  vector<pair<Point, double> > triangmid;
  evalTriang(level-1, triangmid);
  result.insert(result.end(), triangmid.begin(), triangmid.end());
}

//===========================================================================
void Node::addTriangle(Triang* triangle)
//===========================================================================
{
  triang_.push_back(triangle);
}

//===========================================================================
double Node::getVal()
//===========================================================================
{
  if (!val_set_)
    {
      computeVal();
    }
  return val_;
}

//===========================================================================
pair<Point,double> Node::getPosAndVal()
//===========================================================================
{
  double val = getVal();
  return make_pair(pos_, val);
}

//===========================================================================
void Node::computeVal()
//===========================================================================
{
  val_ = 0.0;
  for (size_t ki=0; ki<triang_.size(); ++ki)
    val_ += triang_[ki]->getVal();
  if (triang_.size() > 0)
    val_ /= (double)triang_.size();
  val_set_ = true;
}

//===========================================================================
pair<Point,double> Edge::getPosAndVal(double fac)
//===========================================================================
{
  double val = getVal(fac);
  Point pos = fac*node_[0]->getPos() + (1.0-fac)*node_[1]->getPos();
  return make_pair(pos, val);
}

//===========================================================================
vector<shared_ptr<Edge> > Node::defineEdges()
//===========================================================================
{
  vector<shared_ptr<Edge> > new_edges;
  for (size_t ki=0; ki<triang_.size(); ++ki)
    {
      for (int ka=0; ka<3; ++ka)
	{
	  Node * curr = triang_[ki]->node_[ka];
	  if (curr != this)
	    {
	      // Check if the associated triangle edge is defined
	      size_t kj;
	      for (kj=0; kj<edge_.size(); ++kj)
		if ((edge_[kj]->node_[0] == curr && 
		     edge_[kj]->node_[1] == this) ||
		    (edge_[kj]->node_[0] == this && 
		     edge_[kj]->node_[1] == curr))
		  break;
	      
	      if (kj == edge_.size())
		{
		  shared_ptr<Edge> edg(new Edge(this, curr));
		  new_edges.push_back(edg);
		  edge_.push_back(edg.get());
		}
	    }
	}
    }
  return new_edges;
}

//===========================================================================
Triang::Triang(Node *n1, Node *n2, Node *n3, double val)
//===========================================================================
{
  node_[0] = n1;
  node_[1] = n2;
  node_[2] = n3;
  val_ = val;

  n1->addTriangle(this);
  n2->addTriangle(this);
  n3->addTriangle(this);
}

//===========================================================================
pair<Point,double> Triang::getMid()
//===========================================================================
{
  Point mid = node_[0]->getPos() + node_[1]->getPos() + node_[2]->getPos();
  mid /= 3.0;

  return make_pair(mid, val_);
}

//===========================================================================
vector<pair<Point,double> > Triang::getPosAndVal(int level)
//===========================================================================
{
  vector<pair<Point,double> > result;
  if (level >= 1)
    {
      Point mid = node_[0]->getPos() + node_[1]->getPos() + node_[2]->getPos();
      mid /= 3.0;
      result.push_back(make_pair(mid, val_));
      if (level > 1)
	{
	  vector<pair<Point,double> > nodeval(3);
	  for (int ki=0; ki<3; ++ki)
	    nodeval[ki] = node_[ki]->getPosAndVal();
			      
	  for (int ki=0; ki<3; ++ki)
	    {
	      Point pos = 0.5*(mid+nodeval[ki].first);
	      double val = 0.5*(val_ + nodeval[ki].second);
	      result.push_back(make_pair(pos, val));
	    }
	  for (int ki=0; ki<3; ++ki)
	    {
	      int kj = (ki+1)%3;
	      Point pos1 = 0.5*(nodeval[ki].first+nodeval[kj].first);
	      double val1 = 0.5*(nodeval[ki].second+nodeval[kj].second);
	      Point pos = 0.5*(mid+pos1);
	      double val = 0.5*(val_+val1);
	      result.push_back(make_pair(pos, val));
	    }
	}
    }
  return result;
}

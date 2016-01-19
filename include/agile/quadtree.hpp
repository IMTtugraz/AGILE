// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.

// $Id: quadtree.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef _AGILE_QUADTREE_HPP
#define _AGILE_QUADTREE_HPP

#include <vector>
#include <list>
#include <iostream>

#define QUADTREE_DEFAULT_MAX_DEPTH 10

#define QUADTREE_NUM_CHILDREN 4

#define QUADTREE_QUAD_POSITION_SOUTH_WEST 0
#define QUADTREE_QUAD_POSITION_SOUTH_EAST 1
#define QUADTREE_QUAD_POSITION_NORTH_WEST 2
#define QUADTREE_QUAD_POSITION_NORTH_EAST 3

namespace agile
{

  typedef enum HILBERT_SYMBOL_ENUM {
    HILBERT_SYMBOL_DOWN_RIGHT_UP,
    HILBERT_SYMBOL_UP_LEFT_DOWN,
    HILBERT_SYMBOL_RIGHT_DOWN_LEFT,
    HILBERT_SYMBOL_LEFT_UP_RIGHT
  } HILBERT_SYMBOL;

  typedef enum QUAD_STATE_ENUM {
    QUAD_EMPTY,
    QUAD_LEAF,
    QUAD_NODE
  } QUAD_STATE;


  class vec2_t {
    public:
      float x, y;

      vec2_t() : x(0.0f), y(0.0f) { }

      vec2_t(float tx, float ty) : x(tx), y(ty) { }

      vec2_t operator+(const vec2_t &w) const {
          vec2_t v;
          v.x = x + w.x;
          v.y = y + w.y;
          return v;
      }
      vec2_t operator-(const vec2_t &w) const {
          vec2_t v;
          v.x = x - w.x;
          v.y = y - w.y;
          return v;
      }
      vec2_t operator*(float c) {
          vec2_t v;
          v.x = x * c;
          v.y = y * c;
          return v;
      }/*
      vec2_t operator*(vec2_t &w, float c) {
          vec2_t v;
          v.x = c * w.x;
          v.y = c * w.y;
          return v;
      }*/
      vec2_t& operator+=(const vec2_t &w) {
          x += w.x;
          y += w.y;
          return *this;
      }
      vec2_t& operator-=(const vec2_t &w) {
          x -= w.x;
          y -= w.y;
          return *this;
      }
      vec2_t compmult(const vec2_t &w) {
        vec2_t v;
        v.x = x * w.x;
        v.y = y * w.y;
        return v;
      }
      friend std::ostream& operator<<(std::ostream& os, const vec2_t& v) {
        return os << "(" << v.x << ", " << v.y << ")";
      }
  };


  class QuadtreeItem {
    public:
      QuadtreeItem() : orgIndex(0) {}
      QuadtreeItem(float x, float y, unsigned oIdx) : pos(x, y), orgIndex(oIdx) {}
      vec2_t pos;
      unsigned orgIndex;

      friend std::ostream& operator<<(std::ostream& os, const QuadtreeItem& item) {
        return os << "qItem[" << item.orgIndex << "]";
      }
      friend std::ostream& operator<<(std::ostream& os, const QuadtreeItem *item) {
        return os << "qItem[" << item->orgIndex << "]";
      }
  };


  //  -------
  // | 2 | 3 |
  //  -------    Quad orientation
  // | 0 | 1 |
  //  -------
  class Quadtree
  {
    typedef std::list<QuadtreeItem*> QuadtreeDataListType;

    private:

      vec2_t m_bounds[2];

      unsigned m_max_depth;

      unsigned m_current_depth;

      QUAD_STATE m_state;

      Quadtree* m_quads[QUADTREE_NUM_CHILDREN];

      QuadtreeDataListType m_data;

      // assignment parameter shouldn't
      // be accessable
      Quadtree *operator =(Quadtree &) {
          return NULL;
      };

      // default constructor is not accessable from outside class
      Quadtree() {}

      // internal constructor for building sub-quadtrees
      Quadtree(vec2_t lower_bound, vec2_t upper_bound, unsigned max_depth,
              unsigned current_depth)
      {
        init(lower_bound, upper_bound, max_depth, current_depth);
      }

      void init(vec2_t lower_bound, vec2_t upper_bound, unsigned max_depth,
                unsigned current_depth);

      vec2_t getLowerSubQuadBound(unsigned sub_quad_idx) const;
      vec2_t getUpperSubQuadBound(unsigned sub_quad_idx) const;

      bool isQuadtreeItemEqual(QuadtreeItem *itemOne, QuadtreeItem *itemTwo);

      void getLinearHilbertOrderList(std::vector<QuadtreeItem*> &item_list,
                                     HILBERT_SYMBOL mode);

    public:
      Quadtree(vec2_t bounds[2]) {
        init(bounds[0], bounds[1], QUADTREE_DEFAULT_MAX_DEPTH, 0);
      };
      Quadtree(vec2_t bounds[2], unsigned max_depth) {
        init(bounds[0], bounds[1], max_depth, 0);
      };
      Quadtree(vec2_t lower_bound, vec2_t upper_bound) {
        init(lower_bound, upper_bound, QUADTREE_DEFAULT_MAX_DEPTH, 0);
      }
      Quadtree(vec2_t lower_bound, vec2_t upper_bound, unsigned max_depth) {
        init(lower_bound, upper_bound, max_depth, 0);
      }

      virtual ~Quadtree();

      bool contains(QuadtreeItem *item) const;
      bool addItem(QuadtreeItem *item);
      void getLinearZOrderList(std::vector<QuadtreeItem*> &item_list);
      void getLinearHilbertOrderList(std::vector<QuadtreeItem*> &item_list)
      {
        return getLinearHilbertOrderList(item_list,
                                         HILBERT_SYMBOL_DOWN_RIGHT_UP);
      }
  };

} // namespace agile

#endif // _AGILE_QUADTREE_HPP

// End of $Id: quadtree.hpp 476 2011-06-16 08:54:14Z freiberger $

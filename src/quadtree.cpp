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

// $Id: quadtree.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/quadtree.hpp"

namespace agile
{


  void Quadtree::init(vec2_t lower_bound, vec2_t upper_bound, unsigned max_depth,
                      unsigned current_depth)
  {
    m_bounds[0] = lower_bound;
    m_bounds[1] = upper_bound;
    for (unsigned i = QUADTREE_NUM_CHILDREN; i--; )
      m_quads[i] = NULL;
    m_state = QUAD_EMPTY;
    m_max_depth = max_depth;
    m_current_depth = current_depth;
    m_data.clear();
  }

  Quadtree::~Quadtree()
  {
    for (unsigned i = QUADTREE_NUM_CHILDREN; i--; )
    {
      if (m_quads[i] != NULL)
      {
        delete m_quads[i];
        m_quads[i] = NULL;
      }
    }
  }

  vec2_t Quadtree::getLowerSubQuadBound(unsigned sub_quad_idx) const
  {
    vec2_t half = (m_bounds[1] - m_bounds[0]) * 0.5f;
    vec2_t b = vec2_t(sub_quad_idx & 0x01, (sub_quad_idx & 0x02) >> 1);
    return m_bounds[0] + b.compmult(half);
  }

  vec2_t Quadtree::getUpperSubQuadBound(unsigned sub_quad_idx) const
  {
    vec2_t half = (m_bounds[1] - m_bounds[0]) * 0.5f;
    vec2_t b = vec2_t((sub_quad_idx & 0x01) ^ 0x01, ((sub_quad_idx & 0x02) >> 1) ^ 0x01);
    return m_bounds[1] - b.compmult(half);
  }

  bool Quadtree::contains(QuadtreeItem *item) const
  {
    // check if item is within bounds
    return (item->pos.x >= m_bounds[0].x && item->pos.y >= m_bounds[0].y &&
            item->pos.x < m_bounds[1].x && item->pos.y < m_bounds[1].y);
  }

  void Quadtree::getLinearZOrderList(std::vector<QuadtreeItem*> &item_list)
  {
    if (m_state != QUAD_EMPTY)
    {
      if (m_state == QUAD_LEAF)
      {
        for (QuadtreeDataListType::iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
          item_list.push_back(*iter);
      }
      else if (m_state == QUAD_NODE)
      {
        m_quads[QUADTREE_QUAD_POSITION_NORTH_WEST]->getLinearZOrderList(item_list);
        m_quads[QUADTREE_QUAD_POSITION_NORTH_EAST]->getLinearZOrderList(item_list);
        m_quads[QUADTREE_QUAD_POSITION_SOUTH_WEST]->getLinearZOrderList(item_list);
        m_quads[QUADTREE_QUAD_POSITION_SOUTH_EAST]->getLinearZOrderList(item_list);
      }
    }
  }

  void Quadtree::getLinearHilbertOrderList(std::vector<QuadtreeItem*> &item_list, HILBERT_SYMBOL mode)
  {
    if (m_state != QUAD_EMPTY)
    {
      if (m_state == QUAD_LEAF)
      {
        for (QuadtreeDataListType::iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
          item_list.push_back(*iter);
      }
      else if (m_state == QUAD_NODE)
      {
        switch (mode)
        {
          case HILBERT_SYMBOL_DOWN_RIGHT_UP:
            m_quads[QUADTREE_QUAD_POSITION_NORTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_RIGHT_DOWN_LEFT);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_DOWN_RIGHT_UP);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_DOWN_RIGHT_UP);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_LEFT_UP_RIGHT);
            break;
          case HILBERT_SYMBOL_UP_LEFT_DOWN:
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_LEFT_UP_RIGHT);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_UP_LEFT_DOWN);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_UP_LEFT_DOWN);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_RIGHT_DOWN_LEFT);
            break;
          case HILBERT_SYMBOL_RIGHT_DOWN_LEFT:
            m_quads[QUADTREE_QUAD_POSITION_NORTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_DOWN_RIGHT_UP);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_RIGHT_DOWN_LEFT);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_RIGHT_DOWN_LEFT);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_UP_LEFT_DOWN);
            break;
          case HILBERT_SYMBOL_LEFT_UP_RIGHT:
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_UP_LEFT_DOWN);
            m_quads[QUADTREE_QUAD_POSITION_SOUTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_LEFT_UP_RIGHT);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_WEST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_LEFT_UP_RIGHT);
            m_quads[QUADTREE_QUAD_POSITION_NORTH_EAST]->getLinearHilbertOrderList(item_list, HILBERT_SYMBOL_DOWN_RIGHT_UP);
            break;
        }
      }
    }
  }

  bool Quadtree::isQuadtreeItemEqual(QuadtreeItem *itemOne, QuadtreeItem *itemTwo)
  {
    return itemOne->pos.x == itemTwo->pos.x && itemOne->pos.y == itemTwo->pos.y;
  }

  bool Quadtree::addItem(QuadtreeItem *item)
  {
    // first check if item is within bounds
    if (contains(item))
    {
      if (m_state == QUAD_EMPTY || m_current_depth == m_max_depth
          || (m_data.size() > 0 && isQuadtreeItemEqual(*(m_data.begin()), item)))
      {
        // QUAD_EMPTY => QUAD_LEAF
        m_data.push_back(item);
        m_state = QUAD_LEAF;
      }
      else
      {
        if (m_state == QUAD_LEAF)
        {
          // QUAD_LEAF => QUAD_NODE
          for (unsigned i = QUADTREE_NUM_CHILDREN; i--; )
          {
            if (m_quads[i] == NULL)
            {
              m_quads[i] = new Quadtree(getLowerSubQuadBound(i),
                                        getUpperSubQuadBound(i),
                                        m_max_depth,
                                        m_current_depth + 1);
            }

            // move all data down to child nodes
            for (QuadtreeDataListType::iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
            {
              if (m_quads[i]->contains(*iter)) {
                m_quads[i]->addItem(*iter);
              }
            }
          }
          // clear data in node, already moved to child-quads above
          m_data.clear();
          m_state = QUAD_NODE;
        }

        // add item to child
        for (unsigned i = QUADTREE_NUM_CHILDREN; i--; )
        {
          if (m_quads[i]->contains(item)) {
            m_quads[i]->addItem(item);
            break;
          }
        }
      }
      return true;
    }
    return false;
  }

} // namespace agile

// End of $Id: quadtree.cpp 476 2011-06-16 08:54:14Z freiberger $.

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

// $Id: io_utils.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/io/io_utils.hpp"

#include <algorithm>
#include <sstream>
#include <iterator>

namespace agile
{
  void FileHeader::read(std::ifstream& ifstr)
  {
    AGILE_ASSERT(ifstr.is_open(), StandardException::ExceptionIO());

    std::string line;
    ifstr >> line;
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    if (line.compare("version: 1.0") != 0)
      throw StandardException::ExceptionIO();
    version = "1.0";

    ifstr >> line;
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    if (line.find("type: ") != 0)
      throw StandardException::ExceptionIO();
    file_type = line.substr(6);

    size.clear();
    ifstr >> line;
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    if (line.find("size: ") != 0)
      throw StandardException::ExceptionIO();
    std::istringstream sstr(line.substr(6));
    std::istream_iterator<unsigned> iter(sstr), end_iter;
    while (iter != end_iter)
      size.push_back(*iter++);

    ifstr >> comment;
    line = comment;
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    if (line.find("comment: ") != 0)
      throw StandardException::ExceptionIO();
    comment = comment.substr(9);

    // split the file type using a space as delimiter
    split_file_type.clear();
    std::string::size_type start_pos = file_type.find_first_not_of(" ", 0);
    std::string::size_type end_pos = file_type.find_first_of(" ", start_pos);
    while (start_pos != std::string::npos)
    {
      if (end_pos == std::string::npos)
        split_file_type.push_back(file_type.substr(start_pos));
      else
        split_file_type.push_back(file_type.substr(
                                    start_pos, end_pos - start_pos));
      start_pos = file_type.find_first_not_of(" ", end_pos);
      end_pos = file_type.find_first_of(" ", start_pos);
    }

    // determine the size of a type
    type_size = 0;
    if (split_file_type.empty())
      return;
    if (split_file_type.back() == type_to_file_string<unsigned>::get())
      type_size = sizeof(unsigned);
    else if (split_file_type.back() == type_to_file_string<float>::get())
      type_size = sizeof(float);
    else if (split_file_type.back() == type_to_file_string<double>::get())
      type_size = sizeof(double);

    // if the type is complex, the size has to be doubled
    if (split_file_type.size() > 1
        && split_file_type.at(split_file_type.size() - 2) == "complex")
      type_size *= 2;
  }

  bool FileHeader::checkType(const std::string& check_type) const
  {
    std::istringstream check_sstr(check_type);
    std::istream_iterator<std::string> check_iter(check_sstr), check_iter_end;
    std::istringstream type_sstr(file_type);
    std::istream_iterator<std::string> type_iter(type_sstr), type_iter_end;

    while (check_iter != check_iter_end)
    {
      if (type_iter == type_iter_end
          || check_iter->compare(*type_iter) != 0)
        return false;
      ++check_iter;
      ++type_iter;
    }
    return true;
  }

  void writeHeader(std::ofstream& ofstr, const std::string& type,
                   const std::vector<unsigned>& size,
                   const std::string& comment)
  {
    // make sure the stream is open
    AGILE_ASSERT(ofstr.is_open(), StandardException::ExceptionIO());

    // write version, type, size and comment
    ofstr << "Version: 1.0" << std::endl;
    ofstr << "Type: " << type << std::endl;
    ofstr << "Size: ";
    for (unsigned counter = 0; counter < size.size(); ++counter)
      ofstr << size[counter] << " ";
    ofstr << std::endl;
    ofstr << "Comment: " << comment << std::endl;
  }

} // namespace agile

// End of $Id: io_utils.cpp 476 2011-06-16 08:54:14Z freiberger $.

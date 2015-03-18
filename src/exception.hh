/*
taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

Copyright (C) 2010 Johannes Dröge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef exception_hh_
#define exception_hh_

#include <boost/exception/all.hpp>
#include <exception>
#include <string>

// tags
typedef boost::error_info<struct exception_tag_seqid, const std::string> general_info;
typedef boost::error_info<struct exception_tag_seqid, const std::string> file_info;
typedef boost::error_info<struct exception_tag_seqid, const std::string> seqid_info;
typedef boost::error_info<struct exception_tag_seqid, const TaxonID> taxid_info;
typedef boost::error_info<struct exception_tag_seqid, const uint> line_info;


// exception classes
class Exception : public boost::exception, public std::exception {};

class SequenceNotFound : public Exception {
    const char *what() const noexcept { return "could not find sequence identifier"; }
};

class TaxonNotFound : public Exception {
    const char *what() const noexcept { return "could not find taxon in taxonomy"; }
};

class TaxonIDNotFound : public Exception {
    const char *what() const noexcept { return "could not find taxon identifier in mapping"; }
};

class FileNotFound : public Exception {
    const char *what() const noexcept { return "could not find file"; }
};

class FileError : public Exception {
    const char *what() const noexcept { return "could not access file"; }
};

class ParsingError : public Exception {
  const char *what() const noexcept { return "could not parse record"; }
};

class EOFError : public Exception {
  const char *what() const noexcept { return "could not read: end of file"; }
};

class GeneralError : public Exception {
  const char *what() const noexcept { return "general error"; }
};

#endif // exception_hh_

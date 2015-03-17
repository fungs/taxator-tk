#include <boost/exception/exception.hpp>
#include <exception>
#include <string>

// tags
typedef boost::error_info<struct exception_tag_seqid, std::string> seqid_info;

// exception classes
class SequenceNotFound : public boost::exception, public std::exception {
    const char *what() const noexcept { return "could not find sequence"; }
};

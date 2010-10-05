#ifndef PLY_IO_OPERATORS_HPP_INCLUDED
#define PLY_IO_OPERATORS_HPP_INCLUDED

#include <istream>
#include <limits>
#include <ostream>

namespace ply {

namespace io_operators {

inline std::istream& operator>>(std::istream& istream, int8 &value)
{
  int16 tmp;
  if (istream >> tmp) {
    if (tmp <= std::numeric_limits<int8>::max()) {
      value = static_cast<int8>(tmp);
    }
    else {
      istream.setstate(std::ios_base::failbit);
    }
  }
  return istream;
}

inline std::istream& operator>>(std::istream& istream, uint8 &value)
{
  uint16 tmp;
  if (istream >> tmp) {
    if (tmp <= std::numeric_limits<uint8>::max()) {
      value = static_cast<uint8>(tmp);
    }
    else {
      istream.setstate(std::ios_base::failbit);
    }
  }
  return istream;
}

inline std::ostream& operator<<(std::ostream& ostream, int8 value)
{
  return ostream << static_cast<int16>(value);
}

inline std::ostream& operator<<(std::ostream& ostream, uint8 value)
{
  return ostream << static_cast<uint16>(value);
}

} // namespace io_operators

} // namespace ply

#endif // PLY_IO_OPERATORS_HPP_INCLUDED

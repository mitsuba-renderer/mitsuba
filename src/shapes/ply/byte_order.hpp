#ifndef PLY_BYTE_ORDER_HPP_INCLUDED
#define PLY_BYTE_ORDER_HPP_INCLUDED

namespace ply {

#if defined(PLY_BIG_ENDIAN) || defined(PLY_LITTLE_ENDIAN)
#  error
#endif

#if (defined(__powerpc) || defined(__powerpc__) || defined(__POWERPC__) || defined(__ppc__) || defined(_M_PPC) || defined(__ARCH_PPC))
#  define PLY_BIG_ENDIAN
#elif (defined(i386) || defined(__i386__) || defined(__i386) || defined(_M_IX86) || defined(_X86_) || defined(__THW_INTEL__) || defined(__I86__) || defined(__INTEL__)) \
   || (defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(_M_X64)) \
   || defined(__ARMEL__)
#  define PLY_LITTLE_ENDIAN
#else
#  error
#endif

enum byte_order
{
  little_endian_byte_order = 0,
  big_endian_byte_order = 1,
#if defined(PLY_BIG_ENDIAN)
  host_byte_order = big_endian_byte_order,
#elif defined(PLY_LITTLE_ENDIAN)
  host_byte_order = little_endian_byte_order,
#else
#  error
#endif
  network_byte_order = big_endian_byte_order
};

#undef PLY_BIG_ENDIAN
#undef PLY_LITTLE_ENDIAN

template <std::size_t N>
void swap_byte_order(char* bytes);

template <>
inline void swap_byte_order<1>(char* bytes)
{
}

template <>
inline void swap_byte_order<2>(char* bytes)
{
  std::swap(bytes[0], bytes[1]);
}

template <>
inline void swap_byte_order<4>(char* bytes)
{
  std::swap(bytes[0], bytes[3]);
  std::swap(bytes[1], bytes[2]);
}

template <>
inline void swap_byte_order<8>(char* bytes)
{
  std::swap(bytes[0], bytes[7]);
  std::swap(bytes[1], bytes[6]);
  std::swap(bytes[2], bytes[5]);
  std::swap(bytes[3], bytes[4]);
}

template <typename T>
void swap_byte_order(T& value)
{
  swap_byte_order<sizeof(T)>(reinterpret_cast<char*>(&value));
}

} // namespace ply

#endif // PLY_BYTE_ORDER_HPP_INCLUDED

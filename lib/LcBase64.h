/**
 * @file	LcBase64.h
 *
 * @brief	Class to provide encoding/decoding from base64 format.
 */

#ifndef LC_BASE64_H
#define LC_BASE64_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

// libb64 includes
#include <b64/encode.h>
#include <b64/decode.h>

namespace Lucee
{
/**
 * Class to encode/decode data with base64 encoding.
 */
  template <typename T>
  class Base64
  {
    public:
/**
 * Encode the data in vector as a base64 string. Length of encoded data is returned.
 *
 * @param data Data to encode.
 * @param b64out Pointer to encoded data. This must be an empty buffer on input.
 * @return Length of encoded data.
 */
      int encode(const std::vector<T>& data, char *b64out);

/**
 * Decode the data in base64 string into a vector.
 *
 * @param len Lenght of data buffer (b64in).
 * @param b64in Pointer to data to decode.
 * @param data Data vector with decoded data. This will be cleared on entry.
 */
      void decode(int len, const char *b64in, std::vector<T>& data);
  };
}

#endif // LC_BASE64_H

/**
 * @file	LcBase64.cpp
 *
 * @brief	Class to provide encoding/decoding from base64 format.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBase64.h>


namespace Lucee
{
  template <typename T>
  int
  Base64<T>::encode(const std::vector<T>& data, char *b64out)
  {
    int len = sizeof(T)*data.size(); // size in bytes
    base64::encoder ec;
    const T *dptr = &data[0];
    return ec.encode((const char*) dptr, len, b64out);
  }

  template <typename T>
  void
  Base64<T>::decode(int len, const char *b64in, std::vector<T>& data)
  {
    base64::decoder dc;
    char *dout;
    int outlen = dc.decode(b64in, len, dout);
// now pop stuff into vector
    data.clear();
    T *dptr = (T *) dout;
    unsigned nvals = outlen/sizeof(T);
    data.resize(nvals);
    for (unsigned i=0; i<nvals; ++i)
      data[i] = dptr[i];
    delete [] dptr;
  }

// instantiations
  template class Lucee::Base64<unsigned>;
  template class Lucee::Base64<int>;
  template class Lucee::Base64<float>;
  template class Lucee::Base64<double>;
}

/**
 * @file	LcArrayIo.h
 *
 * @brief	Functions for I/O for various Lucee arrays.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_ARRAY_IO_H
#define LC_ARRAY_IO_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>
#include <LcIoBase.h>
#include <LcMatrix.h>
#include <LcVector.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Write a Lucee::Vector to file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of data as it should appear in output.
 * @param vec Vector to write.
 * @return node to which data was written.
 */
  template <typename T>
  Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node, 
    const std::string& nm, Lucee::Vector<T>& vec)
  {
// construct sizes and shapes to write stuff out
    std::vector<size_t> dataSetSize(1), dataSetBeg(1), dataSetLen(1);
    dataSetSize[0] = vec.getLength();
    dataSetBeg[0] = 0;
    dataSetLen[0] = vec.getLength();
// make sure vector is contiguous
    Lucee::Vector<T> vecDup(vec);
    if (vec.isContiguous() == false)
      vecDup = vec.duplicate(); // not, so allocate fresh vector

// write it out
    return 
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &vecDup[0]);
  }

/**
 * Write a Lucee::Array to file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of data as it should appear in output.
 * @param arr Array to write.
 * @return node to which data was written.
 */
  template <unsigned NDIM, typename T>
  Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
    const std::string& nm, Lucee::Array<NDIM, T>& arr)
  {
// construct sizes and shapes to write stuff out
    std::vector<size_t> dataSetSize(NDIM), dataSetBeg(NDIM), dataSetLen(NDIM);

    for (unsigned i=0; i<NDIM; ++i)
    {
      dataSetSize[i] = arr.getShape(i);
      dataSetBeg[i] = 0;
      dataSetLen[i] = arr.getShape(i);
    }
// make sure vector is contiguous
    Lucee::Array<NDIM, T> arrDup(arr);
    if (arr.isContiguous() == false)
      arrDup = arr.duplicate(); // not, so allocate fresh array

// write it out
    return 
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &arrDup.first());
  }
}

#endif // LC_ARRAY_IO_H

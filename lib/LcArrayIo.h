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
// write it out
    return 
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &vec[0]);
  }
}

#endif // LC_ARRAY_IO_H

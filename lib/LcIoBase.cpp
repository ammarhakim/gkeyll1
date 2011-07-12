/**
 * @file	LcIoBase.cpp
 *
 * @brief	Base class for I/O using hierachical file system.
 */

// lib includes
#include <LcIoBase.h>
#include <LcIoTmpl.h>

// std includes
#include <sstream>

namespace Lucee
{
  IoBase::IoBase() 
  {
  }

  IoBase::~IoBase() 
  {
  }

  void IoBase::addOpenFile(Lucee::IoNodeType node) 
  {
    openFiles.push_back(node);
  }

  void IoBase::removeOpenFile(Lucee::IoNodeType node) 
  {
    std::vector<Lucee::IoNodeType>::iterator itr;
    for (itr = openFiles.begin(); itr != openFiles.end(); ++itr) 
    {
      if (**itr == *node) 
      {
        openFiles.erase(itr);
        break;
      }
    }
  }

  void IoBase::closeOpenFiles() 
  {
    std::vector<Lucee::IoNodeType>::iterator itr;
    for (itr = openFiles.begin(); itr != openFiles.end(); ++itr) 
    {
      delete *itr;
    }
    openFiles.clear();
  }
}

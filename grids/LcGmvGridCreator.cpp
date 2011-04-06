/**
 * @file	LcGmvGridCreator.cpp
 *
 * @brief	Creator for grids defined in GMV format.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

// lucee includes
#include <LcExcept.h>
#include <LcGmvGridCreator.h>

namespace Lucee
{
  static
  void
  splitLine(std::vector<std::string>& strVec, const std::string& str, const std::string& delim)
  {
    strVec.clear();
    size_t  start = 0, end = 0;
    while ( end != std::string::npos )
    {
      end = str.find( delim, start );
      strVec.push_back( str.substr( start,
          (end == std::string::npos) ? std::string::npos : end - start ) );
      start = (   ( end > (std::string::npos - delim.size()) )
        ?  std::string::npos  :  end + delim.size()    );
    }
  }

  static
  void
  convertToIntVec(int s, int e, 
    const std::vector<std::string>& tokens, std::vector<int>& ivec)
  {
    ivec.clear();
    for (unsigned i=s; i<e; ++i)
      ivec.push_back(std::atoi(tokens[i].c_str()));
  }

  template <typename REAL>
  GmvGridCreator<REAL>::GmvGridCreator(unsigned ndim, std::istream& gmvStrm)
    : UnstructGridCreator<REAL>(ndim), ndim(ndim)
  {
    readFromGmvFile(gmvStrm);
  }

  template <typename REAL>
  void
  GmvGridCreator<REAL>::readFromGmvFile(std::istream& gmvStrm)
  {
    std::vector<std::string> tokens;
    std::string line;
    char* linePtr;

// line 1: identification
    std::getline(gmvStrm, line);
    splitLine(tokens, line, " ");
    if (tokens[0] != "gmvinput")
    {
      Lucee::Except lce("GmvGridCreator::readFromGmvFile: File stream is not a GMV file");
      throw lce;
    }
// line 2: number of nodes
    std::getline(gmvStrm, line);
    splitLine(tokens, line, " ");
    if (tokens[0] != "nodes")
    {
      Lucee::Except lce("GmvGridCreator::readFromGmvFile: Second line must specify number of nodes");
      throw lce;
    }
    int nnodes = std::atoi(tokens[1].c_str());
    this->setNumVertices(nnodes);

// read X coordinates
    std::getline(gmvStrm, line);
    linePtr = &line[0];
    for (unsigned i=0; i<nnodes; ++i)
    {
      double x = std::strtod(linePtr, &linePtr);
      this->setVertexXCoord(i, x);
    }

// read Y coordinates
    std::getline(gmvStrm, line);
    linePtr = &line[0];
    for (unsigned i=0; i<nnodes; ++i)
    {
      double y = std::strtod(linePtr, &linePtr);
      this->setVertexYCoord(i, y);
    }

    if (ndim==3)
    { 
// read Z coordinates
      std::getline(gmvStrm, line);
      linePtr = &line[0];
      for (unsigned i=0; i<nnodes; ++i)
      {
        double z = std::strtod(linePtr, &linePtr);
        this->setVertexZCoord(i, z);
      }
    }

// line with number of cells
    std::getline(gmvStrm, line);
    splitLine(tokens, line, " ");
    if (tokens[0] != "cells")
    {
      Lucee::Except lce("GmvGridCreator::readFromGmvFile: After vertices must specify number of cells");
      throw lce;
    }
    int ncells = std::atoi(tokens[1].c_str());
    this->setNumCells(ncells);

    std::vector<int> ivec;
// read lines with cell connectivities
    for (unsigned i=0; i<ncells; ++i)
    {
      std::getline(gmvStrm, line);
      splitLine(tokens, line, " ");
      convertToIntVec(2, tokens.size(), tokens, ivec);

// check cell type
      if (tokens[0] == "tet")
      {
        this->addTet(ivec[0]-1, ivec[1]-1, ivec[2]-1, ivec[3]-1);
      }
      else if (tokens[0] == "hex")
      {
      }
      else if (tokens[0] == "tri")
      {
        this->addTri(ivec[0]-1, ivec[1]-1, ivec[2]-1);
      }
      else if (tokens[0] == "quad")
      {
        this->addQuad(ivec[0]-1, ivec[1]-1, ivec[2]-1, ivec[3]-1);
      }        
    }
  }

// instantiations
  template class GmvGridCreator<float>;
  template class GmvGridCreator<double>;
}

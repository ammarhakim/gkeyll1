/**
 * @file	LcUnstructGrid.h
 *
 * @brief	Unstructured grid class.
 */
#ifndef LC_UNSTRUCT_CELL_ENUM_H
#define LC_UNSTRUCT_CELL_ENUM_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/** Enum for all supported unstructured cell types */
  enum UnstructCellEnum
  {
    LC_LINE, // line
    LC_TRI, // triangle
    LC_QUAD, // quadrilateral
    LC_TET, // tetrahedron
    LC_HEX // hexahedron
  };
}

#endif // LC_UNSTRUCT_CELL_ENUM_H

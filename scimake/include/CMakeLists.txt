######################################################################
#
# CMakeLists.txt for scimake/include
#
# $Id: CMakeLists.txt 657 2014-10-25 17:30:08Z jrobcary $
#
# Copyright &copy; 2005-2013 Tech-X Corporation, Boulder, CO
# See LICENSE file for conditions of use.
#
######################################################################

if (FALSE)

if (NOT DEFINED SCI_CREATE_STD_INCLUDES)
  set(SCI_CREATE_STD_INCLUDES TRUE)
endif ()
if (NOT DEFINED SCI_PACKAGE_PREFIX)
  set(SCI_PACKAGE_PREFIX sci)
endif ()
string(TOUPPER ${SCI_PACKAGE_PREFIX} SCI_PACKAGE_PREFIX_UC)
message(STATUS "SCI_PACKAGE_PREFIX_UC = ${SCI_PACKAGE_PREFIX_UC}")
if (NOT DEFINED SCI_STD_TYPES_INSTALL)
  set(SCI_STD_TYPES_INSTALL TRUE)
endif ()

macro(CreateStdTypeFile type)
  string(TOUPPER ${type} typeuc)
  if (HAVE_ANY_${typeuc})
    message(STATUS "creating ${SCI_PACKAGE_PREFIX}_${type}.")
    set(SCI_INCLUDE ${${typeuc}_INCLUDE})
    get_filename_component(SCI_INCLUDE_NAME ${${typeuc}_INCLUDE} NAME)
    string(REPLACE "." "_" SCI_INCLUDE_NAME ${SCI_INCLUDE_NAME})
    string(TOUPPER ${SCI_INCLUDE_NAME} SCI_INCLUDE_NAME_UC)
    set(SCI_TYPE ${type})
    set(SCI_NAMESPACE ${${typeuc}_NAMESPACE})
    configure_file(sci_type.in
      ${CMAKE_CURRENT_BINARY_DIR}/${SCI_PACKAGE_PREFIX}_${type}
    )
    if (SCI_STD_TYPES_INSTALL)
      install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${SCI_PACKAGE_PREFIX}_${type}
        DESTINATION include
      )
    endif ()
  endif ()
endmacro()

# Create c++11 files
if (SCI_CREATE_STD_INCLUDES)
  CreateStdTypeFile(shared_ptr)
  CreateStdTypeFile(array)
endif ()

# Create derived files
configure_file(sci_array_arith.in sci_array_arith)

endif ()


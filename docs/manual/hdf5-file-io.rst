Reading and writing HDF5 files
------------------------------

HDF5 is the native binary I/O format used in Lucee. The class
``Lucee::Hdf5FileIo`` provides methods to read and write HDF5
files. This class is a wrapper around the HDF5 library and makes
writing data of arbitrary types, in both serial and parallel, uniform
and easy.

``Lucee::Hdf5FileIo``: HDF5 file I/O
++++++++++++++++++++++++++++++++++++

.. class:: Hdf5FileIo

  This class provided methods to open/create HDF5 files, open/create
  groups, datasets and write/read attributes. All I/O methods work
  with ``int``, ``double``, ``std::string`` or ``std::vector`` of
  these. For example::

    Lucee::Hdf5FileIo io;

    // create file and open a group to write data
    Lucee::IoNode fh = io.openFile("ouput.h5", "w");
    Lucee::IoNode root = io.openGroup(fh, "/");
    Lucee::IoNode grp = io.createGroup(root, "lucee");

    // write some data to file
    

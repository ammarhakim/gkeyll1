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

  .. cfunction:: IoNode createFile(const std::string& fileName)
    :noindex:

    Create a new HDF5 file. A handle to the created file is returned.
  
  .. cfunction:: IoNode openFile(const std::string& fileName, const std::string& perms)
    :noindex:

    Open a new HDF5 file with specified permissions. The ``perms`` can
    be one of "r" or "rw". A handle to the opened file is
    returned. This should be used in subsequent calls to perform the
    actual I/O.

  .. cfunction:: void closeFile(IoNode fileNode)
    :noindex:

    Close the HDF5 file. The ``fileNode`` should be a node object
    returned by either a ``createFile`` or ``openFile``. An exception
    is thrown if ``fileNode`` is not returned by one of these methods.

  .. cfunction:: IoNode createGroup(IoNode node, const std::string& grp)
    :noindex:

    Creates a new group named ``grp``. The ``node`` object should be a
    file or group node. A handle to the created group is returned. An
    exception is thrown if the group already exists.

  .. cfunction:: IoNode openGroup(IoNode node, const std::string& dataName)
    :noindex:

    Open an existing group named ``grp``. The ``node`` object should
    be a file or a group node. A handle to the opened group is
    returned. An exception is thrown if the node does not exist.

  .. cfunction:: IoNode createDataSet(IoNode node, const std::string& dataName)
    :noindex:

    Create a new dataset named ``dataName`` for writing array
    data. The ``node`` should be a file or a group node. A handle to
    the created dataset is returned. An exception is thrown if the
    dataset already exists.

  .. cfunction:: IoNode openDataSet(IoNode node, const std::string& dataName)
    :noindex:

    Open an existing dataset named ``dataName`` for reading/writing
    array data. The ``node`` should be a file or a group node. A
    handle to the opened dataset is returned. An exception is thrown
    if the dataset does not exist.
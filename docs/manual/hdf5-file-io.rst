Reading and writing HDF5 files
------------------------------

HDF5 is the native binary I/O format used in Lucee. The class
``Lucee::Hdf5FileIo`` provides methods to read and write HDF5
files. This class is a wrapper around the HDF5 library and makes
writing data of arbitrary types, in both serial and parallel, uniform
and easy.
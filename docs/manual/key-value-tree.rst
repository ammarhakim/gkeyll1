Storing hierarchical data
-------------------------

Lucee provides the class ``Lucee::KeyVal`` to store data as (key,
value) pairs. Hierarchical data is represented by
``Lucee::KeyValTree`` which creates a tree of ``Lucee::KeyVal``
objects.

``Lucee::KeyVal``: Key, value pairs
+++++++++++++++++++++++++++++++++++

.. class:: KeyVal

  This class provides a container to store data as (key, value)
  pairs. The keys name the data as ``std::string``, while the values
  can be of type ``int``, ``double``, ``std::string`` or
  ``std::vector`` of these. For example::

    Lucee::KeyVal kv;
    kv.add("nx", 10);
    kv.add("ny", 20);
    kv.add("length", 12.5);

    std::vector<int> cells(2);
    cells[0] = 10; cells[1] = 20;
    kv.add("cells", cells);

    std::cout << kv.get<int>("nx") << std::endl; // 10
    std::cout << kv.get<double>("length") << std::endl; // 12.5

  .. cfunction:: bool add(const std::string& key, VALUETYPE value)
    :noindex:

    Add a new (key, value) pair to the object. Returns ``true`` if the
    data was added and ``false`` if the ``key`` already exists in the
    object.

  .. cfunction:: VALUETYPE& get(const std::string& key)
    :noindex:

    Return the value associated with ``key``. An exception is thrown
    if the ``key`` does not exists. This method must be called with
    the type of the expected value as a template parameter::

      nx = kv.get<int>("nx");
      length = kv.get<double>("length");

  .. cfunction:: bool has(const std::string& key)
    :noindex:

    Check if ``key`` exists in the object. This method must be called
    with the type of the queried value as a template parameter::

      bool hasNx = kv.has<int>("nx");
      bool hasLength = kv.has<double>("length");

  .. cfunction:: unsigned getNum()
    :noindex:

    Return number of (key, value) pairs stored for the specified
    type. This method must be called with the type of the queried
    value as a template parameter::

      unsigned numInts = kv.getNum<int>();
      unsigned numDbls = kv.getNum<double>();

  .. cfunction:: void setToFirst()
    :noindex:

    Set internal iterator of the object to point to the first (key,
    value) pair of the specified type. This is useful when used in
    conjuction with the ``getAndBump()`` method to iterate over all
    (key, value) pairs of a specified type. This method must be called
    with the type of the value as a template parameter::

      kv.setToFirst<int>(); // set to first int (key, value) pair

  .. cfunction:: std::pair<std::string, VALUETYPE> getAndBump()
    :noindex:

    Return the next (key, value) pair in object and set the internal
    iterator to the next pair. To print out all the ``int`` data
    stored in the set, for example::
     
      kv.setToFirst<int>();
      for (unsigned i=0; i<kv.getNum<int>(); ++i)
      {
        std::pair<std::string, int> p = kv.getAndBump();
	std::cout << p.first << " = " << p.second << std::endl;
      }
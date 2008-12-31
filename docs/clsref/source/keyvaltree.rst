.. -*- rst -*-

:mod:`Lucee::KeyValTree` --- Representing data hierarchically
=============================================================

.. highlight:: c++

.. module:: Lucee::KeyValTree
 :synopsis: Classes for representing data hierarchically

Reference for :class:`Lucee::KeyValTree`
----------------------------------------

Class defined in header file ``lckeyvaltree.h``.

.. class:: Lucee::KeyValTree(const std::string& nm)

  :param nm: Name of tree. Defaults to "KeyValTree".

  The class :class:`Lucee::KeyValTree` allows representing data as a
  tree. It is particularly useful to represent configuration data for
  a simulation. Each node in the tree can store arbitrary key-value
  pairs, where the keys are strings and values can be of any of
  type. Further, each node can have arbitrary number of children
  nodes, each accessible by thier names.

  .. method:: add(const std::string& key, T value)

    :param key: Key identifying of data
    :param value: Value of data
    :rtype: bool True if key-value pair added successfully, false otherwise

    Add a new key-value pair to the node. The method works with any
    data type ``T`` which supports a copy constructor.

  .. method:: get(const std::string& key)

    :param key: Key identifying data
    :rtype: T Value of data

    Get value of data with given key. This is a templated method and
    the expected type must be passed as a template parameter::

      int val = kvt.get<int>("myInt");

  .. method:: has(const std::string& key)

    :param key: Key identifying data
    :rtype: True if data exists in tree, false otherwise

    Check if tree has data with specified key. This is a templated
    method and the expected type must be passed as a template
    parameter::

      bool hasMyInt = kvt.has<int>("myInt");

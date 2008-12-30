.. -*- rst -*-

:mod:`Lucee::Logger` --- Logging messages to multiple streams
=============================================================

.. highlight:: c++

.. module:: Lucee::Logger
 :synopsis: Classes for logging messages to multiple streams.

The Lucee logging system provides an uniform interface to log messages
to multiple output streams. These streams may be directed to the
console or text files. Lucee loggers are non-intrusive, i.e.  class
and function *interfaces* do not need to be modified to use the
logger. A single log call puts the message in a controlled manner into
multiple streams at the same time.

Example usage
-------------

The logging interface is defined in the ``lclogger.h`` header
file. Loggers first need to be created and their verbosity set.::

  Lucee::Logger& logger = Lucee::Logger::create("myLogger");
  logger.setLevel("debug");

The :meth:`setLevel` method sets the verbosity level for the
logger. Messages with lower verbosity level will not be logger. The
levels, in increasing order of verbosity are, ``debug``, ``info``,
``warning``, ``error`` and ``critical``.

Once the logger is created we can attach it to a multiple handlers.::

  Lucee::StreamHandler strmHndlr(std::cout);
  strmHndlr.attachLogger(logger);

Now messages can be logged by getting log-streams::

  Lucee::LogStream dbgStrm = logger.getDebugStream();
  dbgStrm << "This is a debug message" << std::endl;

Reference
---------

.. class:: Lucee::Logger

  The class :class:`Lucee::Logger` provides methods to create new
  loggers and obtain streams for logging messages.

  .. staticmethod:: create(const std::string& nm) -> Lucee::Logger&

    Create a new logger named ``nm`` and return a reference to
    it. Once a logger with a given name is created it can be accessed
    by the :meth:`get` method.
  
  .. staticmethod:: get(const std::string& nm) -> Lucee::Logger&

    Get a reference to an already created logger. An exception is
    thrown if the logger does not exist.

  .. method:: getDebugStream() -> Lucee::LogStream
  .. method:: getInfoStream() -> Lucee::LogStream
  .. method:: getWarningStream() -> Lucee::LogStream

    Return a new stream with proper level to which messages can be
    logged.
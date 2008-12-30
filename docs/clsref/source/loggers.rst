.. -*- rst -*-

.. TODO:
.. 1) Indicate that StreamHandler and FileHandler are children of
.. LogRecordHandler.
.. 2) Ensure examples work.
.. 3) Document the enum LogMsgLevels.

:mod:`Lucee::Logger` --- Logging messages to multiple streams
=============================================================

.. highlight:: c++

.. module:: Lucee::Logger
 :synopsis: Classes for logging messages to multiple streams.

The Lucee logging system provides an uniform interface to log messages
to multiple output streams. These streams may be directed to the
console or text files. Lucee loggers are non-intrusive, i.e.  class
and function interfaces do not need to be modified to use the
logger. A single log call puts the message in a controlled manner into
multiple streams at the same time.

Reference for :class:`Lucee::Logger`
------------------------------------

.. class:: Lucee::Logger

  The class :class:`Lucee::Logger` provides methods to create new
  loggers and obtain streams for logging messages.

  .. staticmethod:: create(const std::string& nm)

    :param nm: name of new logger
    :rtype: Lucee::Logger&

    Create a new logger named ``nm`` and return a reference to
    it. Once a logger with a given name is created it can be accessed
    by the :meth:`get` method. A dot in the logger name creates a
    child logger::

      Lucee::Logger& cl = Lucee::Logger::create("parent.child");

    Child loggers inherit log-handlers and its log-level from the
    parent logger.
  
  .. staticmethod:: get(const std::string& nm)

    :param nm: name of logger to get
    :rtype: Lucee::Logger&

    Get a reference to an already created logger. An exception is
    thrown if the logger does not exist.

  .. method:: getName()

     :rtype: std::string

     Get name of logger.

  .. method:: debug(const std::string& msg)

    :param msg: Message to log

    Log a message with log-level "debug".

  .. method:: info(const std::string& msg)

    :param msg: Message to log

    Log a message with log-level "info".

  .. method:: warning(const std::string& msg)

    :param msg: Message to log

    Log a message with log-level "warning".

  .. method:: error(const std::string& msg)

    :param msg: Message to log

    Log a message with log-level "error".

  .. method:: critical(const std::string& msg)

    :param msg: Message to log

    Log a message with log-level "critical".

  .. method:: getDebugStream()

    :rtype: Lucee::LogStream

    Get a stream to log ``debug`` messages.

  .. method:: getInfoStream()

    :rtype: Lucee::LogStream

    Get a stream to log ``info`` messages.

  .. method:: getWarningStream()

    :rtype: Lucee::LogStream

    Get a stream to log ``warning`` messages.

  .. method:: getErrorStream()

    :rtype: Lucee::LogStream

    Get a stream to log ``error`` messages.

  .. method:: getCriticalStream()

    :rtype: Lucee::LogStream

    Get a stream to log ``critical`` messages.

  .. method:: setLevel(Lucee::LogMsgLevels level)

    :param level: Log-level specified as one of Lucee::LogMsgLevels flags.

    Set the log-level of the logger. The ``level`` parameter must be
    one of the enumerated types defined in the ``Lucee::LogMsgLevels``
    enumeration.

  .. method:: setLevel(const std::string& level)

    :param level: Log-level specified as a string.

    Set the log-level of the logger. The ``level`` parameter must be
    one of "debug", "info", "warning", "error", or "critical".

  .. method:: getLevel()

    :rtype: Lucee::LogMsgLevels

    Get the log-level of the logger.

  .. method:: getLevelStr()

    :rtype: std::string

    Get the log-level of the logger as string.

  .. method:: disable()

    Disable all logging to this logger. After this call no log
    messages will be sent to any of the log-handlers.

  .. method:: enable()

    Enable logging to this logger. This method can be called to
    re-enable logging after a call to the :meth:`disable` method.

Reference for :class:`Lucee::LogRecordHandler`
----------------------------------------------

.. class:: Lucee::LogRecordHandler

  The class :class:`Lucee::LogRecordHandler` is the base class for
  classes that handle writing log messages to output streams. This
  class can not be instantiated. Instead one of its children must be
  used. Manging the lifetime of the handler is left to the user. Once
  a handler goes out of scope it automatically detaches itself from
  all the loggers it is attached to.

  .. method:: attachToLogger(const std::string& name)

    :param name: name of logger to which handler should be attached

    Attach this handler to the logger with the specified name.

  .. method:: attachToLogger(Lucee::Logger& logger)

    :param logger: logger to which handler should be attached

    Attach this handler to the specified logger.

  .. method:: detachFromLogger(const std::string& name)

    :param name: name of logger from which to detach handler
    :rtype: bool

    Remove handler from logger with specified name. Return ``true`` if
    detaching worked, ``false`` otherwise.

  .. method:: detachFromLogger(Lucee::Logger& logger)

    :param logger: logger from which to detach handler
    :rtype: bool

    Remove handler from specified logger. Return ``true`` if detaching
    worked, ``false`` otherwise.

  .. method:: loggerNames()

    :rtype: std::vector<std::string>

    Get list of loggers to which handler is attached.

Reference for :class:`Lucee::StreamHandler`
-------------------------------------------

.. class:: Lucee::StreamHandler(std::ostream& stream)

  :param stream: Standard C++ I/O stream to attach to

  This class is derived from :class:`Lucee::LogRecordHandler` and is
  used to create a handler attached to any Standard C++ I/O stream
  object.

Reference for :class:`Lucee::FileHandler`
-----------------------------------------

.. class:: Lucee::FileHandler(const std::string& fn, std::ios_base::openmode mode)

  :param fn: Name of file to log messages
  :param mode: Mode to open file. By default this is ``std::ios_base::trunc``.

  This class is derived from :class:`Lucee::LogRecordHandler` and is
  used to log messages to the file named ``fn``. By default this file
  is truncated (i.e. its contents discarded) if it already exists and
  is created if it does not exist.
   

Example usage
-------------

The logging interface is defined in the ``lclogger.h`` header file. In
the following example a logger with id ``myLogger`` is created and its
log-level is set to ``info``. This means that all messages which have
log-level equal or higher than ``info`` will be logged, and all other
messages will be ignored. Hence, for this particular logger, all
``debug`` messages will be ignored. The levels, in increasing order
are, ``debug``, ``info``, ``warning``, ``error`` and ``critical``::

  Lucee::Logger& logger = Lucee::Logger::create("myLogger");
  logger.setLevel("info");

Once the logger is created we attach it to multiple handlers::

  Lucee::StreamHandler strmHndlr(std::cout);
  strmHndlr.attachToLogger(logger);

  Lucee::FileHandler fileHndlr("myLogFile");
  fileHndlr.attachToLogger(logger);

With these handlers all log messages will go to the console and a file
named "myLogFile". This completes the logger setup, which needs to be
done, in general, only once at the top level of the application.

Once loggers are created and handlers attached, they can be accessed
from any point in the code using the :meth:`get` method. To log
messages, log-streams are used::

  Lucee::Logger& logger = Lucee::Logger::get("myLogger");

  Lucee::LogStream dbgStrm = logger.getDebugStream();
  dbgStrm << "This is a debug message" << std::endl;

  Lucee::LogStream infoStrm = logger.getInfoStream();
  infoStrm << "This is a informational message" << std::endl;

As the logger's verbosity is set to ``info``, the first message will
not appear in the console or the file, but the second message will.

A hierarchical system of loggers can be created by inheriting from an
existing logger. Child loggers are created by using a dot in the
logger name::

  Lucee::Logger& childLogger = Lucee::Logger("myLogger.child");

This will create a child logger with id ``child``, which inherits all
its handlers and its log-level from its parent. However, additional
handlers can be added and log-level set independently::

  childLogger.setLevel("debug");
  Lucee::FileHandler childFileHndlr("childLogFile");
  childFileHndlr.attachToLogger(childLogger);

When messages are logged to the ``childLogger`` they will appear in
the parent's handlers as well as it own::

  Lucee::LogStream dbgStrm = childLogger.getDebugStrm();
  dbgStrm << "This is a debug message" << std::endl;

Due to the log-level of the ``childLogger`` the debug message will
appear in the ``childLogFile`` but not in the parent's handler or the
console.

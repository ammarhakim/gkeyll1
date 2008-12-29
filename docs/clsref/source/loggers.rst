.. -*- rst -*-

:mod:`Lucee::Logger` --- Logging messages to multiple streams
=============================================================

.. module:: Lucee::Logger
 :synopsis: Classes for logging messages to multiple streams.

The Lucee logging system provides an uniform interface to log messages
to multiple output streams. These streams may be directed to the
console, text files or even GUI elements. A good logging system must
be non-intrusive. This means that class and function interfaces must
not be modified to use the logger. In addition, a single log call
should put the message in a controlled manner into multiple streams at
the same time. For example, only error messages may need to be
displayed in the console or dialog boxes, while debugging messages may
be directed to a file to be examined later. Further, the application
top level should be able to "choke" the logger, i.e., with a single
call disable all log messages or reduce their verbosity. This is
particularly useful during the development and testing phase of the
application when log messages are most useful.

In the following sections the design and use of a simple logging
system is described.

Creating a logger
-----------------

Applications need different loggers for different parts of the
code. At the top level an application creates a set of loggers for use
in its various subsystems. Each of these loggers is given a name and
functions and classes wishing to log use these names to get the
logger. Once the logger is obtained, various streams can be created.

When a non-existent logger is accessed for the first time it is
created and registered by the name used. For example, to create a
logger with the name "root-logger" the application would use::

  Lucee::Logger *logger = Lucee::Logger::get("root-logger");

This creates a log named "root-logger". Once a logger is created, its
verbosity can be set. The highest verbosity level is "debug", in which
case all messages are logged. The other levels, in decreasing levels
of verbosity are "info", "warning", "error" and "critical".::

  logger->setLevel("debug");

Various streams can now be attached to the logger. To add a log file
named "log.txt"::

  Lucee::LogRecordHandler *h = new Lucee::FileHandler("log.txt");
  logger->addHandler(h);

Additional streams can also be added to the same logger. For example,
to log to the console::

  Lucee::LogRecordHandler *h = new Lucee::StreamHandler();
  logger->addHandler(h);

After these handlers are added all debug log messages are directed to
the file "log.txt" and the console at the same time.

Loggers support inheritance. To create a logger which not only logs to
the streams provided by the "root-logger" and also its own specific
streams, one can do::

  Lucee::Logger *myLogger = Lucee::Logger::get("root-logger.child");

This will create a new logger *myLogger* which inherits all streams and
settings from the "root-logger". The child logger can now add its own
streams::

  Lucee::LogRecordHandler *h = new Lucee::FileHandler("childLog.txt");
  myLogger->addHandler(h);

When the the "root-logger.child" is used the log messages go to the
"root-logger" stream as well as the child logger stream. Note that the
child loggers can set their own verbosity level independently of the
parent::

  myLogger->setLevel("error");

In this this case only error (and lower verbosity) messages will
appear in the child logger's stream, while all more verbose messages
will appear in the parent's steams.
  

Using the logging system
------------------------

To illustrate the use of the loggers, consider that a logger named
"root-logger" has been created by the top level application. To use
this logger from some function, we first obtain the logger by its
name.::

  Lucee::Logger *logger = Lucee::Logger::get("root-logger");

Once the is logger obtained various streams can be created. The log
stream is a C++ iostream object and can be used as any other C++
stream.::

  Lucee::LogStream debStrm = logger->getDebugStream();
  debStrm << "Debug message" << std::endl;

This will put the log message, depending on global settings of the
logger, into the various files and output elements attached to the
debug stream. For example, the message will appear on the console and
the "log.txt" file on the disk.

Other log streams can also be obtained::

  Lucee::LogStream wrnStrm = logger->getWarningStream();
  wrnStrm << "Warning message" << std::endl;

The warning messages will be logged only if the logger's level is set
to "warning" or greater verbosity.

The level of the logger can be changed::

  logger->setLevel("warning");

Now the logger will choke all "debug" and "info" messages::

  dbgStrm << "This message is not logged" << std::endl;

All logging can be disabled::

  logger->disable();

As loggers are global objects, a change made in one function affects
all other functions.
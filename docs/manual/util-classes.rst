Utility Classes and Functions
-----------------------------

Lucee provides several utility classes to perform miscellenous
jobs. These are described below.

``Lucee::CmdLineArgs``: Parsing Command Lines
+++++++++++++++++++++++++++++++++++++++++++++

The class ``Lucee::CmdLineArgs`` provides a simple class to parse
command line arguments. Two types of arguments are supported,
*variables* and *switches*. Variable arguments have an associated
value, while switches simply represent a flag. For example, to parse a
command line like::

  program -i file.txt -xml -o ouput

the following code fragment can be used::

  void main(int argc, char *agrv[]) {
    Lucee::CmdLineArgs cmd;
    // set up parser with expected values and help string
    cmd.addSwitch("xml", "Set to output XML");
    cmd.addArg("i", "INPUT", "Input file");
    cmd.addArg("o", "OUTPUT", "Output prefix");

    // parse command line
    cmd.parse(argc, argv);

    // check if "xml" switch was specified
    if (cmd.hasSwitch("xml")) {
      // XML specific code
    }
 
    // check input file name
    std::string inputFile("default");
    if (cmd.hasArg("i")) 
      inputFile = cmd.getArg("i");
  }

The class also provides a method for parsing out *extra arguments*,
i.e. those which are not switches or arguments. For example, in the
command line::
 
   program -xml file-1.c file-2.c file-3.c

the names ``file-1.c``, ``file-2.c`` and ``file-3.c`` are extra
arguments and can be retrieved using::

  std::vector<std::string> extra = cmd.getExtraArgs();

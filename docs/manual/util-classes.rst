Utility classes and functions
-----------------------------

Lucee provides several utility classes to perform miscellenous
jobs. These are described below.

``Lucee::Region``: Half-open regions in N-dimensional space
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The class ``Lucee::Region`` represents half-open regions in
N-dimensional space. For example, consider a region in the space of
two-dimensional reals defined by :math:`[a_0,b_0) \times
[a_1,b_1)`. This is represented as::

  double lower[2] = {a0, a1};
  double upper[2] = {b0, b1};
  Lucee::Region<2, double> rgn(lower, upper);

Region objects allow the computation of volumes::

  double vol = rgn.getVolume();

and is defined as, for the example above, :math:`(b0-a0)(b1-a1)`. Also
provided is a method to check if a given point is inside a region::

  double p1[2] = {0.5, 0.5};
  bool isIn = rng.isInside(p1);

Regions in multi-dimensional integer spaces are particularly useful
for constructing index spaces::

  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> irgn(lower, upper);

This will create a region spanning :math:`[3,13) \times [4,14)`,
i.e.::

  int pa[2] = {3, 10};
  bool paIn = irgn.isInside(pa); // true;

  int pb[2] = {13, 10};
  bool pbIn = irgn.isInside(pa); // false;

Regions can be extended. Consider the region :math:`[1,10) \times
[2,10)`. This can be expanded by specifying the amount to expand
by. For example::

  int lower[2] = {1,2};
  int upper[2] = {10, 10};
  Lucee::Region<2, int> rgn(lower, upper);

  int lowerExp[2] = {1, 2};
  int upperExp[2] = {1, 2};
  Lucee::Region<2, int> expRgn = rgn.extend(lowerExp, upperExp);

will create a new box ``expRgn`` :math:`[0,11) \times [0,12)`.

Regions can be inflated, i.e. their dimension can be increased by one
and the lower and upper bound of the newly added dimension
specified. For example::

  int lower[2] = {1,2};
  int upper[2] = {10, 10};
  Lucee::Region<2, int> rgn(lower, upper);

  Lucee::Region<3, int> infRgn = rgn.inflate(0, 3);

will create a 3-dimensional region :math:`[1,10) \times [2,10) \times
[0,3)` from a 2-dimensional region :math:`[1,10) \times [2,10)`.

Intersection of two regions can also be computed. For example::

  int lo1[2] = {0, 0};
  int up1[2] = {10, 10};
  Lucee::Region<2, int> ibox1(lo3, up3);

  int lo2[2] = {5, 5};
  int up2[2] = {15, 15};
  Lucee::Region<2, int> ibox2(lo4, up4);

  Lucee::Region<2, int> ibox12 = ibox1.intersect(ibox2);

This will create the region :math:`[5,10)\times[5,10)`.


``Lucee::CmdLineArgs``: Parsing command line arguments
++++++++++++++++++++++++++++++++++++++++++++++++++++++

The class ``Lucee::CmdLineArgs`` provides a simple class to parse
command line arguments. Two types of arguments are supported,
*variables* and *switches*. Variable arguments have an associated
value, while switches simply represent a flag. For example, to parse a
command line like::

  program -i file.txt -xml -o ouput

the following code fragment can be used::

  void main(int argc, char *agrv[]) 
  {
    Lucee::CmdLineArgs cmd;
    // set up parser with expected values and help string
    cmd.addSwitch("xml", "Set to output XML");
    cmd.addArg("i", "INPUT", "Input file");
    cmd.addArg("o", "OUTPUT", "Output prefix");

    // parse command line
    cmd.parse(argc, argv);

    // check if "xml" switch was specified
    if (cmd.hasSwitch("xml")) 
    {
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

/**
 * $Id: gendeclstatics.cxx 58 2012-09-15 13:43:53Z jrobcary $
 */

template <class TYPE> class X {
  public:
    static int r;
};
template <class TYPE> int X<TYPE>::r = 0;

int main (int argc, char* argv[]) {
  X<double> x;
  int rr = x.r + X<float>::r;
};


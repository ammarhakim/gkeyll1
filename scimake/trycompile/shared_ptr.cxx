// $Id: shared_ptr.cxx 195 2013-02-03 18:18:39Z jrobcary $

#include <sci_shared_ptr>

struct S {
  int i;
};

int main(int argc, char** argv) {
  sci_shared_ptr<S> sptr;
  return 0;
}

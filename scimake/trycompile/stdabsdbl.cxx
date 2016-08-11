/**
 * $Id: stdabsdbl.cxx 58 2012-09-15 13:43:53Z jrobcary $
 *
 * Determine whether the compiler knows std::abs<double>.
 */

#include <cmath>

int main(int argc, char** argv) {
  double a = 0;
  double b = std::abs(a);
  return 0;
}


/**
 * @file	LcTest.h
 *
 * @brief	Macros for unit testing
 *
 */

#ifndef LC_TEST_H
#define LC_TEST_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
#  include <mpi.h>
#endif

/**
 * Code in this file provides a simple test framework to rapidly test
 * various bits of code. Most functionality is provided via macros.
 */

// std includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/**
 * Keeps track of how many tests have passed or failed
 */
class LcTestCounter
{
  public:
    LcTestCounter()
      : passed(0), failed(0) {
    }

    int passed;
    int failed;
    std::ofstream msgFile;
    std::ofstream msgFileFail;

    std::vector<std::string> failedTests;

    void addFailedTest(int line, const char *file, const std::string& s) {
      std::ostringstream lf;
      lf << s << " In file " << file;
      lf << " at line " << line;
      failedTests.push_back(std::string(lf.str()));
    }

    void showFailedTests() {
      std::vector<std::string>::iterator i;
      for (i=failedTests.begin(); i!=failedTests.end(); ++i) {
        //std::cout << "FAILED TEST: " << *i << std::endl;
        msgFileFail   << "FAILED TEST: " << *i << std::endl;
      }
    }
    
    void clearFailedTests() {
      failedTests.erase( failedTests.begin(), failedTests.end() );
    }
};
// `tc` is global to ease counting process
static LcTestCounter __tc;

/**
 * The LC_BEGIN_TESTS and LC_END_TEST macros can be used to sandwich a
 * set of LC_ASSERT and LC_RAISES macros defined below. These macros
 * count the number of passed and failed tests, printing them out
 * after the tests are done
 */
#ifdef HAVE_MPI

#define _LC_MPI_BEGIN_TESTS(file)                                        \
    MPI_Init(&argc, &argv);                                             \
    do {                                                                \
      __tc.passed = 0;                                                  \
      __tc.failed = 0;                                                  \
      int myrank;                                                       \
      std::ostringstream fname, fnameF;                                 \
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);                           \
      fname << file << "-" << myrank << ".log";                         \
      __tc.msgFile.open(fname.str().c_str(), std::fstream::out);        \
      fnameF << file << "-failed-" << myrank << ".log";                     \
      __tc.msgFileFail.open(fnameF.str().c_str(), std::fstream::out);        \
    } while (0)

#define _LC_MPI_END_TESTS                                                \
    int __tflag;                                                        \
    do {                                                                \
      MPI_Barrier(MPI_COMM_WORLD);                                      \
      unsigned totalPassed, totalFailed;                                \
      MPI_Reduce(&__tc.passed, &totalPassed, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); \
      MPI_Reduce(&__tc.failed, &totalFailed, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); \
      int __myrank;                                                     \
      MPI_Comm_rank(MPI_COMM_WORLD, &__myrank);                         \
      if (__myrank == 0) {                                              \
        std::cout << "PASSED = " << totalPassed <<  ". FAILED = " << totalFailed << std::endl; \
      }                                                                 \
      __tc.showFailedTests();                                           \
      __tflag = totalFailed > 0 ? 1 : 0;                                \
    } while (0);                                                        \
    MPI_Finalize();                                                     \
    return __tflag;

#else
# define _LC_MPI_BEGIN_TESTS(file)
# define _LC_MPI_END_TESTS
#endif // HAVE_MPI

#define _LC_BEGIN_TESTS(file)                                            \
    do {                                                                \
      __tc.passed = 0;                                                  \
      __tc.failed = 0;                                                  \
      std::ostringstream fname, fnameF;                                 \
      fname << file << "-" << "0" << ".log";                            \
      __tc.msgFile.open(fname.str().c_str(), std::fstream::out);        \
      fnameF << file << "-failed-" << "0" << ".log";                     \
      __tc.msgFileFail.open(fnameF.str().c_str(), std::fstream::out);     \
    } while (0)

#define _LC_END_TESTS                                                    \
    int __tflag;                                                        \
    do {                                                                \
      std::cout << "PASSED = " << __tc.passed <<  ". FAILED = " << __tc.failed << std::endl; \
      __tc.showFailedTests();                                           \
      __tflag = __tc.failed > 0 ? 1 : 0;                                \
      __tc.passed = 0;                                                  \
      __tc.failed = 0;                                                  \
      __tc.clearFailedTests();                                          \
    } while (0);                                                        \
    return __tflag;

#ifdef HAVE_MPI
# define LC_BEGIN_TESTS _LC_MPI_BEGIN_TESTS
# define LC_END_TESTS _LC_MPI_END_TESTS  
#else
# define LC_BEGIN_TESTS _LC_BEGIN_TESTS
# define LC_END_TESTS _LC_END_TESTS  
#endif

/** 
 * The following macro can be used for running test cases. To test
 * some statement use, for example,
 *
 * LC_ASSERT("Running test 3", 2==2);
 *
 * The first macro parameter is the message while the second parameter
 * is the actual statement to test. This statement should be a
 * predicate, i.e. should evaluate to either true or false. If a test
 * fails the program does not halt.
 *
 */
#define LC_ASSERT(msg,expr)                                             \
    do {                                                                \
      __tc.msgFile << msg << "\n";                                      \
      __tc.msgFile << "   Testing if " << #expr << "\n";                \
      if ((expr) == true)                                               \
      {                                                                 \
        __tc.msgFile << "     PASSED\n";                                \
        __tc.passed++;                                                  \
      }                                                                 \
      else                                                              \
      {                                                                 \
        __tc.msgFile << "     FAILED\n";                                \
        __tc.failed++;                                                  \
        __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
      }                                                                 \
    } while (0)

/**
 * The following macro can be used for running test cases. To test
 * some statement use, for example,
 *
 *  LC_RAISES("Running test 3", 2==2, const char *);
 *
 * The first macro parameter is the message while the second parameter
 * is the actual statement to test.  This statement should be a
 * predicate, i.e. should evaluate to either true or false. If a test
 * fails the program does not halt. The third parameter is a type
 * indicating the type of exception which needs to be caught.
 */
#define LC_RAISES(msg,expr,expc)                                        \
    do {                                                                \
      try                                                               \
      {                                                                 \
        __tc.msgFile << msg << std::endl;                               \
        __tc.msgFile << "   Executing " << #expr << "...\n";            \
        __tc.msgFile << "   Expecting exception of type " << #expc << "...\n"; \
        expr;                                                           \
        __tc.msgFile << "     No exception thrown...\n";                \
        __tc.msgFile << "     FAILED\n";                                \
        __tc.failed++;                                                  \
        __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
      }                                                                 \
      catch(expc e)                                                     \
      {                                                                 \
        __tc.msgFile << "     Caught exception of type " << #expc << "...\n"; \
        __tc.msgFile << "     PASSED\n";                                \
        __tc.passed++;                                                  \
      }                                                                 \
      catch(...)                                                        \
      {                                                                 \
        __tc.msgFile << "     Caught exception of unexpected type\n";   \
        __tc.msgFile << "     FAILED\n";                                \
        __tc.failed++;                                                  \
        __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
      }                                                                 \
    } while (0)

/**
 * Compares two arrays and returns true if they are equal.
 *
 * @param a first array to compare
 * @param b second array to compare
 * @param length size of arrays
 * @param true if arrays are the same
 */
template<typename T>
bool arraycmp(T a[], T b[], unsigned length)
{
  for (unsigned i = 0; i < length; ++i) 
  {
    if (a[i] != b[i]) 
      return false;
  }
  return true;
}

/**
 * Compares two arrays and returns true if they are equal.
 *
 * @param a first array to compare
 * @param b second array to compare
 * @return true if vectors are the same
 */
template<typename T>
bool arraycmp(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size())
    return false;
  for (unsigned i=0; i<a.size(); ++i)
  {
    if (a[i] != b[i]) 
      return false;
  }
  return true;
}

/**
 * Compares two floating point numbers to floating-point precision.
 *
 * @param a first value to compare.
 * @param b second value to compare.
 * @param fact Factor of floating-point epsilon.
 */
template <typename T>
bool epsCmp(const T& a, const T& b, double fact=5)
{
  if (fabs(a) <= fact*std::numeric_limits<T>::epsilon())
    return fabs(a-b) <= fact*std::numeric_limits<T>::epsilon();
  return fabs(1-b/a) <= fact*std::numeric_limits<T>::epsilon();
}

/**
 * Compares two floating point numbers to floating-point precision.
 *
 * @param a first value to compare.
 * @param b second value to compare.
 * @param eps Precision to compare
 */
template <typename T>
bool diffCmp(T a, T b, double eps)
{
  return fabs(a-b) <= eps;
}

#endif // LC_TEST_H

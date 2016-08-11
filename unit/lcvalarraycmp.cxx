#include <iostream>
#include <valarray>
#include <vector>
#include <ctime>
#include <iomanip>

#include <LcVector.h>

using namespace std;

static unsigned numruns=5;

void rndm(float& a);
float getRndm();
void zerocheck(float& a);
double duration(vector<clock_t>::iterator start, vector<clock_t>::iterator stop);

int test_1()
{
  float addVal = 5.3235;
  int size = 4000000;
  float mulVal = 7.3452;
  valarray<float> valfloat(size);

  clock_t start = clock();
  for (unsigned n=0; n<numruns; ++n)
    for (unsigned i=0;i<size;i++)
      rndm(valfloat[i]);
  clock_t end = clock();
  cout << right << fixed << difftime(end,start)  << " (valarray init.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    valfloat += addVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (valarray add.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    valfloat *= mulVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (valarray mul.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    valfloat /= mulVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (valarray div.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    for (unsigned i=1; i<size-1; i++)
      valfloat[i] = (valfloat[i+1]-2*valfloat[i]+valfloat[i-1])/2.25;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (valarray diff.)" << endl;
}

int test_2()
{
  float addVal = 5.3235;
  float mulVal = 7.3452;
  int size = 4000000;
  Lucee::Vector<float> lvecfloat(size);

  clock_t start = clock();
  for (unsigned n=0; n<numruns; ++n)
    for (unsigned i=0;i<size;i++)
      lvecfloat[i]= getRndm();
  clock_t end = clock();
  cout << right << fixed << difftime(end,start)  << " (Lucee::Vector init.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    lvecfloat += addVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (Lucee::Vector add.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    lvecfloat *= mulVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (Lucee::Vector mul.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    lvecfloat /= mulVal;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (Lucee::Vector div.)" << endl;

  start = clock();
  for (unsigned n=0; n<numruns; ++n)
    for (unsigned i=1; i<size-1; i++)
      lvecfloat[i] = (lvecfloat[i+1]-2*lvecfloat[i]+lvecfloat[i-1])/2.25;
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (Lucee::Vector diff.)" << endl;
}

void rndm(float& a)
{
  a=rand()%1000000;
  a/=(1000);
}

float getRndm()
{
  float a=rand()%1000000;
  a/=(1000);
  return a;
}

void zerocheck(float& a)
{
  if (a<0) { a*=-1; } 
}

double duration(
  vector<clock_t>::iterator start, 
  vector<clock_t>::iterator stop
                )
{
  return difftime(*stop,*start) ;
}

int
main (void) 
{
  std::cout << "Valarray times ..." << std::endl;
  test_1();
  std::cout << std::endl << "Lucee::Vector times ..." << std::endl;
  test_2();

  return 0;
}

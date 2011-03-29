#include <iostream>
#include <valarray>
#include <vector>
#include <ctime>
#include <iomanip>

#include <LcVector.h>

using namespace std;

void rndm(float& a);
float getRndm();
void zerocheck(float& a);
double duration(vector<clock_t>::iterator start, vector<clock_t>::iterator stop);

int main()
{
  int i(0);
  int p(0);
  float k(5.3235),l(9.212);
  int size(4000000);

  vector<clock_t> vstart;
  vector<clock_t> vstop;
  vector<clock_t>::iterator istartclock;
  vector<clock_t>::iterator istopclock;
  vector<float> vfloat(size);
  vector<float>::iterator ifloat;
  valarray<float> valfloat(size);
  Lucee::Vector<float> lvecfloat(size);
  float *vecraw = new float[size];

  cout << "Filling float array with " << size << " random elements" << endl;
  cout << "(each element needs " << sizeof(float) << " bytes of memory, so " \
       << "we need to make use of " << size*sizeof(float)/1000000<< "MBytes.)"<< endl;

  // ======================= ORDINARY OPERATIONS =======================
  // fill array

  vstart.push_back(clock());
  //  vfloat.resize(size);
  ifloat=vfloat.begin();
  for_each(vfloat.begin(),vfloat.end(),rndm);
  vstop.push_back(clock());

  cout << "****** Phase 1: Math. operations on ordinary array" << endl;  

  // Add k to each element
  cout << "*** Adding " << k << " to each element " ;
  cout << "(10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)+=k; }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")" << endl;


  // Multiply each element with k
  cout << "*** Multiplying each element with " << k ;
  cout << " (10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)*=k; }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")"<< endl;

  // Divide each element by k
  cout << "*** Dividing each element by " << l ;
  cout << " (10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)/=l; }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")" << endl;

  // Cos of each
  cout << "*** Calculate cos for each element "  ;
  cout << " (10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)= cosf(*(ifloat+i)); }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")" << endl;

  for_each(vfloat.begin(),vfloat.end(),zerocheck);

  // Expand each element by itself
  cout << "*** Calculate element ^ element "  ;
  cout << " (10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)= pow(*(ifloat+i),*(ifloat+i)); }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")" << endl;

  // log
  cout << "*** Calculate log of each element "  ;
  cout << " (10th element was: " << *(ifloat+10) ;
  vstart.push_back(clock());
  for (i=0; i<vfloat.size(); i++) { *(ifloat+i)= log(*(ifloat+i)); }
  vstop.push_back(clock());
  cout << " now is: " << *(ifloat+10) << ")" << endl;


  // Summary
  istartclock=vstart.begin(); istopclock=vstop.begin();
  cout << "****** Clock ticks for operations: " << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (random fill)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (add.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (mult.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (div.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (cos)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (pow)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock,istopclock) << " (log)" << endl;
  cout << right << fixed << setw(13) << *(vstop.rbegin())-*(vstart.begin()) << " (total) " << endl;


  vstart.clear(); vstop.clear();
  // ======================= VALARRAY OPERATIONS =======================
  // fill valarray

  cout << "****** Phase 2: Math. operations on valarray" << endl;  
  vstart.push_back(clock());
  for (i=0;i<size;i++)
  {
    rndm(valfloat[i]);
  }
  vstop.push_back(clock());


  // Add k to each element
  cout << "*** Adding " << k << " to each element " ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat+=k;
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Multiply
  cout << "*** Multiplying each element with " << k ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat*=k;
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Divide
  cout << "*** Dividing each element by " << l ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat/=l;
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Cosinus
  cout << "*** Calculate cos for each element "  ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat=cos(valfloat);
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;


  for (i=0;i<size;i++)
  {
    zerocheck(valfloat[i]);
  }

  // Expand each element by itself
  cout << "*** Calculate element ^ element "  ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat=pow(valfloat,valfloat);
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // log
  cout << "*** Calculate log of each element "  ;
  cout << "(10th element was: " << valfloat[10] ;
  vstart.push_back(clock());
  valfloat=log(valfloat);
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Summary

  cout << "****** Clock ticks for operations: " << endl;
  istartclock=vstart.begin(); istopclock=vstop.begin();
  
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (random fill)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (add.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (mult.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (div.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (cos)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++)  << " (pow)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock,istopclock)      << " (log)" << endl;
  cout << right << fixed << setw(13) << vstop[vstop.size()-1]-vstart[0] << " (total) " << endl;

  // ======================= Lucee::Vector OPERATIONS =======================
  // fill array

  vstart.clear(); vstop.clear();

  vstart.push_back(clock());
  //  vfloat.resize(size);
  ifloat=vfloat.begin();
  cout << "****** Phase 3: Math. operations on Lucee::Vector" << endl;
  vstart.push_back(clock());
  for (i=0;i<size;i++)
  {
    lvecfloat[i]= getRndm();
  }
  vstop.push_back(clock());

  cout << "****** Phase 1: Math. operations on ordinary array" << endl;  

  // Add k to each element
  cout << "*** Adding " << k << " to each element " ;
  cout << "(10th element was: " << lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) { lvecfloat[i]+=k; }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")" << endl;


  // Multiply each element with k
  cout << "*** Multiplying each element with " << k ;
  cout << " (10th element was: " <<  lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) { lvecfloat[i]*=k; }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")"<< endl;

  // Divide each element by k
  cout << "*** Dividing each element by " << l ;
  cout << " (10th element was: " << lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) { lvecfloat[i]/=l; }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")" << endl;

  // Cos of each
  cout << "*** Calculate cos for each element "  ;
  cout << " (10th element was: " << lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) { lvecfloat[i]= cosf(lvecfloat[i]); }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")" << endl;

//   for_each(vfloat.begin(),vfloat.end(),zerocheck);
  for (i=0; i<lvecfloat.getLength(); i++)
    if (lvecfloat[i] < 0.0)
      lvecfloat[i] *= -1;

  // Expand each element by itself
  cout << "*** Calculate element ^ element "  ;
  cout << " (10th element was: " << lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) 
  { 
    lvecfloat[i]= pow(lvecfloat[i], lvecfloat[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")" << endl;

  // log
  cout << "*** Calculate log of each element "  ;
  cout << " (10th element was: " << lvecfloat[10];
  vstart.push_back(clock());
  for (i=0; i<lvecfloat.getLength(); i++) 
  { 
    lvecfloat[i]= log(lvecfloat[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << lvecfloat[10] << ")" << endl;

  // Summary
  istartclock=vstart.begin(); istopclock=vstop.begin();
  cout << "****** Clock ticks for operations: " << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (random fill)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (add.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (mult.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (div.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (cos)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (pow)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock,istopclock) << " (log)" << endl;
  //cout << right << fixed << setw(13) << *(vstop.rbegin())-*(vstart.begin()) << " (total) " << endl;
  cout << right << fixed << setw(13) << vstop[vstop.size()-1]-vstart[0] << " (total) " << endl;

  // ======================= raw OPERATIONS =======================
  // fill array

  vstart.clear(); vstop.clear();

  vstart.push_back(clock());
  //  vfloat.resize(size);
  ifloat=vfloat.begin();
  cout << "****** Phase 3: Math. operations on raw array" << endl;
  vstart.push_back(clock());
  for (i=0;i<size;i++)
  {
    vecraw[i]= getRndm();
  }
  vstop.push_back(clock());

  cout << "****** Phase 1: Math. operations on ordinary array" << endl;  

  // Add k to each element
  cout << "*** Adding " << k << " to each element " ;
  cout << "(10th element was: " << vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { vecraw[i]+=k; }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")" << endl;


  // Multiply each element with k
  cout << "*** Multiplying each element with " << k ;
  cout << " (10th element was: " <<  vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { vecraw[i]*=k; }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")"<< endl;

  // Divide each element by k
  cout << "*** Dividing each element by " << l ;
  cout << " (10th element was: " << vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { vecraw[i]/=l; }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")" << endl;

  // Cos of each
  cout << "*** Calculate cos for each element "  ;
  cout << " (10th element was: " << vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { vecraw[i]= cosf(vecraw[i]); }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")" << endl;

//   for_each(vfloat.begin(),vfloat.end(),zerocheck);
  for (i=0; i<size; i++)
    if (vecraw[i] < 0.0)
      vecraw[i] *= -1;

  // Expand each element by itself
  cout << "*** Calculate element ^ element "  ;
  cout << " (10th element was: " << vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) 
  { 
    vecraw[i]= pow(vecraw[i], vecraw[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")" << endl;

  // log
  cout << "*** Calculate log of each element "  ;
  cout << " (10th element was: " << vecraw[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) 
  { 
    vecraw[i]= log(vecraw[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << vecraw[10] << ")" << endl;

  // Summary
  istartclock=vstart.begin(); istopclock=vstop.begin();
  cout << "****** Clock ticks for operations: " << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (random fill)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (add.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (mult.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (div.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (cos)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (pow)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock,istopclock) << " (log)" << endl;
  //cout << right << fixed << setw(13) << *(vstop.rbegin())-*(vstart.begin()) << " (total) " << endl;
  cout << right << fixed << setw(13) << vstop[vstop.size()-1]-vstart[0] << " (total) " << endl;

  // ======================= more valarrary OPERATIONS =======================
  // fill array

  vstart.clear(); vstop.clear();

  vstart.push_back(clock());
  //  vfloat.resize(size);
  ifloat=vfloat.begin();
  cout << "****** Phase 3: DIRECT Math. operations on Valarray" << endl;
  vstart.push_back(clock());
  for (i=0;i<size;i++)
  {
    valfloat[i]= getRndm();
  }
  vstop.push_back(clock());

  cout << "****** Phase 1: Math. operations on ordinary array" << endl;  

  // Add k to each element
  cout << "*** Adding " << k << " to each element " ;
  cout << "(10th element was: " << valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { valfloat[i]+=k; }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;


  // Multiply each element with k
  cout << "*** Multiplying each element with " << k ;
  cout << " (10th element was: " <<  valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { valfloat[i]*=k; }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")"<< endl;

  // Divide each element by k
  cout << "*** Dividing each element by " << l ;
  cout << " (10th element was: " << valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { valfloat[i]/=l; }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Cos of each
  cout << "*** Calculate cos for each element "  ;
  cout << " (10th element was: " << valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) { valfloat[i]= cosf(valfloat[i]); }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

//   for_each(vfloat.begin(),vfloat.end(),zerocheck);
  for (i=0; i<size; i++)
    if (valfloat[i] < 0.0)
      valfloat[i] *= -1;

  // Expand each element by itself
  cout << "*** Calculate element ^ element "  ;
  cout << " (10th element was: " << valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) 
  { 
    valfloat[i]= pow(valfloat[i], valfloat[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // log
  cout << "*** Calculate log of each element "  ;
  cout << " (10th element was: " << valfloat[10];
  vstart.push_back(clock());
  for (i=0; i<size; i++) 
  { 
    valfloat[i]= log(valfloat[i]); 
  }
  vstop.push_back(clock());
  cout << " now is: " << valfloat[10] << ")" << endl;

  // Summary
  istartclock=vstart.begin(); istopclock=vstop.begin();
  cout << "****** Clock ticks for operations: " << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (random fill)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (add.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (mult.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (div.)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (cos)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock++,istopclock++) << " (pow)" << endl;
  cout << right << fixed << setw(13) << duration(istartclock,istopclock) << " (log)" << endl;
  //cout << right << fixed << setw(13) << *(vstop.rbegin())-*(vstart.begin()) << " (total) " << endl;
  cout << right << fixed << setw(13) << vstop[vstop.size()-1]-vstart[0] << " (total) " << endl;
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


/**
 * @file	lcpartests.cxx
 *
 * @brief	Testing various parallel ideas/things
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// txbase includes
#ifdef HAVE_MPI
# include <TxMpiBase.h>
#else
# include <TxSelfBase.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcDecompRegion.h>
#include <LcGlobals.h>
#include <LcRectCartGrid.h>
#include <LcTest.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <iostream>
#include <vector>

bool isValid(TxCommBase *c)
{ return (bool) c; }

class CommHolder
{
  public:
    CommHolder(TxCommBase *gc, std::vector<TxCommBase*>& lc)
      : globalComm(gc)
    {
      commList.resize(lc.size());
      for (unsigned i=0; i<lc.size(); ++i)
        commList[i] = lc[i];
      for (unsigned i=0; i<lc.size(); ++i)
        if (isValid(lc[i]))
          localComm = lc[i];
    }

    TxCommBase* getComm() const
    { return localComm; }

    unsigned getGlobalRank() const
    { return globalComm->getRank(); }

    unsigned getLocalRank() const
    { return localComm->getRank(); }

  private:
    TxCommBase *globalComm;
    std::vector<TxCommBase*> commList;
    TxCommBase *localComm;
};

void
test_0(TxCommBase *comm)
{
// split communicator
  std::vector<int> subRanks(2);
  std::vector<TxCommBase*> subComm(2);

  subRanks[0] = 0; subRanks[1] = 1;
  subComm[0] = comm->createSubComm(subRanks);

  subRanks[0] = 2; subRanks[1] = 3;
  subComm[1] = comm->createSubComm(subRanks);

  if (comm->getRank() == 0)
    LC_ASSERT("Testing if communicator 0 is valid", isValid(subComm[0]));
  if (comm->getRank() == 1)
    LC_ASSERT("Testing if communicator 0 is valid", isValid(subComm[0]));

  if (comm->getRank() == 2)
    LC_ASSERT("Testing if communicator 0 is valid", !isValid(subComm[0]));
  if (comm->getRank() == 3)
    LC_ASSERT("Testing if communicator 0 is valid", !isValid(subComm[0]));

  if (comm->getRank() == 0)
    LC_ASSERT("Testing if communicator 0 is valid", !isValid(subComm[1]));
  if (comm->getRank() == 1)
    LC_ASSERT("Testing if communicator 0 is valid", !isValid(subComm[1]));

  if (comm->getRank() == 2)
    LC_ASSERT("Testing if communicator 0 is valid", isValid(subComm[1]));
  if (comm->getRank() == 3)
    LC_ASSERT("Testing if communicator 0 is valid", isValid(subComm[1]));
}

void
test_1(TxCommBase *comm)
{
// split communicator
  std::vector<int> subRanks(2);
  std::vector<TxCommBase*> subComm(2);

  subRanks[0] = 0; subRanks[1] = 1;
  subComm[0] = comm->createSubComm(subRanks);

  subRanks[0] = 2; subRanks[1] = 3;
  subComm[1] = comm->createSubComm(subRanks);

  CommHolder commHolder(comm, subComm);
  LC_ASSERT("Testing if localComm is always valid", isValid(commHolder.getComm()));
}

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcpartests");

// get hold of global communicator
  TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
    ::Instance().comm;

  test_0(comm);
  test_1(comm);

  LC_END_TESTS;
  return 0;
}

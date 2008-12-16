/**
 * @file	lcobjcreator.h
 *
 * @brief	Class for handling object creation.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_OBJ_CREATOR
#define LC_OBJ_CREATOR

// lib includes
#include <lcexcept.h>

// std includes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
// forward declare register class
  template <typename B> class ObjRegisterBase;

/**
 * Base class for object creators. All action is in the derive class
 * Lucee::ObjCreator.
 */
  template <class B>
  class ObjCreatorBase
  {
    public:
/** Type definition for map of names to register objects */
      typedef std::map<std::string, Lucee::ObjRegisterBase<B>*, std::less<std::string> > CreatorMap_t;
/** Type definition for pair of names to register objects */
      typedef std::pair<std::string, Lucee::ObjRegisterBase<B>*> CreatorPair_t;

/**
 * Add a new creator object into the list of available creators.
 *
 * @param nm Name of the creator
 * @param b Pointer to the creator base object
 */
      static void addCreator(const std::string& nm, Lucee::ObjRegisterBase<B>* b) 
      {
        // if doing this first time, make a new map
        if (!creators) creators = new CreatorMap_t();
        // insert object
        creators->insert(CreatorPair_t(nm, b));
      }

/**
 * Remove a creator from the list
 *
 * @param nm Name of creator to remove
 */
      static void removeCreator(const std::string& nm) 
      {
        creators->erase(nm);
      }

    protected:
/** Map of names to creator objects */
      static CreatorMap_t *creators;
  };

/**
 * Class to create objects. This class is used in conjunction with the
 * Lucee::ObjRegister class. Once an object has been registered this
 * class can be used to create derived classes using their names.
 */
  template <class B>
  class ObjCreator : public ObjCreatorBase<B>
  {
    public:
      typedef std::map<std::string, Lucee::ObjRegisterBase<B>*, std::less<std::string> > CreatorMap_t;

/**
 * Get a new object whose creator has the given name. The returned
 * object points to the base class.
 *
 * @param nm Name of the creator.
 * @return pointer to newly created object.
 */
      static B* getNew(const std::string& nm) 
      {
        if (ObjCreatorBase<B>::creators)
        {
          typename CreatorMap_t::iterator i = ObjCreatorBase<B>::creators->find(nm);
          if (i != ObjCreatorBase<B>::creators->end())
            return i->second->getNew();
        }
        Lucee::Except ex;
        ex << "Creator for class " << nm << " not found";
        throw ex;
      }

/**
 * Get a list of registered names.
 *
 * @return List of registered names.
 */
      static std::vector<std::string> registeredNames() 
      {
        std::vector<std::string> names;
        if (ObjCreatorBase<B>::creators)
        { 
          typename CreatorMap_t::const_iterator i;
          for (i=ObjCreatorBase<B>::creators->begin(); i!=ObjCreatorBase<B>::creators->end(); ++i)
            names.push_back( (*i).first );
        }
        return names;
      }

/**
 * Check if creator with given name is registered.
 *
 * @param nm Name of creator to check.
 * @return true if creator exisits, false otherwise.
 */
      static bool has(const std::string& nm) 
      {
        if (ObjCreatorBase<B>::creators)
        {
          typename CreatorMap_t::iterator i = ObjCreatorBase<B>::creators->find(nm);
          if (i != ObjCreatorBase<B>::creators->end())
            return true;
        }
        return false;
      }
  };

  // initialize map
  template <class B>
  std::map<std::string, Lucee::ObjRegisterBase<B>*, std::less<std::string> >*
  Lucee::ObjCreatorBase<B>::creators = 0;
}

#endif // LC_OBJ_CREATOR

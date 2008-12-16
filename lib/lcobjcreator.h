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
      static CreatorMap_t *creators;
  };

  template <class B>
  class ObjCreator : public ObjCreatorBase<B>
  {
    public:
      typedef std::map<std::string, Lucee::ObjRegisterBase<B>*, std::less<std::string> > CreatorMap_t;

/**
 * Get a new object whose creator has the given name. The returned
 * object points to the base class.
 *
 *   @param nm Name of the creator.
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

/**
 * 
 */
  template<class B>
  class ObjRegisterBase
  {
    public:
/**
 * Register a new object which will be created by its given name.
 *
 * @param nm Name by which object will be created.
 */
      ObjRegisterBase(const std::string& nm) 
        : name(nm) 
      {
        Lucee::ObjCreator<B>::addCreator(name, this);
      }

/**
 * Delete the object, unregistering the creator class.
 */
      virtual ~ObjRegisterBase() 
      {
        Lucee::ObjCreator<B>::removeCreator(name);
      }

/**
 * Return a newly allocated object. Must be provided by derived
 * classes.
 *
 * @return Newly allocated object.
 */
      virtual B* getNew() = 0;

    private:
/** Name by which this object is to be created */
      std::string name;
  };

  // initialize map
  template <class B>
  std::map<std::string, Lucee::ObjRegisterBase<B>*, std::less<std::string> >*
  Lucee::ObjCreatorBase<B>::creators = 0;
}

#endif // LC_OBJ_CREATOR

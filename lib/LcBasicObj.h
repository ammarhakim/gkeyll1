/**
 * @file	LcBasicObj.h
 *
 * @brief	Interface class for basic Lucee objects.
 */

#ifndef LC_BASIC_OBJ_H
#define LC_BASIC_OBJ_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaFuncMap.h>
#include <LcLuaTable.h>

// txbase includes
#include <TxCommBase.h>
#include <TxIoBase.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Represents a simple object in Lucee that can be initialized from a
 * Lua script. Its main purpose is to provide a method of holding
 * pointers to derived type objects in maps. It also provides a
 * default implementation of a static function to add Lua callable
 * methods to Lucee.
 */
  class BasicObj
  {
    public:
/**
 * Create an object.
 */
      BasicObj();

/**
 * Create an object with specified name.
 *
 * @param nm Name of object.
 */
      BasicObj(const std::string& nm);

/**
 * Destroy object.
 */
      virtual ~BasicObj();

/**
 * Get name of object.
 *
 * @return Name of object.
 */
      std::string getName() const;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize object.
 */
      virtual void initialize();

/**
 * Set the communicator
 *
 * @param ci valid communicator
 */
      void setComm(TxCommBase* ci);

/**
 * Set the communicator
 *
 * @param ci data communicator
 */
      void setDataComm(TxCommBase* ci);

/**
 * Get the communicator
 *
 * @return ci the communicator
 */
      TxCommBase* getComm() const;

/**
 * Get the communicator for I/O. Note that this communicator IS NOT
 * valid on all processors. Hence, it should be tested to ensure that
 * it is valid before using it.
 *
 * @return ci the communicator
 */
      TxCommBase* getDataComm() const;

/**
 * Check if this object is valid on this rank.
 *
 * @return true if object is valid.
 */
      bool isValidOnRank() const;

/**
 * Get type of object.
 *
 * @return Type of object.
 */
      std::string getType() const {
        return typeid(*this).name();
      }

/**
 * Set the type of the base class for object.
 */
      template <typename B>
      void setBaseType()
      {
        baseType = typeid(B).name();
      }

/**
 * Set the type of the derived class for object.
 */
      template <typename D>
      void setDerivedType()
      {
        derivedType = typeid(D).name();
      }

/**
 * Get base type of object.
 *
 * @param base type ID of object.
 */
      std::string getBaseType() const 
      { return baseType; }

/**
 * Get derived type of object.
 *
 * @param derived type ID of object.
 */
      std::string getDerivedType() const 
      { return derivedType; }
      
/**
 * Default method that performs registration of Lua functions. This
 * function does nothing: if derived classes need to register Lua
 * callable functions they must provide this method. Methods should be
 * added to the lfm object by calling the appendFunc() method.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    protected:
/**
 * Set name of object. This should be called by derived classes to set
 * their names.
 *
 * @param nm Name of solver.
 */
      void setName(const std::string& nm);
      
    private:
/** Name of the object */      
      std::string nm;
/** Base type of object */
      std::string baseType;
/** Derived type of object */
      std::string derivedType;

/** Communicator on which object is valid */
      TxCommBase* comm;
/** Communicator for performing I/O */
      TxCommBase* dataComm;
  };
}

#endif // LC_BASIC_OBJ_H

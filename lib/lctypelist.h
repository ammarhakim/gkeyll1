/**
 * @file	lctypelist.h
 *
 * @brief	Set of classes for handling typelists.
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_TYPE_LIST_H
#define LC_TYPE_LIST_H

namespace Lucee
{

  class NullType; // does nothing, just serves as a sentinel

/**
 * TypeList class provides a means of defining a list of types
 * (hence the name "typelist"). These are very useful for generic
 * programming, for example can be used to create abstract factories
 * using the GenScatterHier template class. The macros below make it
 * easy to generate typelist of a given length quite easily.
 */
  template<class T, class U>
  struct TypeList 
  {
      typedef T Head;
      typedef U Tail;
  };

// The following macros were generated automatically. They make
// constructing new typelists quite easy. For example, a four element
// typelist can be created using
//
// LC_TYPELIST_4(int, usigned int, char, unsigned char)
//
#define LC_TYPELIST_1(t1)                       \
  TypeList<t1, NullType>
#define LC_TYPELIST_2(t1, t2)                   \
  TypeList<t1, LC_TYPELIST_1(t2) >
#define LC_TYPELIST_3(t1, t2, t3)               \
  TypeList<t1, LC_TYPELIST_2(t2, t3) >
#define LC_TYPELIST_4(t1, t2, t3, t4)           \
  TypeList<t1, LC_TYPELIST_3(t2, t3, t4) >
#define LC_TYPELIST_5(t1, t2, t3, t4, t5)       \
  TypeList<t1, LC_TYPELIST_4(t2, t3, t4, t5) >
#define LC_TYPELIST_6(t1, t2, t3, t4, t5, t6)           \
  TypeList<t1, LC_TYPELIST_5(t2, t3, t4, t5, t6) >
#define LC_TYPELIST_7(t1, t2, t3, t4, t5, t6, t7)       \
  TypeList<t1, LC_TYPELIST_6(t2, t3, t4, t5, t6, t7) >
#define LC_TYPELIST_8(t1, t2, t3, t4, t5, t6, t7, t8)           \
  TypeList<t1, LC_TYPELIST_7(t2, t3, t4, t5, t6, t7, t8) >
#define LC_TYPELIST_9(t1, t2, t3, t4, t5, t6, t7, t8, t9)       \
  TypeList<t1, LC_TYPELIST_8(t2, t3, t4, t5, t6, t7, t8, t9) >
#define LC_TYPELIST_10(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)         \
  TypeList<t1, LC_TYPELIST_9(t2, t3, t4, t5, t6, t7, t8, t9, t10) >
#define LC_TYPELIST_11(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)    \
  TypeList<t1, LC_TYPELIST_10(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11) >
#define LC_TYPELIST_12(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12) \
  TypeList<t1, LC_TYPELIST_11(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12) >
#define LC_TYPELIST_13(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13) \
  TypeList<t1, LC_TYPELIST_12(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13) >
#define LC_TYPELIST_14(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14) \
  TypeList<t1, LC_TYPELIST_13(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14) >
#define LC_TYPELIST_15(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15) \
  TypeList<t1, LC_TYPELIST_14(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15) >
#define LC_TYPELIST_16(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16) \
  TypeList<t1, LC_TYPELIST_15(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16) >
#define LC_TYPELIST_17(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17) \
  TypeList<t1, LC_TYPELIST_16(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17) >
#define LC_TYPELIST_18(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18) \
  TypeList<t1, LC_TYPELIST_17(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18) >
#define LC_TYPELIST_19(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19) \
  TypeList<t1, LC_TYPELIST_18(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19) >
#define LC_TYPELIST_20(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20) \
  TypeList<t1, LC_TYPELIST_19(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20) >

/**
 * TypeMap can be used to generate a whole class hierachy at compile
 * time. It is used in conjunction with a typelist and a class Unit
 * which takes a single template argument. An object of the TypeMap
 * inherits from all classes Unit<T> where T is a member of the
 * supplied typelist.
 */
  template<class TList, template <class> class Unit> class TypeMap;

/**
 * Specialization 1: Inherit from TypeMap generated from the
 * elements of the typelist.
 */
  template <class T1, class T2, template <class> class Unit>
  class TypeMap<TypeList<T1, T2>, Unit> :
        public TypeMap<T1, Unit> , public TypeMap<T2, Unit> 
  {
    public:
      template <typename T> struct Rebind 
      {
          typedef Unit<T> Result;
      };
  };

/**
 * Specialization 2: Inherit from Unit<AtomicType>
 */
  template <class AtomicType, template <class> class Unit>
  class TypeMap : public Unit<AtomicType> 
  {
    public:
      template <typename T> struct Rebind 
      {
          typedef Unit<T> Result;
      };
  };

/**
 * Specialization 3: For NullType do nothing
 */
  template <template <class> class Unit>
  class TypeMap<NullType, Unit> 
  {
    public:
      template <typename T> struct Rebind 
      {
          typedef Unit<T> Result;
      };
  };

/**
 * Say one has created TypeMap from a typelist and class Unit. Now,
 * given a type T, the ScatterHierField function returns a reference
 * to the class Unit<T> portion of the obj. Obj is of type TypeMap.
 *
 * @param obj Object to extract.
 * @return extracted object.
 */
  template <class T, class H>
  typename H::template Rebind<T>::Result& typeMapExtract(H& obj)
  {
    return obj;
  }

/**
 * Say one has created TypeMap from a typelist and class Unit. Now,
 * given a type T, the ScatterHierField function returns a reference
 * to the class Unit<T> portion of the obj. Obj is of type TypeMap.
 *
 * @param obj Object to extract.
 * @return extracted object.
 */
  template <class T, class H>
  const
  typename H::template Rebind<T>::Result& typeMapExtract(const H& obj)
  {
    return obj;
  }
}

#endif // LC_TYPE_LIST_H

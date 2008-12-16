/**
 * @file	lckeyvaltreeparser.cc
 *
 * @brief	Class to hold parse input file to create Lucee::KeyValTree object
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcexcept.h>
#include <lckeyvaltreeparser.h>

namespace Lucee
{
  template <typename REAL>
  KeyValTreeParser<REAL>::KeyValTreeParser(std::istream& istr, Lucee::KeyValTree& kvt)
    : fl(istr), kvt(kvt)
  {
    // initialize parsing by reading first token
    sym = fl.YYLex();
    tokenStr = fl.YYText();
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::parse()
  {
    input_file();
  }

  template <typename REAL>
  KeyValTreeParser<REAL>::KeyValTreeParser(KeyValTreeParser<REAL>& kvp, Lucee::KeyValTree& kvt)
    : fl(kvp.fl), kvt(kvt), sym(kvp.sym), tokenStr(kvp.tokenStr)
  {
  }

  template <typename REAL>
  bool
  KeyValTreeParser<REAL>::accept(int lsym)
  {
    if (lsym == sym)
    {
      // copy text of token 'sym' and advance the token stream
      tokenStr = fl.YYText();
      sym = fl.YYLex();
      return true;
    }
    return false;
  }

  template <typename REAL>
  bool
  KeyValTreeParser<REAL>::expect(int lsym)
  {
    if ( accept(lsym) )
      return true;

    Lucee::Except ex;
    ex << "Error: Unexpected symbol " << tokenStr
       << "found  on line no "
       << fl.getLineno();
    throw ex;
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_top()
  {
    expect(LC_LEFT_ANGLE);
    expect(LC_ID);
    tokens.push_back( tokenStr );
    expect(LC_RIGHT_ANGLE);
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_end()
  {
    expect(LC_LEFT_ANGLE_FRONT_SLASH);
    expect(LC_ID);
    expect(LC_RIGHT_ANGLE);
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value()
  {
    if (accept(LC_INT))
    {
      tokens.push_back( fl.getInteger() );
      if (!isParsingList)
        kvt.add<int>(idName, fl.getInteger());
    }
    else if (accept(LC_REAL))
    {
      tokens.push_back( fl.getReal() );
      if (!isParsingList)
        kvt.add<REAL>(idName, fl.getReal());
    }
    else if (accept(LC_ID))
    {
      tokens.push_back( tokenStr );
      if (!isParsingList)
        kvt.add<std::string>(idName, tokenStr);
    }
    else if (accept(LC_VALUE))
    {
      tokens.push_back( tokenStr );
      if (!isParsingList)
        kvt.add<std::string>(idName, tokenStr);
    }
    else if (accept(LC_STRING))
    {
      tokens.push_back( tokenStr );
      if (!isParsingList)
        kvt.add<std::string>(idName, tokenStr);
    }
    else
    {
      Lucee::Except ex;
      ex << "Error: Unknown token " << tokenStr << " encountered on line no "
         << fl.getLineno() << std::endl;
      throw ex;
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value_sub_list_rest()
  {
    if ( accept(LC_COMMA) )
    {
      value();
      // add value ...
      valueList.push_back( tokens.back() );
      tokens.pop_back(); // .. and remove it
      value_sub_list_rest();
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value_sub_list()
  {
    value();
    // add value ...
    valueList.push_back( tokens.back() );
    tokens.pop_back(); // .. and remove it
    value_sub_list_rest();
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value_list()
  {
    if (accept(LC_LEFT_BOX))
    {
      isParsingList = true;
      // clear valueList to insert list
      valueList.erase( valueList.begin(), valueList.end() );
      value_sub_list();
      expect(LC_RIGHT_BOX);
      // add it to parsed data
      if (valueList[0].type() == typeid(int))
      {
        std::vector<int> ivals;
        for (unsigned i=0; i<valueList.size(); ++i)
          ivals.push_back( Lucee::any_cast<int>(valueList[i]) );
        kvt.add(idName, ivals);
      }
      else if (valueList[0].type() == typeid(REAL))
      {
        std::vector<REAL> rvals;
        for (unsigned i=0; i<valueList.size(); ++i)
          rvals.push_back( Lucee::any_cast<REAL>(valueList[i]) );
        kvt.add(idName, rvals);
      }
      else if (valueList[0].type() == typeid(std::string))
      {
        std::vector<std::string> svals;
        for (unsigned i=0; i<valueList.size(); ++i)
          svals.push_back( Lucee::any_cast<std::string>(valueList[i]) );
        kvt.add(idName, svals);
      }
      else
      {
        // this can never happen
      }

      tokens.push_back( valueList );
      //kvt.add(idName, valueList);
      isParsingList = false;
    }
    else
    {
      isParsingList = false;
      value();
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_value()
  {
    if (accept(LC_ID))
    {
      std::string name = tokenStr;
      idName = name;
      expect(LC_EQUAL);
      // get value
      value_list();

      // add it to cryptset
      //kvt.add(name, tokens.back());
      tokens.pop_back();
    }
    else if (sym == LC_LEFT_ANGLE)
    {
      Lucee::KeyValTree kvt_c;
      // parse it
      KeyValTreeParser<REAL> wxl(*this, kvt_c);
      wxl.decl();
      // add it to current cryptset
      kvt.addSet(kvt_c);
      // update values of _sym and _tokenStr. This may be better
      // handled if we make _sym and _tokenStr into some sort of
      // per-class (rather than per-object) variables
      sym = wxl.sym; tokenStr = wxl.tokenStr;
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_list_rest()
  {
    if ((sym == LC_ID) || (sym == LC_LEFT_ANGLE))
    {
      decl_value();
      decl_list_rest();
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_list()
  {
    decl_value();
    decl_list_rest();
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl()
  {
    // tag open
    decl_top();

    // get hold of name and set it
    kvt.setName( (tokens.back()).template to_value<std::string>() );
    tokens.pop_back();

    // inside stuff
    decl_list();

    // tag close
    decl_end();
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::input_file()
  {
    decl();
  }

  // instantiations
  template class KeyValTreeParser<float>;
  template class KeyValTreeParser<double>;
}

/**
 * @file	lckeyvaltreeparser.h
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
    : _fl(istr), _kvt(kvt)
  {
    // initialize parsing by reading first token
    _sym = _fl.YYLex();
    _tokenStr = _fl.YYText();
  }

  template <typename REAL>
  KeyValTreeParser<REAL>::KeyValTreeParser(KeyValTreeParser<REAL>& kvp, Lucee::KeyValTree& kvt)
    : _fl(kvp._fl), _kvt(kvt), _sym(kvp._sym), _tokenStr(kvp._tokenStr)
  {
  }

  template <typename REAL>
  bool
  KeyValTreeParser<REAL>::accept(int sym)
  {
    if (sym == _sym)
    {
      // copy text of token 'sym' and advance the token stream
      _tokenStr = _fl.YYText();
      _sym = _fl.YYLex();
      return true;
    }
    return false;
  }

  template <typename REAL>
  bool
  KeyValTreeParser<REAL>::expect(int sym)
  {
    if ( accept(sym) )
      return true;

    Lucee::Except ex;
    ex << "Error: Unexpected symbol " << _tokenStr
       << "found  on line no "
       << _fl.lineno();
    throw ex;
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_top()
  {
    expect(LC_LEFT_ANGLE);
    expect(LC_ID);
    _tokens.push_back( Lucee::Any(_tokenStr) );
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
      _tokens.push_back( Lucee::Any(_fl.getInteger) );
      if (!isParsingList)
        _kvt.add<int>(_idName, _fl.getInteger);
    }
    else if (accept(LC_REAL))
    {
      _tokens.push_back( Lucee::Any(_fl.getReal()) );
      if (!isParsingList)
        _kvt.add<REAL>(_idName, _fl.getReal());
    }
    else if (accept(LC_ID))
    {
      _tokens.push_back( Lucee::Any(_tokenStr) );
      if (!isParsingList)
        _kvt.add<std::string>(_idName, _tokenStr);
    }
    else if (accept(LC_VALUE))
    {
      _tokens.push_back( Lucee::Any(_tokenStr) );
      if (!isParsingList)
        _kvt.add<std::string>(_idName, _tokenStr);
    }
    else if (accept(LC_STRING))
    {
      _tokens.push_back( Lucee::Any(_tokenStr) );
      if (!isParsingList)
        _kvt.add<std::string>(_idName, _tokenStr);
    }
    else
    {
      Lucee::Except ex;
      ex << "Error: Unknown token " << _tokenStr << " encountered on line no "
         << _fl.lineno() << std::endl;
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
      _valueList.push_back( _tokens.back() );
      _tokens.pop_back(); // .. and remove it
      value_sub_list_rest();
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value_sub_list()
  {
    value();
    // add value ...
    _valueList.push_back( _tokens.back() );
    _tokens.pop_back(); // .. and remove it
    value_sub_list_rest();
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::value_list()
  {
    if (accept(LC_LEFT_BOX))
    {
      isParsingList = true;
      // clear _valueList to insert list
      _valueList.erase( _valueList.begin(), _valueList.end() );
      value_sub_list();
      expect(LC_RIGHT_BOX);
      // add it to parsed data
      _tokens.push_back( _valueList );
      _kvt.add(_idName, _valueList);
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
      std::string name = _tokenStr;
      _idName = name;
      expect(LC_EQUAL);
      // get value
      value_list();

      // add it to cryptset
      //_kvt.add(name, _tokens.back());
      _tokens.pop_back();
    }
    else if (_sym == LC_LEFT_ANGLE)
    {
      Lucee::KeyValTree _kvt_c;
      // parse it
      KeyValTreeParser<REAL> wxl(this, &_kvt_c);
      wxl.decl();
      // add it to current cryptset
      _kvt.addSet(_kvt_c);
      // update values of _sym and _tokenStr. This may be better
      // handled if we make _sym and _tokenStr into some sort of
      // per-class (rather than per-object) variables
      _sym = wxl._sym; _tokenStr = wxl._tokenStr;
    }
  }

  template <typename REAL>
  void
  KeyValTreeParser<REAL>::decl_list_rest()
  {
    if ((_sym == LC_ID) || (_sym == LC_LEFT_ANGLE))
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
    _kvt.setName( (_tokens.back()).template to_value<std::string>() );
    _tokens.pop_back();

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
}

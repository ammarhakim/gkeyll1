/**
 * @file	lckeyvaltreelexer.cc
 *
 * @brief	Class to hold a tree of key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lin includes
#include <lcexcept.h>
#include <lckeyvaltreelexer.h>

// std includes
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>

#include <ctype.h>
#include <string.h>

namespace Lucee
{
  template <typename REAL>
  KeyValTreeLexer<REAL>::KeyValTreeLexer(std::istream& is)
    : is(is), lineno(1), integer(0), real(0.0) 
  {
  }

  template <typename REAL>
  int
  KeyValTreeLexer<REAL>::YYLex() 
  {
    lastSym = yylex();
    return lastSym;
  }

  template <typename REAL>
  int
  KeyValTreeLexer<REAL>::yylex()
  {
    char c;
    // skip white spaces, new lines  and tabs
    while ((c=is.get()) == ' ' || c == '\t')
      ;

    // end-of-file
    if (c == EOF)
      return LC_DONE;
    // line continuation
    if (c == '\\')
    {
      c = is.get();
      if (c == '\n')
      {
        lineno++;
        return YYLex();
      }
    }
    // comment
    if (c == '#')
    {
      while ((c=is.get()) != '\n' && c != EOF)
        ;
      if (c=='\n')
        lineno++;
      return yylex();
    }
    // number
    if (c == '-' || c == '+' || c == '.' || isdigit(c))
    {
      yytext = ""; // clear token string
      is.putback(c);
      return scanNumber();
    }
    // name
    if (isalpha(c) || c == '_' || c == '/' || c == '~')
    {
      yytext = ""; // clear token string
      do
      {
        yytext.push_back(c);
      } while ((c=is.get()) != EOF && (isidchar(c)));
      is.putback(c);
      return LC_KEY;
    }
    // string
    if (c == '"') 
    {
      yytext = ""; // clear token string
      while ((c=is.get()) != '"')
      {
        if (c == '\n' || c == EOF) {
          Lucee::Except ex;
          ex << "Missing quote for string";
          throw ex;
        }
        yytext.push_back(this->backslash(c));
      }
      return LC_STRING;
    }
    // operators, assignment, definitions and declaration
    switch (c) 
    {
      case '>':	return LC_RIGHT_ANGLE;
      case '<':	return follow('/', LC_LEFT_ANGLE_FRONT_SLASH, LC_LEFT_ANGLE);
      case '=':	return LC_EQUAL;
      case ',':	return LC_COMMA;
      case '[':	return LC_LEFT_BOX;
      case ']':	return LC_RIGHT_BOX;
      case ';': return LC_SEMI_COLON;
      case '\n': lineno++; return yylex(); // skip over a newline
      default: return c;
    }
  }

  template <typename REAL>
  char
  KeyValTreeLexer<REAL>::backslash(char c) 
  {
    static char transtab[] = "b\bf\fn\nr\rt\t";
    if (c != '\\')
      return c;
    c = is.get();
    if (islower(c) && strchr(transtab, c))
      return strchr(transtab, c)[1];
    return c;
  }

  template <typename REAL>
  unsigned
  KeyValTreeLexer<REAL>::follow(char expect, unsigned ifyes, unsigned ifno) 
  {
    int c = is.get();
    // check next character character
    if (c == expect) return ifyes;
    // unexpected one
    is.putback(c);
    return ifno;
  }

  template <typename REAL>
  bool
  KeyValTreeLexer<REAL>::isidchar(char c)
  {
    return (isalnum(c) ||
      c == '_' || 
      c == '-' ||
      c == '.' ||
      c == '/' ||
      c == '~' ||
      c == ':');
  }

// NOTE: this function needs much improvement
  template <typename REAL>
  int
  KeyValTreeLexer<REAL>::scanNumber()
  {
    int state = 0; // start state
    bool isint = true; // true of number scanned is an integer
    char c;

    // following only recognizes if the token is a number. the actual
    // conversion is handled by atof and atoi (FIX THIS: MAIN ISSUE IS
    // THAT AN OVERFLOW MAY OCCUR AND HANDLING THIS IS TRICKY)
    while( 1 )
    {
      c = is.get();
      switch(state)
      {
        case 0:
          // +,-,.,[0-9]
          if (c == '+')
            // ignore +ve sign
            state = 1;
          else if(c == '-')
            state = 1;
          else if (isdigit(c))
            state = 1;
          else if(c == '.')
            state = 2;
          else
            // error
            state = -1;

          break;

        case 1:
          // [0-9], ., [e,E]
          if (isdigit(c))
            state = 1;
          else if (c == '.')
            state = 2;
          else if((c == 'e') || (c == 'E'))
            state = 3;
          else
            // error
            state = -1;

          break;

        case 2:
          isint = false;
          // [0-9], [e,E]
          if (isdigit(c))
            state = 1;
          else if ((c == 'e') || (c == 'E'))
            state = 3;
          else
            state = -1;

          break;

        case 3:
          isint = false;
          // [+,-], [0-9]
          if(isdigit(c))
            state = 4;
          else if( c == '+')
            state = 4;
          else if (c == '-')
            state = 4;
          else
            state = -1;

          break;

        case 4:
          isint = false;
          // [0-9]
          if(isdigit(c))
            state = 4;
          else
            state = -1;
      }
      if(state == -1)
        break; // finish parsing
      yytext.push_back(c); // push character to current token string
    }
    is.putback(c);
    // compute number
    if (isint)
    {
      integer = atoi(yytext.c_str());
      return LC_INT;
    }
    else
    {
      real = atof(yytext.c_str());
      return LC_REAL;
    }
  }

  // instantiations
  template class KeyValTreeLexer<float>;
  template class KeyValTreeLexer<double>;
}

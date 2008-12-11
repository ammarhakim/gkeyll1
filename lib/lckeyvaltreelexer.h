/**
 * @file	lckeyvaltreelexer.h
 *
 * @brief	Class to hold a tree of key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_KEYVAL_TREE_LEXER_H
#define LC_KEYVAL_TREE_LEXER_H

// std includes
#include <iostream>
#include <istream>
#include <string>

namespace Lucee
{
/** Token definitions */
  enum 
  {
    LC_ERROR = -2,
    LC_DONE = -1,
    LC_LEFT_BOX = 256,
    LC_RIGHT_BOX,
    LC_LEFT_ANGLE,
    LC_RIGHT_ANGLE,
    LC_ID,
    LC_KEY = LC_ID, // key and value are just IDs
    LC_VALUE = LC_ID,
    LC_EQUAL,
    LC_FRONT_SLASH,
    LC_COMMA,
    LC_INT,
    LC_REAL,
    LC_STRING,
    LC_SEMI_COLON,
    LC_LEFT_ANGLE_FRONT_SLASH
  };

/**
 * Lexical analysis class for tokenizing Lucee input files.
 */
  template <typename REAL>
  class KeyValTreeLexer
  {
    public:
/**
 * Set input stream to one supplied. This defaults to standard input.
 *
 * @param is Input stream to scan for characters
 */
      KeyValTreeLexer(std::istream& is)
        : _is(is), _lineno(1), _integer(0), _real(0.0) {
      }

/**
 * Scans input stream and returns a single token.
 */
      int YYLex() {
        _lastSym = yylex();
        return _lastSym;
      }

/**
 * Returns a character pointer representing the current token scanned.
 */
      std::string YYText() const {
        return _yytext;
      }

/**
 * Returns the current line number being scanned
 */
      unsigned lineno() const {
        return _lineno;
      }

/**
 * Returns integer scanned
 */
      int integer() const {
        return _integer;
      }

/**
 * Returns real number scanned
 */
      REAL real() const {
        return _real;
      }

      int lastSym() const {
        return _lastSym;
      }

    private:
      std::istream& _is; // input stream to tokenize
      unsigned _lineno;
      std::string _yytext; // text of last token read
      int _integer; // scanned integer
      REAL _real; // scanned real
      unsigned _lastSym; // last symbol scanned

      // do the actual tokenization
      int yylex();

      // get next char with \'s interpreted
      char backslash(char c);

      // look ahead for two character operators
      unsigned follow(char expect, unsigned ifyes, unsigned ifno);

      // check if character can belong to LC_ID
      bool isidchar(char c);

      // scans a number
      int scanNumber();
  };
}

#endif // LC_KEYVAL_TREE_LEXER_H

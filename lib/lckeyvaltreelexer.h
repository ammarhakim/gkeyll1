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
/** Token definitions for use in lexer */
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
 * Create a new lexer class with a input stream.
 *
 * @param is Input stream to scan for characters
 */
      KeyValTreeLexer(std::istream& is)
        : _is(is), _lineno(1), _integer(0), _real(0.0) {
      }

/**
 * Scans input stream and returns a single token.
 *
 * @return token of last symbol read.
 */
      int YYLex() {
        _lastSym = yylex();
        return _lastSym;
      }

/**
 * String representing last token scanned.
 *
 * @return string representation of last token read.
 */
      std::string YYText() const {
        return _yytext;
      }

/**
 * Current line number being scanned.
 *
 * @return line number being scanned.
 */
      unsigned lineno() const {
        return _lineno;
      }

/**
 * Last integer scanned
 *
 * @return last integer scanned.
 */
      int integer() const {
        return _integer;
      }

/**
 * Last real number scanned.
 *
 * @return last real number scanned.
 */
      REAL real() const {
        return _real;
      }

/**
 * Last symbol scanned.
 *
 * @return last symbol scanned.
 */
      int lastSym() const {
        return _lastSym;
      }

    private:
/** Input stream of characters */
      std::istream& _is;
/** Current line number */
      unsigned _lineno;
/** Text for the last token read */
      std::string _yytext;
/** Last integer read */
      int _integer;
/** Last real scanned */
      REAL _real;
/** Last symbols scanned */
      unsigned _lastSym;

/**
 * Tokenize and return token scanned 
 *
 * @return token scanned
 */
      int yylex();

/**
 * Get next char with \'s interpreted
 *
 * @return next character
 */
      char backslash(char c);

/**
 * Look ahead for two character operators
 *
 * @param expect character to expect
 * @param ifyes token to return if next character is 'expect'
 * @param ifno token to return if next character is not 'expect'
 * @return either 'ifyes' if next character is 'expect', 'ifno' otherwise
 */
      unsigned follow(char expect, unsigned ifyes, unsigned ifno);

/**
 * Is character an ID character?
 *
 * @param c Character to check
 * @return true if character is ID character, false otherwise.
 */
      bool isidchar(char c);

/**
 * Scan stream for a number. Does not return the number, but this can
 * be gotten by the integer() or real() methods, depending on the
 * token returned.
 *
 * @return token, either LC_INT or LC_REAL
 */
      int scanNumber();
  };
}

#endif // LC_KEYVAL_TREE_LEXER_H

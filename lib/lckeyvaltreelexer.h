/**
 * @file	lckeyvaltreelexer.h
 *
 * @brief	Class to tokenize input file.
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
      KeyValTreeLexer(std::istream& is);

/**
 * Scans input stream and returns a single token.
 *
 * @return token of last symbol read.
 */
      int YYLex();

/**
 * String representing last token scanned.
 *
 * @return string representation of last token read.
 */
      std::string YYText() const {
        return yytext;
      }

/**
 * Current line number being scanned.
 *
 * @return line number being scanned.
 */
      unsigned getLineno() const {
        return lineno;
      }

/**
 * Last integer scanned
 *
 * @return last integer scanned.
 */
      int getInteger() const {
        return integer;
      }

/**
 * Last real number scanned.
 *
 * @return last real number scanned.
 */
      REAL getReal() const {
        return real;
      }

/**
 * Last symbol scanned.
 *
 * @return last symbol scanned.
 */
      int getLastSym() const {
        return lastSym;
      }

    private:
/** Input stream of characters */
      std::istream& is;
/** Current line number */
      unsigned lineno;
/** Text for the last token read */
      std::string yytext;
/** Last integer read */
      int integer;
/** Last real scanned */
      REAL real;
/** Last symbols scanned */
      unsigned lastSym;

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

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

#ifndef LC_KEYVAL_TREE_PARSER_H
#define LC_KEYVAL_TREE_PARSER_H

// lib includes
#include <lcany.h>
#include <lckeyvaltree.h>
#include <lckeyvaltreelexer.h>

// std includes
#include <iostream>
#include <string>

namespace Lucee
{

/**
 * Class to parse input files and create a Lucee::KeyValTree
 * corresponding to the input file.
 */
  template <typename REAL>
  class KeyValTreeParser
  {
    public:
/**
 * Create a new parser, which reads from an input stream. As the
 * parser reads from the stream it updates the KeyValTree supplied.
 *
 * @param istr Input stream for parsing
 * @param kvt KeyValTree object to update as the stream is parsed
 */
      KeyValTreeParser(std::istream& istr, Lucee::KeyValTree& kvt);

/**
 * Parse input file and initialize KeyValTree set at construction
 * time.
 */
      void parse();

    private:
/**
 * Internal constructor
 */
      KeyValTreeParser(KeyValTreeParser<REAL>& kvp, Lucee::KeyValTree& kvt);

/**
 * Accept supplied symbol.
 *
 * @param sym Symbol to accept.
 * @return true if symbol was accepted, false otherwise.
 */
      bool accept(int sym);

/**
 * Expect supplied symbol.
 *
 * @param sym Symbol to expect.
 * @return true if expected symbol was found, false otherwise.
 */
      bool expect(int sym);

/** Parse non-terminal in grammar */
      void input_file();

/** Parse non-terminal in grammar */
      void decl();

/** Parse non-terminal in grammar */
      void decl_top();

/** Parse non-terminal in grammar */
      void decl_end();

/** Parse non-terminal in grammar */
      void decl_list();

/** Parse non-terminal in grammar */
      void decl_list_rest();

/** Parse non-terminal in grammar */
      void decl_value();

/** Parse non-terminal in grammar */
      void value();

/** Parse non-terminal in grammar */
      void value_list();

/** Parse non-terminal in grammar */
      void value_sub_list();

/** Parse non-terminal in grammar */
      void value_sub_list_rest();
      
/** Lexer to get tokens from stream */
      Lucee::KeyValTreeLexer<REAL> _fl;
/** Reference to tree being modified */
      Lucee::KeyValTree& _kvt; // 
/** Current symbol read */
      int _sym;
/** Current token string */
      std::string _tokenStr;
/** List of tokens */
      std::vector<Lucee::Any> _tokens;
/** Values read from a list */
      std::vector<Lucee::Any> _valueList;
/** Name os key being parser */
      std::string _idName;
/** Flag to indicate is we are parsing a list */
      bool isParsingList;
  };
}

#endif // LC_KEYVAL_TREE_PARSER_H

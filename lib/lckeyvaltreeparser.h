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

      // symbol lookup functions
      bool accept(int sym);
      bool expect(int sym);

      // non-terminals 
      void input_file();
      void decl();
      void decl_top();
      void decl_end();
      void decl_list();
      void decl_list_rest();
      void decl_value();
      void value();
      void value_list();
      void value_sub_list();
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

/*
 * tokenBySpace.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef TOKENBYSPACE_H_
#define TOKENBYSPACE_H_

#include<iostream>
#include<vector>
#include<string>
#include<sstream>

using namespace std;

vector<string>
tokenBySpace( const string& str );

int
string2int( const string& str );

//! Break a string (supplied as the second argument) into tokens, returned
//! in the first argument. Tokens are determined by the delimiters supplied
//! (defaults to whitespace (i.e., spaces, tabs, newlines)
bool tokenize(std::vector<std::string> &vcr, const char *buf,
                    const char *delimstr);

#endif /* TOKENBYSPACE_H_ */

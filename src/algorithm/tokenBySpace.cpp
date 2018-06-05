/*
 * tokenBySpace.cpp
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#include"tokenBySpace.h"

vector<string>
tokenBySpace( const string& str ){
	if( str.empty() ){
		string						s("");
		vector<string>		vs;
		vs.push_back( s );
		return vs;
	}
	string buf ="";  // Have a buffer string

	stringstream ss(str);  // Insert the string into a stream
	vector<string> tokens;  // Create vector to hold our words
	while (ss >> buf){
		tokens.push_back(buf);
		//     cout<<"buf :"<<buf<<endl;
	}
	if( tokens.empty() ){
		string						s("");
		vector<string>		vs;
		vs.push_back( s );
		return vs;
	}
	buf ="";
	buf.clear();
	return tokens;
}

int
string2int( const string& str ){
	int b;
	stringstream 	ss(str);
	ss>>b;
	return b;
}

//! Break a string (supplied as the second argument) into tokens, returned
//! in the first argument. Tokens are determined by the delimiters supplied
//! (defaults to whitespace (i.e., spaces, tabs, newlines)
bool tokenize(std::vector<std::string> &vcr, const char *buf,
                    const char *delimstr)
{
  vcr.clear();
  if (!buf || !delimstr)
    return false;

  string s(buf);
  s += delimstr[0]; //forces last token to be parsed
  size_t startpos=0,endpos=0;

  for (;;)
    {
      startpos = s.find_first_not_of(delimstr,startpos);
      endpos = s.find_first_of(delimstr,startpos);

      if (endpos <= s.size() && startpos <= s.size())
        vcr.push_back(s.substr(startpos,endpos-startpos));
      else
        break;

      startpos = endpos+1;
    }

  return(true);
}

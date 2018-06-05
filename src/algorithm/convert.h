/*
 * convert.h
 *
 *  Created on: Jul 28, 2014
 *      Author: stan
 */

#ifndef CONVERT_H_
#define CONVERT_H_

#include <sstream>
std::string float2string(float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}


#endif /* CONVERT_H_ */

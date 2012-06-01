/*
 * File:   utils.h
 * Author: mvyvermn
 *
 * Created on 27 januari 2012, 16:16
 */

#ifndef UTILS_H
#define	UTILS_H

#include <string>

namespace Utils{

    std::string convertInt(const int number);
    bool contains(const std::string & str,const long begin,const long end,const char c);

}

#endif	/* UTILS_H */


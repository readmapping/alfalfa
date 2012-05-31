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

    std::string convertInt(int number);
    bool contains(std::string & str, long begin, long end, char c);

}

#endif	/* UTILS_H */


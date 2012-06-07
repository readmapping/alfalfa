/*
 * Copyright 2012, Michael Vyverman <michael.vyverman@ugent.be>
 *
 * This file is part of ALFALFA.
 *
 * ALFALFA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ALFALFA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ALFALFA.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "utils.h"
#include <sstream>
#include <assert.h>

namespace Utils{

    std::string convertInt(const int number){
       std::stringstream ss;//create a stringstream
       ss << number;//add number to the stream
       return ss.str();//return a string with the contents of the stream
    }

    bool contains(const std::string & str,const long begin,const long end, const char c){
        bool contains = false;
        long i = begin;
        assert(begin >= 0 && end < str.length());
        while(!contains && i <= end){
            contains = str[i] == c;
            i++;
        }
        return contains;
    }

}

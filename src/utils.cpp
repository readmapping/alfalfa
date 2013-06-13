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
#include <algorithm>

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
    
    // Return the reverse complement of sequence. This allows searching
    // the plus strand of instances on the minus strand.
    void reverse_complement(std::string &seq_rc, bool nucleotides_only) {//TODO optimize this for char arrays
        // Reverse in-place.
        std::reverse(seq_rc.begin(), seq_rc.end());
        for(long i = 0; i < (long)seq_rc.length(); i++) {
            // Adapted from Kurtz code in MUMmer v3.
            switch(seq_rc[i]) {
                case 'a': seq_rc[i] = 't'; break;
                case 'c': seq_rc[i] = 'g'; break;
                case 'g': seq_rc[i] = 'c'; break;
                case 't': seq_rc[i] = 'a'; break;
                case 'r': seq_rc[i] = 'y'; break; /* a or g */
                case 'y': seq_rc[i] = 'r'; break; /* c or t */
                case 's': seq_rc[i] = 's'; break; /* c or g */
                case 'w': seq_rc[i] = 'w'; break; /* a or t */
                case 'm': seq_rc[i] = 'k'; break; /* a or c */
                case 'k': seq_rc[i] = 'm'; break; /* g or t */
                case 'b': seq_rc[i] = 'v'; break; /* c, g or t */
                case 'd': seq_rc[i] = 'h'; break; /* a, g or t */
                case 'h': seq_rc[i] = 'd'; break; /* a, c or t */
                case 'v': seq_rc[i] = 'b'; break; /* a, c or g */
                default:
                if(!nucleotides_only) seq_rc[i] = 'n';
                break; /* anything */
            }
        }
    }

}

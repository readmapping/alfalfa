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

#ifndef UTILS_H
#define	UTILS_H

#include <string>
#include <stdint.h>

namespace Utils{
    
    static const char RC_TABLE[256] = {  'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//0-9
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//10-19
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//20-29
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//30-39
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//40-49
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//50-59
                                        'n', 'n', 'n', 'n', 'n', 'T', 'V', 'G', 'H', 'N',//60-69 65:A, 67:C
                                        'N', 'C', 'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N',//70-79 71:G
                                        'N', 'N', 'Y', 'S', 'A', 'N', 'B', 'W', 'N', 'R',//80-89 84:T
                                        'N', 'N', 'N', 'N', 'N', 'N', 'N', 't', 'v', 'g',//90-99 97:a, 99: c
                                        'h', 'n', 'n', 'c', 'd', 'n', 'n', 'm', 'n', 'k',//100-109 103:g
                                        'n', 'n', 'n', 'n', 'y', 's', 'a', 'n', 'b', 'w',//110-119 116:t
                                        'n', 'r', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//120-129
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//130-139
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//140-149
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//150-159
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//160-169
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//170-179
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//180-189
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//190-199
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//200-209
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//210-219
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//220-229
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//230-239
                                        'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n',//240-249
                                        'n', 'n', 'n', 'n', 'n', 'n' };//250-255
    
static const uint8_t   ORDVALUE[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//0-9
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//10-19
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//20-29
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//30-39
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//40-49
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//50-59
                                        4, 4, 4, 4, 4, 0, 4, 1, 4, 4,//60-69 65:A, 67:C
                                        4, 2, 4, 4, 4, 4, 4, 4, 4, 4,//70-79 71:G
                                        4, 4, 4, 4, 3, 4, 4, 4, 4, 4,//80-89 84:T
                                        4, 4, 4, 4, 4, 4, 4, 0, 4, 1,//90-99 97:a, 99: c
                                        4, 4, 4, 2, 4, 4, 4, 4, 4, 4,//100-109 103:g
                                        4, 4, 4, 4, 4, 4, 3, 4, 4, 4,//110-119 116:t
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//120-129
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//130-139
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//140-149
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//150-159
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//160-169
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//170-179
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//180-189
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//190-199
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//200-209
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//210-219
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//220-229
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//230-239
                                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//240-249
                                        4, 4, 4, 4, 4, 4 };//250-255

    std::string convertInt(const int number);//not used anymore
    bool contains(const std::string & str,const long begin,const long end,const char c);//not used anymore
    void reverse_complement(std::string &seq_rc, bool nucleotides_only);

}

#endif	/* UTILS_H */


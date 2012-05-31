
#include "utils.h"
#include <sstream>
#include <assert.h>

namespace Utils{

    std::string convertInt(int number){
       std::stringstream ss;//create a stringstream
       ss << number;//add number to the stream
       return ss.str();//return a string with the contents of the stream
    }

    bool contains(std::string & str, long begin, long end, char c){
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


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

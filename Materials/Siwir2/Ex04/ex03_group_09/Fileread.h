#pragma once
#include <map>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <stdexcept>
 
class FileReader{
public:
    template< typename Type>
    Type getParameter(const std::string& key) const{
        Type val;
        std::map<std::string,std::string>::const_iterator it;
        it = mymap.find(key);
        if(it == mymap.end()){
            std::cerr<<"Key does not exist "<<std::endl;
            throw std::runtime_error("Key does not exist ");
        }
        else{
            std::stringstream str;
            str<< (*it).second; 
            str>> val;
            return val;
        }
    }
    void readParameterFile(std::ifstream& i);
 
private:
    std::map<std::string,std::string> mymap;
}; 

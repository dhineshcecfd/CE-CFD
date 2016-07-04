#include "Fileread.h"
 
void FileReader::readParameterFile(std::ifstream& i){
    std::string temp1, temp2;
    while(!i.eof()){
        i>>temp1;
        i>>temp2;
        mymap[temp1]= temp2;
    }
}

#include <iostream>
#include "utils.h"

using namespace std;

Utils::Utils()
{

}

string Utils::strToUpper(string s)
{
    string saida;
    for(int i=0;i<s.size();i++){
        saida+=toupper(s[i]);
    }
    return saida;
}


vector<string> Utils::splitString(string str, char delimitador)
{
    vector<string> tokens;
    string token = "";
    for(unsigned int i = 0; i < str.size(); i++){
        if(str.at(i) == delimitador){
            tokens.push_back(token);
            token = "";
        }else{
            token.push_back(str.at(i));
        }
        if(i == str.size()-1){
            tokens.push_back(token);
        }
    }
    return tokens;
}

void Utils::replaceAll(string& str, string from, string to) {
    size_t start_pos = str.find(from);
    while(start_pos != string::npos){
        str.replace(start_pos, from.length(), to);
        start_pos = str.find(from);
    }

}

bool Utils::replace(string& str, string from, string to) {
    size_t start_pos = str.find(from);
    if(start_pos == string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

#define EPSILON 1E-6
bool Utils::doubleEqual(double a, double b) {
    return fabs(a - b) < EPSILON;
}

bool Utils::isBinary(double a) {
    return Utils::doubleEqual(a, 1.0) || Utils::doubleEqual(a, 0.0);
}

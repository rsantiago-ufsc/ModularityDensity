#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;
class Utils
{
public:
    Utils();
    static string strToUpper(string s);
    static vector<string> splitString(string str, char delimitador);
    static bool replace(string& str, string from, string to);
    static void replaceAll(string& str, string from, string to) ;
    static bool doubleEqual(double a, double b);
    static bool isBinary(double a);
};

#include "utils.cpp"

#endif // UTILS_H

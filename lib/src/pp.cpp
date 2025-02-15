#include "pp.h"

/** Pretty print for testing purposes */
bool pp_mon(stringstream& ss, zz_p c, long dx, bool p){
    if (c == 0) return false;
    if (p) ss << " + ";

    bool pLast = true;
    if (c != 1) ss << c;
    pLast = false;
    
    if (dx != 0) {
        if (pLast) ss << "*";
        pLast = true;
        
        ss << "x";
    }
    if (dx > 1) ss << "^" << dx;

    return true;
}

string pp(const zz_pX& a){
    if (a == 0)
        return "0";
    
    stringstream ss;
    bool p = false;
    
    for (long j = 0; j < a.rep.length(); j++)
        p |= pp_mon(ss, a.rep[j], j, p);

    return ss.str();
}

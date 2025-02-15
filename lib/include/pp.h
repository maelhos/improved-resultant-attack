
#ifndef PP__H
#define PP__H

#include <NTL/lzz_pX.h>
#include <sstream>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <string>

NTL_CLIENT

inline string pp(const zz_p& a)
{
    return to_string(a._zz_p__rep);
}

template <typename T>
string pp(const Vec<T>& v)
{
    string out("[");
    for (long i = 0; i < v.length(); i++)
    {
        out += pp(v[i]);
        if (i != v.length() - 1) out += ", ";
    }
    return out + "]";
}


template <typename T>
inline string pp(const Mat<T>& v)
{
    return pp(v._mat__rep);
}

string pp(const zz_pX& a);


#endif
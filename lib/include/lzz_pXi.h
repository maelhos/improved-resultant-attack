#ifndef LZZ_PXI__H
#define LZZ_PXI__H

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <utility>
#include <array>
#include <sstream>

NTL_CLIENT

/*------------------------------------------------------------*/
/* Class for Multivariate polynomials                         */
/*------------------------------------------------------------*/

template <auto N>
class zz_pXi{
public:
    Vec<zz_pXi<N-1>> rep;

    zz_pXi() {};
    zz_pXi(Vec<zz_pXi<N-1>>& r) : rep(r) { normalize(); };
    zz_pXi(zz_pXi<N - 1>& r) : rep{}
    {
        rep.SetLength(1);
        rep[0] = r;
    }

    long degXn() const { return rep.length() - 1;}
    bool is_zero() const { return rep.length() == 0; }
    void zero() { rep.SetLength(0); }
    void one() { zz_pXi<N-1> one{}; one.one(); rep.SetLength(1); rep[0] = one; }
    void Xn() { zz_pXi<N-1> one{}; one.one(); rep.SetLength(2); rep[1] = one; }
    void normalize()
    {
        long idx = rep.length() - 1;
        while (idx >= 0 && rep[idx].is_zero())
            idx--;
        rep.SetLength(idx+1);
    }

    zz_pXi<N + 1> extend() { return zz_pXi<N + 1>(*this); }
};

template <>
class zz_pXi<1> : public zz_pX {
public:
    using zz_pX::zz_pX;
    zz_pXi(zz_pX r) : zz_pX(r) {};
    long degXn() const { return rep.length() - 1;}
    bool is_zero() const { return rep.length() == 0; }
    void zero() { rep.SetLength(0); }
    void one() { rep.SetLength(1); rep[0] = 1; }
    void Xn() { rep.SetLength(2); rep[0] = 0; rep[1] = 1; }
    void normalize()
    {
        long idx = rep.length() - 1;
        while (idx >= 0 && (rep[idx] == 0))
            idx--;
        rep.SetLength(idx+1);
    }
    zz_pXi<2> extend() { return zz_pXi<2>(*this); }
};

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/

template <auto N>
ostream& operator<<(ostream& s, const zz_pXi<N>& a)
{
    return s << a.rep;
}

template <auto N>
istream& operator>>(istream& s, zz_pXi<N>& a)
{
    return s >> a.rep;
}

/*------------------------------------------------------------*/
/* Pretty print for debug purposes                            */
/*------------------------------------------------------------*/
/** Pretty print for testing purposes */
template <auto N>
bool pp_mon(stringstream& ss, zz_p c, array<long, N>& di, bool p)
{
    if (c == 0) return false;
    if (p) ss << " + ";

    bool pLast = false;
    if (c != 1) 
    {
        ss << c;
        pLast = true;
    }

    for (unsigned int i = 0; i < N; i++)
    {
        long dxi = di[i];
        if (dxi != 0) 
        {
            if (pLast) ss << "*";
            pLast = true;
            ss << "x" << i + 1;
        }
        if (dxi > 1) ss << "^" << dxi;
    }

    return true;
}

template <auto N, auto K>
void rec_mon(stringstream& ss, const zz_pXi<K>& a, array<long, N>& di, bool& p)
{
    if constexpr (K == 1)
    {
        for (long i = 0; i < a.rep.length(); i++)
        {
            di[0] = i;
            p |= pp_mon<N>(ss, a.rep[i], di, p);
        }
    }
    else
    {
        for (long i = 0; i < a.rep.length(); i++)
        {
            zz_pXi<K - 1> ap{};
            ap.rep = move(a.rep[i].rep);
            di[K - 1] = i;
            //cout << ap << " : " << di << endl;
            rec_mon<N, K - 1>(ss, ap, di, p);
        }
    }
    
}

template <auto N>
string pp(const zz_pXi<N>& a)
{
    if (a.is_zero())
        return "0";
    
    stringstream ss;
    bool p = false;
    array<long, N> current_mon{};

    rec_mon(ss, a, current_mon, p);
    return ss.str();
}

/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
template <auto N>
inline void add(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    long dyA = a.degXn();
    long dyB = b.degXn();
    long dy = max(dyA, dyB);
    long dy0 = min(dyA, dyB);
    c.rep.SetLength(dy + 1);
    for (long i = 0; i <= dy0; i++)
        add(c.rep[i], a.rep[i], b.rep[i]);

    for (long i = dy0 + 1; i <= dyB; i++)
        c.rep[i] = b.rep[i];

    for (long i = dy0 + 1; i <= dyA; i++)
        c.rep[i] = a.rep[i];

    c.normalize();
}

/** Returns `a` + `b` */
template <auto N>
inline zz_pXi<N> add(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    zz_pXi<N> c;
    add(c, a, b);
    return c;
}

/** Returns `a` + `b` */
template <auto N>
inline zz_pXi<N> operator+(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    return add(a, b);
}

template <auto N>
inline zz_pXi<N> operator+=(zz_pXi<N>& a, const zz_pXi<N>& b)
{
    add(a, a, b);
    return a;
}

/*------------------------------------------------------------*/
/* negation                                                   */
/*------------------------------------------------------------*/
template <auto N>
inline void neg(zz_pXi<N>& c, const zz_pXi<N>& a){
    c.rep.SetLength(a.rep.length());
    for (long i = 0; i < a.rep.length(); i++)
        neg(c.rep[i], a.rep[i]);
}

template <>
inline void neg(zz_pXi<1>& c, const zz_pXi<1>& a){
    c.rep.SetLength(a.rep.length());
    for (long i = 0; i < a.rep.length(); i++)
        c.rep[i] = -a.rep[i];
}

template <auto N>
inline zz_pXi<N> neg(const zz_pXi<N>& a){
    zz_pXi<N> c;
    neg(c, a);
    return c;
}

template <auto N>
inline zz_pXi<N> operator-(const zz_pXi<N>& a)
{
   return neg(a);
}

/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
template <auto N>
inline void sub(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    long dyA = a.degXn();
    long dyB = b.degXn();
    long dy = max(dyA, dyB);
    long dy0 = min(dyA, dyB);
    c.rep.SetLength(dy + 1);
    for (long i = 0; i <= dy0; i++)
        sub(c.rep[i], a.rep[i], b.rep[i]);
        
    for (long i = dy0 + 1; i <= dyB; i++)
        c.rep[i] = -b.rep[i];

    for (long i = dy0 + 1; i <= dyA; i++)
        c.rep[i] = a.rep[i];

    c.normalize();
}

/** Returns `a` - `b` */
template <auto N>
inline zz_pXi<N> sub(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    zz_pXi<N> c;
    sub(c, a, b);
    return c;
}

/** Returns `a` - `b` */
template <auto N>
inline zz_pXi<N> operator-(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    return sub(a, b);
}

/*------------------------------------------------------------*/
/* equality test                                              */
/*------------------------------------------------------------*/

/**  Returns non-zero if `a` == `b`, 0 otherwise  */ 
template <auto N>
inline long operator==(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
   return a.rep == b.rep;
}

/**  Returns non-zero if `a` != `b`, 0 otherwise  */ 
template <auto N>
inline long operator!=(const zz_pXi<N>& a, const zz_pXi<N>& b)
{ 
    return !(a == b); 
}

/*------------------------------------------------------------*/
/* random element f with deg(f,x_i) < d and for i = 1..N      */
/*------------------------------------------------------------*/
template <auto N>
inline void random_zz_pXi(zz_pXi<N>& f, long d)
{
    f.rep.SetLength(d);
    for (long i = 0; i < d; i++)
        random_zz_pXi(f.rep[i], d);
}

template <>
inline void random_zz_pXi(zz_pXi<1>& f, long d)
{
    f.rep.SetLength(d);
    random(f.rep, d);
}

template <auto N>
inline zz_pXi<N> random_zz_pXi(long d){
    zz_pXi<N> ret;
    random_zz_pXi(ret, d);
    return ret;
}


/*------------------------------------------------------------*/
/* degrees                                                    */
/*------------------------------------------------------------*/
template <auto K, auto N>
inline void __degrees_rec(const zz_pXi<K>& f, array<long, N>& arr)
{
    if constexpr (K == 1)
        arr[0] = max(arr[0], deg(f));
    else
    {
        arr[K - 1] = max(arr[K - 1], f.rep.length() - 1);
        for (long i = 0; i < f.rep.length(); i++)
            __degrees_rec(f.rep[i], arr);
    }
}

template <auto N>
inline array<long, N> degrees(const zz_pXi<N>& f)
{
    array<long, N> ret{};
    __degrees_rec(f, ret);
    return ret;
}

template <auto N>
inline long prodp1(const array<long, N>& degs)
{
    long p = 1;
    for (long i = 0; i < N; i++)
        p *= degs[i] + 1;
    
    return p;    
}

/*------------------------------------------------------------*/
/* to kronecker substitution                                  */
/*------------------------------------------------------------*/
template <auto K, auto N>
inline void rec_to_kro(zz_p*& out_e, const Vec<zz_pXi<K>>& b, const array<long, N>& degs, long prod_p)
{
    if constexpr (K == 1)
    {
        for (long i = 0; i < b.length(); i++)
        {
            long d = deg(b[i]);
            if (d >= 0)
            {
                const zz_p * be = b[i].rep.elts();
                for (long j = 0; j <= d; j++)
                    out_e[j] = be[j];
            }
            out_e += degs[0] + 1;
        }
    }
    else
    {
        long np = prod_p / (degs[K - 1] + 1);
        for (long k = 0; k < b.length(); k++)
        {
            Vec<zz_pXi<K - 1>> a = move(b[k].rep);
            long dznm1 = a.length() - 1;
            rec_to_kro<K - 1, N>(out_e, a, degs, np);

            if (dznm1 + 1 <= degs[K - 1]) out_e += np*(degs[K - 1] - dznm1);
        }
    }
}

template <auto N>
void to_kronecker(zz_pX& out, const zz_pXi<N>& b, const array<long, N - 1>& degs)
{
    if (b.is_zero())
    {
        out = 0;
        return;
    }

    long dXn = b.degXn();
    long deg_p = prodp1<N - 1>(degs);
    out.rep.SetLength(deg_p * (dXn + 1));
    zz_p* out_e = out.rep.elts();
    rec_to_kro<N - 1, N - 1>(out_e, b.rep, degs, deg_p);
    out.normalize();
}

/*------------------------------------------------------------*/
/* from kronecker substitution                                */
/*------------------------------------------------------------*/
template <auto K, auto N>
inline void rec_from_kro(Vec<zz_pXi<K>>& out, const zz_pX& a, long& idx, const array<long, N>& degs)
{
    if constexpr (K == 1)
    {
        const zz_p * ae = a.rep.elts();
        for (long i = 0; i < out.length() ; i++)
        {
            out[i].rep.SetLength(degs[0] + 1);
            zz_p * ce = out[i].rep.elts();
            long j;
            for (j = 0; (j <= degs[0]) && (idx <= deg(a)); j++)
                ce[j] = ae[idx++];
            out[i].normalize();
        }
    }   
    else
    {
        for (long k = 0; k < out.length(); k++)
        {
            out[k].rep.SetLength(degs[K - 1] + 1);
            rec_from_kro<K - 1, N>(out[k].rep, a, idx, degs);
            out[k].normalize();
        }
    }
    
}

template <auto N>
void from_kronecker(zz_pXi<N>& out, const zz_pX& a, const array<long, N - 1>& degs)
{
    if (a == 0)
    {
        out.rep.SetLength(0);
        return;
    }

    long deg_p = prodp1<N - 1>(degs);
    long dXn = deg(a) / deg_p;

    long idx = 0;
    out.rep.SetLength(dXn + 1);
    rec_from_kro<N - 1, N - 1>(out.rep, a, idx, degs);
    out.normalize();
}

/*------------------------------------------------------------*/
/* multiplication                                             */
/*------------------------------------------------------------*/

template <auto N>
inline void mul_naive(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    if (a.is_zero() || b.is_zero())
    {
        c.zero();
        return;
    }

    zz_pXi<N> c_tmp;
    long da = a.degXn();
    long db = b.degXn();

    c_tmp.rep.SetLength(da + db + 1, zz_pXi<N-1>{});
    for (long i = 0; i <= da; i++)
    {
        for (long j = 0; j <= db; j++)
        {
            zz_pXi<N-1> prod_tmp;
            mul_naive(prod_tmp, a.rep[i], b.rep[j]);
            add(c_tmp.rep[i + j], c_tmp.rep[i + j], prod_tmp);
        }
    }

    c_tmp.normalize();
    c = c_tmp;
}

template <>
inline void mul_naive<1>(zz_pXi<1>& c, const zz_pXi<1>& a, const zz_pXi<1>& b)
{
    zz_pX ret{};
    mul(ret, zz_pX(a), zz_pX(b));
    c = ret;
}

/** Returns `a` * `b`, computed using the naive algorithm */
template <auto N>
inline zz_pXi<N> mul_naive(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    zz_pXi<N> c;
    mul_naive(c, a, b);
    return c;
}

/*------------------------------------------------------------*/
/* kronecker multiplication                                   */
/*------------------------------------------------------------*/

/** Returns `a` * `b`, computed using Kronecker substitution */
template <auto N>
inline zz_pXi<N> mul_kronecker(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    zz_pXi<N> c;
    mul_kronecker(c, a, b);
    return c;
}

template <auto N>
void mul_kronecker(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    if (&c == &a || &c == &b)
    {
        c = mul_kronecker(a, b);
        return;
    }

    if (a.is_zero() || b.is_zero())
    {
        c.zero();
        return;
    }
    
    array<long, N> degsa = degrees(a);
    array<long, N> degsb = degrees(b);

    array<long, N - 1> kronecker_basis{};
    for (long i = 0; i < N - 1; i++)
        kronecker_basis[i] = degsa[i] + degsb[i];
    

    zz_pX apX, bpX, cpX;

    to_kronecker(apX, a, kronecker_basis);
    to_kronecker(bpX, b, kronecker_basis);
    cpX = apX * bpX;

    from_kronecker(c, cpX, kronecker_basis);
}

template <auto N>
void sqr(zz_pXi<N>& c, const zz_pXi<N>& a)
{
    if (&c == &a)
    {
        zz_pXi<N> ret{};
        sqr(ret, a);
        c = ret;
        return;
    }

    if (a.is_zero())
    {
        c.zero();
        return;
    }
    
    array<long, N> degsa = degrees(a);

    array<long, N - 1> kronecker_basis{};
    for (long i = 0; i < N - 1; i++)
        kronecker_basis[i] = 2*degsa[i];
    
    zz_pX apX, cpX;

    to_kronecker(apX, a, kronecker_basis);
    sqr(cpX, apX);
    from_kronecker(c, cpX, kronecker_basis);
}

template <>
inline void sqr(zz_pXi<1>& c, const zz_pXi<1>& a)
{
    zz_pX ax(a);
    zz_pX cx{};
    sqr(cx, ax);
    
    c = zz_pXi<1>(cx);
}

template <>
inline void mul_kronecker(zz_pXi<1>& c, const zz_pXi<1>& a, const zz_pXi<1>& b)
{
    zz_pX ret{};
    mul(ret, zz_pX(a), zz_pX(b));
    c = ret;
}

template <auto N>
inline void mul(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    mul_kronecker(c, a, b);
}

/** Returns `a` * `b`. For the time being, we only use
 *  using Kronecker substitution 
 */
template <auto N>
inline zz_pXi<N> mul(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    zz_pXi<N> c;
    mul(c, a, b);
    return c;
}

template <auto N>
inline zz_pXi<N> operator*(const zz_pXi<N>& a, const zz_pXi<N>& b)
{
    return mul(a, b);
}

/*------------------------------------------------------------*/
/* mul cst                                                    */
/*------------------------------------------------------------*/
template <auto N>
void mul_cst(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_p& b){
    if (b == 0)
    {
        c.zero();
        return;
    }

    c = a;
    long dy = c.rep.length();

    for (long i = 0; i < dy; i++)
        mul_cst(c.rep[i], c.rep[i], b);
    c.normalize();
}

template <>
inline void mul_cst(zz_pXi<1>& c, const zz_pXi<1>& a, const zz_p& b){
    if (b == 0)
    {
        c.zero();
        return;
    }

    c = a;
    long dy = c.rep.length();

    for (long i = 0; i < dy; i++)
        c.rep[i] *= b;
    c.normalize();
}

template <auto N>
inline zz_pXi<N> operator*(const zz_pXi<N>& a, const zz_p& b)
{
    zz_pXi<N> c;
    mul_cst(c, a, b);
    return c;
}

template <auto N>
inline zz_pXi<N> operator*(const zz_p& b, const zz_pXi<N>& a)
{
    zz_pXi<N> c;
    mul_cst(c, a, b);
    return c;
}

/*------------------------------------------------------------*/
/* add cst                                                    */
/*------------------------------------------------------------*/
template <auto N>
void add_cst(zz_pXi<N>& c, const zz_pXi<N>& a, const zz_p& b){
    if (b == 0)
        return;

    c = a;
    if (c.rep.length() == 0) c.rep.SetLength(1);
    add_cst(c.rep[0], c.rep[0], b);
    c.normalize();
}

template <>
inline void add_cst(zz_pXi<1>& c, const zz_pXi<1>& a, const zz_p& b){
    if (b == 0)
        return;

    c = a;
    if (c.rep.length() == 0) c.rep.SetLength(1);
    
    c.rep[0] += b;
    c.normalize();
}

template <auto N>
inline zz_pXi<N> operator+(const zz_pXi<N>& a, const zz_p& b)
{
    zz_pXi<N> c;
    add_cst(c, a, b);
    return c;
}

template <auto N>
inline zz_pXi<N> operator+(const zz_p& b, const zz_pXi<N>& a)
{
    zz_pXi<N> c;
    add_cst(c, a, b);
    return c;
}

template <auto N>
inline zz_pXi<N> operator+=(zz_pXi<N>& a, const zz_p& b)
{
    add_cst(a, a, b);
    return a;
}

#endif
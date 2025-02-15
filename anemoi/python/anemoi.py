from sage.all import *
from hashlib import sha256

Fp = GF(0x64ec6dd0392073)
max_rounds = 21
alpha = 3

bSEED = b"ANEMOI_11_R"
a = Fp(int.from_bytes(sha256(bSEED).digest()[:8], "little"))
ci = [Fp(int.from_bytes(sha256(bSEED + str(r).encode()).digest()[:8], "little")) for r in range(max_rounds)]
di = [Fp(int.from_bytes(sha256(str(r).encode() + bSEED).digest()[:8], "little")) for r in range(max_rounds)]

a_inv = a**(-1)
beta = pow(alpha, -1, Fp.order() - 1)
assert gcd(alpha, Fp.order() - 1) == 1

def H(S):
    x, y = S
    x -= a*y**2 + a_inv
    y -= x**beta
    x += a*y**2

    return x, y

def M(S):
    x, y = S
    return 2*x + y, x + y

def R(S, r):
    return H(M((S[0] + ci[r], S[1] + di[r])))

def anemoi(x, y, r):
    S = x, y
    for ri in range(r):
        S = R(S, ri)
    return M(S)

if __name__ == "__main__":
    print(anemoi(15158938831522922, 0, 11))
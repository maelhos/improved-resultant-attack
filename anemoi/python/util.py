from sage.all import *
from anemoi import *

def print_mat(M):
    v = ord('A')
    x, y = M.dimensions()
    nM = [["0" for _ in range(y)] for __ in range(x)]

    dic = {}
    for i in range(y):
        for j in range(x):
            k = M[i][j]
            if k == 1 or k == 0: nM[i][j] = str(M[i][j])
            else:
                if not k in dic:
                    dic[k] = chr(v)
                    v += 1
                nM[i][j] = dic[k]

        print(str(nM[i]).replace(",", "").replace("'", ""))
    print(dic)

def sym_vertical(M):
    return matrix([k[::-1] for k in M])

## inverse anemoi
def H_inv(S):
    x, y = S
    x -= a*y**2
    y += x.nth_root(alpha)
    x += a*y**2 + a_inv
    return x, y

def M_inv(S):
    x, y = S
    return x - y, 2*y - x

def R_inv(S, r):
    S = M_inv(H_inv(S))
    return S[0] - ci[r], S[1] - di[r]

def anemoi_inv(x, y, r, do_M=True):
    S = x, y
    if do_M: S = M_inv(S)
    for ri in range(r - 1, -1, -1):
        S = R_inv(S, ri)
    return S

if __name__ == "__main__":
    print(anemoi_inv(69, 42, 4, False))
from sage.all import *
from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list
from math import *
from hashlib import sha256
import itertools

p = 0x64ec6dd0392073
d = 3
t = 3 
bSEED = b"RESCUE_PRIME2"
R = 5

F = GF(p)
assert gcd(d, p - 1) == 1
n = p.bit_length() # bit
inv_d = int(1 / Zmod(p-1)(d))


round_constants = [ F(int.from_bytes(sha256(bSEED + str(_).encode()).digest()[:8], "little")) for _ in range(2*R*t)]
MDS_matrix_field = matrix(F, [[2, 1, 1],[1, 2, 1],[1, 1, 2]])

###################################################################

def rescue_prime(input_words, round_constants_counter = 0):

    state_words = list(input_words)

    for r in range(R):
        for i in range(t):
            state_words[i] = (state_words[i])**d
        
        state_words = list(MDS_matrix_field * vector(state_words))

        for i in range(t):
            state_words[i] += round_constants[round_constants_counter]
            round_constants_counter += 1

        for i in range(t):
            state_words[i] = state_words[i]**inv_d
        
        state_words = list(MDS_matrix_field * vector(state_words))

        for i in range(t):
            state_words[i] += round_constants[round_constants_counter]
            round_constants_counter += 1

    return state_words


if __name__ == "__main__":
    C = True

    print("mat =", list(MDS_matrix_field ))
    print("round_constants =", round_constants)
    print(f"{p = }")
    print(f"{d = }")
    print(f"{t = }")
    print(f"{R = }")
        
    #print(MATRIX_FULL)
    test = (23342216854837565, 26176428846319551, 0)
    print("input = ", test)

    print("output = ", rescue_prime(test))




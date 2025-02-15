from math import log2
from sage.all import *
from rescue import *

Fp  = GF(p)

constants_0 = vector(Fp,round_constants[:t]) 
constants_1 = vector(Fp,round_constants[t:2*t])

M = MDS_matrix_field

M_inv = M.inverse()

P1 = PolynomialRing(F, "a, b, c, d")
a, b, c, d = P1.gens()
P2 = PolynomialRing(P1, "X")
X = P2.gen()
round_constants_counter = 0

M_small = matrix(F, [[2, 1], [1, 2]])

csts12 = list(M_small**(-1) * -vector(round_constants[0:2]))
s = [X + csts12[0], -X + csts12[1], 0]

s = list(MDS_matrix_field * vector(s))

for i in range(t):
    s[i] += round_constants[round_constants_counter]
    round_constants_counter += 1

# inverst only the linear factor
s[0] = list(s[0])[1]**inv_d * X
s[1] = list(s[1])[1]**inv_d * X
s[2] **= inv_d

s = list(MDS_matrix_field * vector(s))

for i in range(t):
    s[i] += round_constants[round_constants_counter]
    round_constants_counter += 1

alphas = [0] * 3
betas = [0] * 3
(betas[0], alphas[0]), (betas[1], alphas[1]), (betas[2],) = list(map(list, s))

########### test
x_test = Fp.random_element()
round_constants_counter = 0
state_test = [(x_test + csts12[0])**inv_d, (-x_test + csts12[1])**inv_d, 0]

for i in range(t):
    state_test[i] = (state_test[i])**3

state_test = list(MDS_matrix_field * vector(state_test))

for i in range(t):
    state_test[i] += round_constants[round_constants_counter]
    round_constants_counter += 1

for i in range(t):
    state_test[i] = state_test[i]**inv_d

state_test = list(MDS_matrix_field * vector(state_test))

for i in range(t):
    state_test[i] += round_constants[round_constants_counter]
    round_constants_counter += 1

assert state_test == [alphas[i] * x_test**inv_d + betas[i] for i in range(3)]

####################################
# BACKWARD ROUND ::
x = 1811782577516868
x **= 3

state = vector(F, [(x + csts12[0])**inv_d, (-x + csts12[1])**inv_d, 0])


print("SOLVE VECTOR test :", state)

print("output = ", rescue_prime(state))

#####################################################""
#exit()

print("\nC++ PARAMETERS")

print(f"""// CONSTANTS :
#define RESCUE_P      {hex(int(p))}ULL
#define RESCUE_D      {3}
#define RESCUE_N      {p.bit_length()}
#define RESCUE_T      {t}
#define RESCUE_R      {R - 1}

#define RESCUE_ALPHA0 {alphas[0]}ULL
#define RESCUE_ALPHA1 {alphas[1]}ULL
#define RESCUE_ALPHA2 {alphas[2]}ULL

#define RESCUE_BETA0  {betas[0]}ULL
#define RESCUE_BETA1  {betas[1]}ULL
#define RESCUE_BETA2  {betas[2]}ULL


""")

print(f"static string MDSS(\"{str(list(MDS_matrix_field)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string CSTS(\"{str(list(round_constants)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print("// END CONSTANTS \n")
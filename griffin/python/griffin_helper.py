from hashlib import shake_256
from math import log2,ceil
from sage.all import *

set_verbose(-1)
p = 0x64ec6dd0392073
d = 3
nbr_round_bypass = 3
R = 9
b = 12
seed = b"Griffin"

print(f"Prime p: {p}")
p_bits = p.bit_length()
Fp = GF(p)
assert((p-1) % 3 != 0)
inv_3 = pow(3, -1, p-1)
print(inv_3)

alphas = [0,0] + [4*(i+1) for i in range(b - 2)]
gammas = [0,0] + [i+1 for i in range(b - 2)]
betas  = [0,0] + [(i+1)**2 * 7 for i in range(b - 2)]

set_random_seed(int.from_bytes(seed, "little"))


RC = []
for r in range(R):
    RC.append([])
    for branch in range(b):
        RC[r].append(Fp.random_element())

print(f"Round constants: {RC}")

#Non-linear layer:

def F(i,s,new_s,leading_term=False):
    Li = gammas[i]*new_s[0] + new_s[1] + (s[i-1] if i > 2 else 0)
    if not leading_term:
        return (Li**2 + alphas[i]*Li + betas[i])
    else:
        return Li**2

def non_linear_layer(s,zero=False,leading_term=False):
    new_s = []
    # Supposing that s[0] is only a constant, else the next step is impossible
    # The condition s[0][1] == 0 will be posed later
    if zero:
        new_s = [0,0]
    else:
        new_s.append(s[0]**(inv_3))
        new_s.append(s[1]**3)
    for i in range(2, b):
        new_s.append(F(i,s,new_s,leading_term=leading_term)*s[i])

    return new_s

def non_linear_layer_with_new_variable(s,z):
    new_s = []
    # Supposing that s[0] is only a constant, else the next step is impossible
    # The condition s[0][1] == 0 will be posed later
    new_eq = (s[0][0] - z**3)
    new_s.append(z)
    new_s.append(s[1][0]**3)
    for i in range(2, b):
        new_s.append(F(i,s,new_s)*s[i])

    return new_s,new_eq


def linear_layer(s):
    return list(M*vector(s))

def affine_layer(s, RCi):
    s_2 = linear_layer(s)
    print("POST MUL :", s_2)
    return  [s+c for (s,c) in zip(s_2,RCi)]

# Linear Matrix
M = matrix(Fp, [[10,14,2,6,5,7,1,3,5,7,1,3],
                [8,12,2,2,4,6,1,1,4,6,1,1],
                [2,6,10,14,1,3,5,7,1,3,5,7],
                [2,2,8,12,1,1,4,6,1,1,4,6],
                [5,7,1,3,10,14,2,6,5,7,1,3],
                [4,6,1,1,8,12,2,2,4,6,1,1],
                [1,3,5,7,2,6,10,14,1,3,5,7],
                [1,1,4,6,2,2,8,12,1,1,4,6],
                [5,7,1,3,5,7,1,3,10,14,2,6],
                [4,6,1,1,4,6,1,1,8,12,2,2],
                [1,3,5,7,1,3,5,7,2,6,10,14],
                [1,1,4,6,1,1,4,6,2,2,8,12]])

M_inv = M.inverse()

# Polynomial field representing 
# Ai are linear components after 1st linear layer
# Bi are constant components after 1st linear layer
# z is a new variables needed for 2nd non-linear layer



def compute_states_from_Ai(A2,B2,A4,B4,A6,B6,A8,B8,A11,B11,z):
    x_polynomial_ring = PolynomialRing(constant_ring,'X')
    X, = x_polynomial_ring.gens()


    M_in_x = matrix(x_polynomial_ring,M)
    # Compute the round function for the first few rounds

    # constant components
    state_b = [0,0,B2,0,B4,0,B6,0,B8,0,0,B11] 
    # linear components (multiplied by x)
    state_a = [0,0,A2,0,A4,0,A6,0,A8,0,0,A11] 


    #state after first linear layer
    state_x = [state_a[i]*X + state_b[i] for i in range(12)]


    # state at input (first element is 0 for CICO constaints)
    state_input =  M_inv * vector(state_x)

    #state after first non linear layer
    state_x_after_first_S = non_linear_layer(state_x)


    #state after first affine layer (second application of matrix M)
    state_x_after_first_affine_layer = affine_layer(state_x_after_first_S, RC[0])


    #state after second non linear layer

    state_x_after_second_S,eq_z = non_linear_layer_with_new_variable(state_x_after_first_affine_layer,z)


    #State after second linear layer

    state_x_after_second_affine_layer = affine_layer(state_x_after_second_S, RC[1])
    return state_input,state_x, state_x_after_first_S,state_x_after_first_affine_layer,state_x_after_second_S,state_x_after_second_affine_layer,eq_z




#First step, compute the A2,A4,A6,A8,A10
constant_ring = PolynomialRing(Fp,'A2,A4,A6,A8,A11',order="lex")

A2,A4,A6,A8,A11 = constant_ring.gens()

betas = [0,0] + [(i+1)**2 * 7 for i in range(b - 2)]
state_a = [0,0,A2,0,A4,0,A6,0,A8,0,0,A11] 

state_input = M_inv * vector(state_a)
state_after_first_affine_layer = M*vector([b*s for (b,s) in zip(betas,state_a)])

state_after_second_non_linear_layer = non_linear_layer(state_after_first_affine_layer,zero=True,leading_term=True)



state_after_second_non_linear_layer = vector(state_after_second_non_linear_layer)

state_after_second_affine_layer = M* vector(state_after_second_non_linear_layer)


#state_input,state_x, state_x_after_first_S,state_x_after_first_affine_layer,state_x_after_second_S,state_x_after_second_affine_layer,eq_z = compute_states_from_Ai(A2,0,A4,0,A6,0,A8,0,A10,0,0)


leading_coefficient = []

#First equation: the first input word is 0 :
eq1 = state_input[0]
#Second and third equations, the linear component after first affine layer is 0
eq2 = state_after_first_affine_layer[0]
eq3 = state_after_first_affine_layer[1]

#Fourth and fifth equations, the degree 3 component after second affine layer is 0:

eq4 = state_after_second_affine_layer[0]
eq5 = A11 - 2#state_after_second_affine_layer[1]




eq_system = [eq1,eq2,eq3,eq4,eq5]

#for eq in eq_system:
#    print(eq)
Id = constant_ring.ideal(eq_system)

#id_gb = Id.groebner_basis()
#print(id_gb)
#print(id_gb[-1])
#print(id_gb[-1].factor())
id_var = Id.variety()

#print(len(id_var))
#print(id_var)
a2,a4,a6,a8,a11 = id_var[-1]["A2"],id_var[-1]["A4"],id_var[-1]["A6"],id_var[-1]["A8"],id_var[-1]["A11"]
#print(f"a2, a4, a6, a8, a11 = {a2}, {a4}, {a6}, {a8}, {a11}")
inp = M_inv * vector(Fp, [0, 0, a2, 0, a4, 0, a6, 0, a8, 0, 0, a11])

assert(inp[0] == 0)




# Second step: compute the B2,B4,B6,B8,B10,z

constant_ring = PolynomialRing(Fp,'B2,B4,B6,B8,B11,z',order="lex")


B2,B4,B6,B8,B11,z = constant_ring.gens()

state_input,state_x, state_x_after_first_S,state_x_after_first_affine_layer,state_x_after_second_S,state_x_after_second_affine_layer,eq_z = compute_states_from_Ai(a2,B2,a4,B4,a6,B6,a8,B8,a11,B11,z)


# First equation: the first input word is 0 (constant components = 0)
eq1 = state_input[0][0]

#2nd, 3rd, 4th, 5th equations: the degree 1 and 2 components of 1st and 2nd wire after second affine layer is 0

eq2 = state_x_after_second_affine_layer[0][1]
eq3 = state_x_after_second_affine_layer[0][2]
eq4 = B2 - 1
eq5 = B4 - 1

#print(eq3)
#6th equation: eqz
eq6 = eq_z

eq_system = [eq1,eq2,eq3,eq4,eq5,eq6]

#print([e.degree() for e in state_x_after_second_S])
#print([e.degree() for e in state_x_after_second_affine_layer])

Id = constant_ring.ideal(eq_system)

id_var = Id.variety()

#print(id_var)
#exit()

b2,b4,b6,b8,b11,z = id_var[-1]["B2"], id_var[-1]["B4"],id_var[-1]["B6"],id_var[-1]["B8"],id_var[-1]["B11"], id_var[-1]["z"]


# print(A.right_kernel())
# Vector space of degree 5 and dimension 2 over Finite Field of size 28407454060060787
# Basis matrix:
# [                1                 0 24856522302553193                 0  3550931757507597]
# [                0                 1 27786299816764997 27827710099651383  5452353913374155]

# a2, a4, a6, a8, a10 = A.right_kernel().random_element()
#print(f"A2,B2,A4,B4,A6,B6,A8,B8,A11,B11 = {a2}, {b2}, {a4}, {b4}, {a6}, {b6}, {a8}, {b8}, {a11}, {b11}")
inp_a = M_inv * vector(Fp, [0, 0, a2, 0, a4, 0, a6, 0, a8, 0, 0, a11])

inp_b = M_inv * vector(Fp, [0, 0, b2, 0, b4, 0, b6, 0, b8, 0, 0, b11])


# Given the a_i, compute the state after 3 rounds of Griffin.

x_polynomial_ring = PolynomialRing(Fp,'X')
X, = x_polynomial_ring.gens()

state_inp  = [a*X + b for (a,b) in zip(inp_a,inp_b)]

assert(state_inp[0]==0)

print(f"\nInput State: {state_inp}")

state_after_first_linear_layer = linear_layer(state_inp)


state_after_first_non_linear_layer = non_linear_layer(state_after_first_linear_layer)
print(f"\nLINEAR State: {state_after_first_non_linear_layer}")
state_after_first_affine_layer = affine_layer(state_after_first_non_linear_layer,RC[0])

print(f"\nState after 1st affine layer: {state_after_first_affine_layer}")

state_after_second_non_linear_layer = non_linear_layer(state_after_first_affine_layer)
print(f"\nLINEAR State: {state_after_second_non_linear_layer}")
state_after_second_affine_layer = affine_layer(state_after_second_non_linear_layer,RC[1])

print(f"\nState after 2nd affine layer: {state_after_second_affine_layer}")

state_after_third_non_linear_layer = non_linear_layer(state_after_second_affine_layer)

print(f"\nState after 3rd non-linear layer: {state_after_third_non_linear_layer}")


print("\nC++ PARAMETERS")

print(f"""// CONSTANTS :
#define GRIFFIN_P {hex(int(p))}ULL
#define GRIFFIN_D {d}
#define GRIFFIN_N {p_bits} // NBIT OF p
#define GRIFFIN_R {R}
#define GRIFFIN_B {b}
#define GRIFFIN_BYPASS {nbr_round_bypass}
""")

print(f"static string GRIFFIN_CSTS(\"{str(list(RC)[:-1]).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string GRIFFIN_BYPASS_MUL(\"{str(list(inp_a)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string GRIFFIN_BYPASS_ADD(\"{str(list(inp_b)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string GRIFFIN_ALPHAS(\"{str(list(alphas)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string GRIFFIN_BETAS(\"{str(list(betas)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print(f"static string GRIFFIN_GAMMAS(\"{str(list(gammas)).replace(",", "").replace("(", "[").replace(")", "]")}\");\n")

print("// END CONSTANTS \n")
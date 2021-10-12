from sympy import *
from math import trunc
from random import *
from sympy.physics.quantum import TensorProduct
init_printing(use_latex=True)




##########  MATRICES  ##########

Id16 = eye(16)

X01 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
])

X02 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
])

X03 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
])

X12 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
])

X13 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
])

X23 = Matrix([
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]
])

Id2 = eye(2)

H = 1/sqrt(2)*Matrix([[1,1],[1,-1]])
H0 = TensorProduct(H, Id2, Id2, Id2)
H1 = TensorProduct(Id2, H, Id2, Id2)
H2 = TensorProduct(Id2, Id2, H, Id2)
H3 = TensorProduct(Id2, Id2, Id2, H)

X10 = H0*H1*X01*H1*H0
X20 = H0*H2*X02*H2*H0
X30 = H0*H3*X03*H3*H0
X21 = H2*H1*X12*H1*H2
X31 = H3*H1*X13*H1*H3
X32 = H3*H2*X23*H2*H3


P =  Matrix([[1,0],[0,I]]);
P0 = TensorProduct(P, Id2, Id2, Id2)
P1 = TensorProduct(Id2, P, Id2, Id2)
P2 = TensorProduct(Id2, Id2, P, Id2)
P3 = TensorProduct(Id2, Id2, Id2, P)

X = Matrix([[0,1],[1,0]])
Y = Matrix([[0,-I],[I,0]])
Z = Matrix([[1,0],[0,-1]])

X0 = TensorProduct(X, Id2, Id2, Id2)
X1 = TensorProduct(Id2, X, Id2, Id2)
X2 = TensorProduct(Id2, Id2, X, Id2)
X3 = TensorProduct(Id2, Id2, Id2, X)

Y0 = TensorProduct(Y, Id2, Id2, Id2)
Y1 = TensorProduct(Id2, Y, Id2, Id2)
Y2 = TensorProduct(Id2, Id2, Y, Id2)
Y3 = TensorProduct(Id2, Id2, Id2, Y)

Z0 = TensorProduct(Z, Id2, Id2, Id2)
Z1 = TensorProduct(Id2, Z, Id2, Id2)
Z2 = TensorProduct(Id2, Id2, Z, Id2)
Z3 = TensorProduct(Id2, Id2, Id2, Z)

T =  Matrix([[1,0],[0,exp(I*pi/4)]])
T0 = TensorProduct(T, Id2, Id2, Id2)
T1 = TensorProduct(Id2, T, Id2, Id2)
T2 = TensorProduct(Id2, Id2, T, Id2)
T3 = TensorProduct(Id2, Id2, Id2, T)

Id4 = eye(4)

# transvection matrix [ij]
def transvection(i,j) :
    M = eye(4)
    M[i,j] = 1
    return M

# transposition matrix (ij)
def transpose(i,j) :
    mat = transvection(i,j) * transvection(j,i) * transvection(i,j)
    for i in range(0,4) :
        for j in range(0,4) :
            mat[i,j] = mat[i,j] % 2
    return mat

# rotation around the y axes
def R_y(theta) :
    return cos(theta/2)*Id2 - I*sin(theta/2)*Y

# rotation around the z axes
def R_z(theta) :
    return cos(theta/2)*Id2 - I*sin(theta/2)*Z

# compute a single-qubit unitary operator from a list of three parameters
def single_unitary(param_row) :
    alpha = param_row[0]
    beta = param_row[1]
    gamma = param_row[2]
    return R_z(alpha) * R_y(beta) * R_z(gamma)

# compute a fully factorized unitary operator from a matrix of parameters param and a phase phi 
def factorized_unitary(param, phi) :
    U0 = single_unitary(param.row(0))
    U1 = single_unitary(param.row(1))
    U2 = single_unitary(param.row(2))
    U3 = single_unitary(param.row(3))
    return exp(I*phi)*TensorProduct(U0, U1, U2, U3)

# definition of fac_U that represents any factorized unitary operator ####
phi = symbols('phi')
alpha_0, beta_0, gamma_0 = symbols('alpha_0, beta_0, gamma_0')
alpha_1, beta_1, gamma_1 = symbols('alpha_1, beta_1, gamma_1')
alpha_2, beta_2, gamma_2 = symbols('alpha_2, beta_2, gamma_2')
alpha_3, beta_3, gamma_3 = symbols('alpha_3, beta_3, gamma_3')
param = Matrix([[alpha_0, beta_0, gamma_0],
                [alpha_1, beta_1, gamma_1],
                [alpha_2, beta_2, gamma_2],
                [alpha_3, beta_3, gamma_3]])
fac_U = factorized_unitary(param, phi)

# GL4 matrices of type A^(i,j)_k = [ij][jk][ki][il][lj]
def A_matrix(i,j,k,l) :
    mat = transvection(i,j) * transvection(j,k) * transvection(k,i) * transvection(i,l) * transvection(l,j)
    for i in range(0,4) :
        for j in range(0,4) :
            mat[i,j] = mat[i,j] % 2
    return mat

########################################




########## DICTIONARIES ##########

# dictionary of cnot gates
char_to_cnot = dict(I = Id16, a = X01, b = X02, c = X03, d = X12, e = X13, f = X23, g = X10, h = X20, i = X30, j = X21, k = X31, l = X32)

# dictionary of transvections
char_to_trans = dict(I = Id4, a = transvection(0,1), b = transvection(0,2), c = transvection(0,3), d = transvection(1,2), e = transvection(1,3), f = transvection(2,3), g = transvection(1,0), h = transvection(2,0), i = transvection(3,0), j = transvection(2,1), k = transvection(3,1), l = transvection(3,2))

# dictionary of permutation matrices 
perm = {'Id4' : eye(4),
        '(01)' : transpose(0,1),
        '(02)' : transpose(0,2),
        '(03)' : transpose(0,3),
        '(12)' : transpose(1,2),
        '(13)' : transpose(1,3),
        '(23)' : transpose(2,3),
        '(01)(23)' : transpose(0,1) * transpose(2,3),
        '(02)(13)' : transpose(0,2) * transpose(1,3),
        '(03)(12)' : transpose(0,3) * transpose(1,2),
        '(012)' : transpose(0,1) * transpose(1,2),
        '(021)' : transpose(0,2) * transpose(2,1),
        '(013)' : transpose(0,1) * transpose(1,3),
        '(031)' : transpose(0,3) * transpose(3,1),
        '(123)' : transpose(1,2) * transpose(2,3),
        '(132)' : transpose(1,3) * transpose(3,2),
        '(023)' : transpose(0,2) * transpose(2,3),
        '(032)' : transpose(0,3) * transpose(3,2),
        '(0123)' : transpose(0,1) * transpose(1,2) * transpose(2,3),
        '(0132)' : transpose(0,1) * transpose(1,3) * transpose(3,2), 
        '(0231)' : transpose(0,2) * transpose(2,3) * transpose(3,1),
        '(0213)' : transpose(0,2) * transpose(2,1) * transpose(1,3),
        '(0312)' : transpose(0,3) * transpose(3,1) * transpose(1,2),
        '(0321)' : transpose(0,3) * transpose(3,2) * transpose(2,1)}

# dictionary of cnot gate symbols
char_to_symb = dict(I = 'Id16', a = 'X01', b = 'X02', c = 'X03', d = 'X12', e = 'X13', f = 'X23', g = 'X10', h = 'X20', i = 'X30', j ='X21', k = 'X31', l = 'X32')

# conversion of a string to a unitary matrix
def string_to_unitary(str) :
    unitary = Id16 
    for char in str :
        unitary = unitary * char_to_cnot[char]
    return unitary

# conversion of a string to a GL 4x4 matrix
def string_to_GL_matrix(str) :
    mat = eye(4)
    for char in str :
        mat = mat * char_to_trans[char]
    for i in range(0,4) :
        for j in range(0,4) :
            mat[i,j] = mat[i,j] % 2
    return mat

# conversion of a string to another string representing the cnot product 
def string_to_cnot_symb(str) :
    cnot_string = ''
    for char in str[:-1] :
        cnot_string  += char_to_symb[char]
        cnot_string += '.'
    return cnot_string[:-1]

########################################




########## RIGHT COSETS OF THE SYMMETRIC GROUP IN GL4 ##########

# Each GL4 matrix is encoded as an integer : each line of 4 bits corresponds to an integer L_i between 0 and 15
# and the matrix corresponds to the integer L0 * 2**0 + L1 * 2**4 + L2 * 2**8 + L3 * 2**12
def matrix_to_int(mat) :
    sum = 0
    for i in range(0,4) :
        for j in range(0,4) :
            sum += mat[i,j] * 2**(4*i) * 2**j
    return sum

# Compute the right cosets of the subgroup of GL_4 consisting in the permutation matrices.
# For each GL4 matrix, computes a decomposition in the form of a product of transvection of minimal length.
# This product is encoded as a string following the dictionary char_to_trans =
# dict(I = Id4, a = transvection(0,1), b = transvection(0,2), c = transvection(0,3), d = transvection(1,2), e = transvection(1,3), f = transvection(2,3), g = transvection(1,0), h = transvection(2,0), i = transvection(3,0), j = transvection(2,1), k = transvection(3,1), l = transvection(3,2))
# The string associated to each coset representative is written to the file file_name.
def right_cosets_perm_GL4(file_name) :
    global Id4, char_to_trans, perm
    print("\nComputing the elements of the group GL4 on the field F2 ... (takes about two minutes)")
    GL4 = [('I', Id4, matrix_to_int(Id4))]
    order_of_GL4 = 1
    marks = [0 for i in range(0,2**16)] # marks[i] = 1 iff the matrix corresponding to the integer i is already in GL4
    marks[GL4[0][2]] = 1
    last_elements = [GL4[0]]
    new_elements = []
    while True :
        for elt in last_elements :
            for char in char_to_trans :
                candidate_matrix = char_to_trans[char] * elt[1]
                for i in range(0,4) :
                    for j in range(0,4) :
                        candidate_matrix[i,j] = candidate_matrix[i,j] % 2
                matrix_value = matrix_to_int(candidate_matrix)
                if marks[matrix_value] == 0 : #found new element of GL4
                    new_elt = (char + elt[0], candidate_matrix, matrix_value)
                    GL4.append(new_elt)
                    new_elements.append(new_elt)
                    #print(order_of_GL4)
                    order_of_GL4 += 1
                    marks[matrix_value] = 1
        if len(new_elements) == 0 :
            break
        last_elements = new_elements
        new_elements = []
    print("\nFound %d elements in GL4."%order_of_GL4)
    print("\nComputing the right cosets of the subgroup of permutation matrices in GL4 ...")
    marks = [0 for i in range(0,2**16)] # marks[i] = 1 iff the matrix corresponding to the integer i is already in one coset
    cosets_representatives = []
    order_of_cosets = 0
    for elt in GL4 :
        if marks[elt[2]] == 0 :
            cosets_representatives.append(elt)
            order_of_cosets += 1
            for p in perm :
                matrix_to_mark = perm[p] * elt[1]
                matrix_value = matrix_to_int(matrix_to_mark)
                marks[matrix_value] = 1
    print("\nFound %d cosets."%order_of_cosets)
    f_output = open(file_name, 'w')
    for elt in cosets_representatives :
        f_output.write(elt[0] + '\n')
    f_output.close()

########################################




########## STATES ##########

# definition of the state |L>
phi_plus = 1/sqrt(2)*Matrix([1,0,0,1])
phi_minus = 1/sqrt(2)*Matrix([1,0,0,-1])
psi_plus = 1/sqrt(2)*Matrix([0,1,1,0])
psi_minus = 1/sqrt(2)*Matrix([0,1,-1,0])
u0 = TensorProduct(phi_plus, phi_plus)
u1 = TensorProduct(phi_minus, phi_minus)
u2 = TensorProduct(psi_plus, psi_plus)
L = 1/sqrt(3)*(u0 + exp(I*pi/3)*u1 + exp(-I*pi/3)*u2)

# definition of the state |Phi_5>
Phi5 = 1/sqrt(6) * Matrix([0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,sqrt(2)])

# definition of the state |M_2222> 
v1 = 1/sqrt(6) * Matrix([1,0,0,0,0,1,-1,0,0,-1,1,0,0,0,0,1])
v2 = 1/sqrt(2) * Matrix([0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0])
v3 = 1/sqrt(8) * Matrix([0,-1,1,0,-1,0,0,1,1,0,0,-1,0,1,-1,0])
M2222 = 1/sqrt(8)*v1 + sqrt(6)/4*v2 + 1/sqrt(2)*v3

# definition of the state |0>^4 
zero_state = Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

# definition of a fully factorized four-qubit state
a0, a1, b0, b1, c0, c1, d0, d1 = symbols('a_0 a_1 b_0 b_1 c_0 c_1 d_0 d_1')
a = Matrix([a0, a1])
b = Matrix([b0, b1])
c = Matrix([c0, c1])
d = Matrix([d0, d1])
fac_state = TensorProduct(a, b, c, d)

# definition of the matrices of parameters P_max and P_max' ####
P_max = Matrix([
        [pi/2, pi/2, 0],
        [pi/2, pi/2, 0],
        [pi/4, acos(sqrt(3)/3), 0],
        [pi/4, acos(sqrt(3)/3), 0]])

P_max_prime = Matrix([
        [pi/2, pi/2, 0],
        [pi/2, pi/2, 0],
        [3*pi/4, acos(sqrt(3)/3), 0],
        [3*pi/4, acos(sqrt(3)/3), 0]])

# utilitaries for parameters
def parameters_to_0_2pi(param) :
    for i in range(0,4) :
        for j in range(0,3):
            if param[i,j] >= 2*pi :
                angle = (float(param[i,j]/(2*pi)) - trunc(float(param[i,j]/(2*pi))))*2*pi
            elif param[i,j] <= - 2*pi :
                angle = (param[i,j]/(2*pi) - trunc(float(param[i,j]/(2*pi))))*2*pi + 2*pi
            elif param[i,j] >=0 :
                angle = param[i,j]
            else :
                angle = param[i,j] + 2*pi
            param[i,j]= angle
    return(param)

def parameters_to_minus_pi_pi(param) :
    param=parameters_to_0_2pi(param)
    for i in range(0,4) :
        for j in range(0,3):
            if param[i,j] > pi :
                param[i,j] -= 2*pi
    return(param)
                
# definition of the states |psi_max> and |psi_max'>
# unitary_M01_2 is the CNOT product X01.X12.X20.X03.X31 corresponding to the string 'adhckI'
unitary_M01_2 = string_to_unitary('adhckI')
psi_max = simplify(unitary_M01_2 * factorized_unitary(P_max, 0))
psi_max_prime = simplify(unitary_M01_2 * factorized_unitary(P_max_prime, 0))
    
########################################





########## HYPERDETERMINANT ##########

def hyper_det(state): 
    x0, x1, y0, y1, z0, z1, t0, t1 = symbols('x_0 x_1 y_0 y_1 z_0 z_1 t_0 t_1')
    a0000=state[0]
    a0001=state[1]
    a0010=state[2]
    a0011=state[3]
    a0100=state[4]
    a0101=state[5]
    a0110=state[6]
    a0111=state[7]
    a1000=state[8]
    a1001=state[9]
    a1010=state[10]
    a1011=state[11]
    a1100=state[12]
    a1101=state[13]
    a1110=state[14]
    a1111=state[15]
    A = a0000*x0*y0*z0*t0+    a0001*x0*y0*z0*t1+    a0010*x0*y0*z1*t0+    a0011*x0*y0*z1*t1+    a0100*x0*y1*z0*t0+    a0101*x0*y1*z0*t1+    a0110*x0*y1*z1*t0+    a0111*x0*y1*z1*t1+    a1000*x1*y0*z0*t0+    a1001*x1*y0*z0*t1+    a1010*x1*y0*z1*t0+    a1011*x1*y0*z1*t1+    a1100*x1*y1*z0*t0+    a1101*x1*y1*z0*t1+    a1110*x1*y1*z1*t0+    a1111*x1*y1*z1*t1
    B=a0000*a1111 - a1000*a0111 + a0100*a1011 + a1100*a0011 - a0010*a1101 + a1010*a0101 - a0110*a1001 + a1110*a0001
    L = Matrix([
        [a0000, a0010, a0001, a0011],
        [a1000, a1010, a1001, a1011],
        [a0100, a0110, a0101,a0111],
        [a1100, a1110, a1101, a1111]
    ]).det()
    M = Matrix([
        [a0000, a0001, a0100, a0101],   
        [a1000, a1001, a1100, a1101],
        [a0010, a0011, a0110, a0111],
        [a1010, a1011, a1110, a1111]
    ]).det()
    N = -L-M
    bxy = diff(A,z0,t0)*diff(A,z1,t1)-diff(A,z1,t0)*diff(A,z0,t1)
    Rx = diff(bxy,y0,y0)*diff(bxy,y1,y1)-diff(bxy,y0,y1)*diff(bxy,y1,y0)
    # Rx= C0 * x0**4 + 4*C1 * x0**3 * x1 + 6*C2 * x0**2 * x1**2 + 4*C_3 * x0 * x1**3 + C_4 * x1**4
    C0 = diff(Rx, x0, x0, x0, x0)/24
    C1 = diff(Rx, x0, x0, x0, x1)/24
    C2 = diff(Rx, x0, x0, x1, x1)/24
    C3 = diff(Rx, x0, x1, x1, x1)/24
    C4 = diff(Rx, x1, x1, x1, x1)/24 
    S = C0*C4 - 4*C1*C3 + 3*C2**2
    T = C0*C2*C4 - C0*C3**2 + 2*C1*C2*C3 - C1**2*C4 - C2**3
    return S**3 - 27*T**2

# return True iff the hyperdeterminant is null for a state resulting from a cnot circuit (corresponding to the string str) applied on any factorized state.
def is_HD_null(str) :
    unitary = string_to_unitary(str)
    image = unitary * fac_state
    HD = hyper_det(image)
    return expand(HD) == 0

# write to a file the strings str such that is_HD_null(str) is false.
def non_zero_HD_strings(input_file_name, output_file_name) :
    cosets_count = 0
    non_zero_HD_count = 0
    f_input = open(input_file_name, 'r')
    f_output = open(output_file_name, 'w')
    str = f_input.readline()
    while str != '' :
        cosets_count += 1
        print("\nComputing the hyperdeterminant for coset %d ..."%cosets_count)
        if is_HD_null(str[:-1]) == False :
            print(str[:-1])
            f_output.write(str)
            non_zero_HD_count += 1
            print("\nNon zero HD_count = %d\n"%non_zero_HD_count)
        else :
            print("\nHD = 0\n")
        str = f_input.readline()
    f_input.close()
    f_output.close()

########################################




########## OPTIMIZATION PROGRAMS ##########
# In all functions, we use a random walk heurisitic, that is a type of random search on the space of parameters.


# Search for a matrix of parameters param (where the last column is null) such that the state resulting from the cnot circuit (defined by the string str) applied on the factorized state (defined by param) maximizes the absolute value of the hyperdeterminant HD.
# In the first loop, the program tries to find parameters such that the hyperdeterminant is greater than 8.8e-8, because there are some local maxima around 8.7e-8.
# If such parameters have been found, the program tries (in the second loop) to optimize the matrix param to make |HD| maximal (around 1.98e-7)
def search_max_HD(str) :
    global fac_state
    cnot_prod = string_to_unitary(str);
    image = cnot_prod * fac_state
    HD = hyper_det(image)
    max_try = 4
    nb_try = 0
    shots = 50
    ampli_limit = 1.e-3
    factor = 0.5
    max_HD = 8.8e-8
    found_param = False
    # first loop : try to make |HD| greater than 8.8e-8
    while found_param == False and nb_try <= max_try : 
        amplitude = 1
        nb_try += 1
        print('\n********** Attempt %d/%d to make |HD| greater than %e  **********'%(nb_try, max_try, max_HD))
        param = Matrix(4, 3, lambda i,j : (2*random()-1)*pi)
        for i in range(0,4) :
            param[i,2] = 0.
        U0 = single_unitary(param.row(0))
        U1 = single_unitary(param.row(1))
        U2 = single_unitary(param.row(2))
        U3 = single_unitary(param.row(3))
        HD_eval = (HD.subs([(a0, U0[0,0]), (a1, U0[1,0]), (b0, U1[0,0]), (b1, U1[1,0]), (c0, U2[0,0]), (c1, U2[1,0]), (d0, U3[0,0]), (d1, U3[1,0])])).evalf()
        current_max = abs(HD_eval)
        print('\ninitial max= %.15e'%current_max)
        print("\ninitial parameters = \n")
        pprint(param)
        while current_max < max_HD and amplitude > ampli_limit : 
            found_better = False
            print("\n\nStarting a new loop of %d shots to increase |HD|."%shots)
            print('\namplitude = %e (limit = %e)\n'%(amplitude, ampli_limit))
            for i in range(0, shots) :
                print('shot %d/%d\n'%(i+1, shots))
                direction = Matrix(4, 3, lambda i,j : (2*random()-1)*pi)
                for j in range(0,4) :
                    direction[j,2] = 0
                candidate = param + amplitude * direction
                U0 = single_unitary(candidate.row(0))
                U1 = single_unitary(candidate.row(1))
                U2 = single_unitary(candidate.row(2))
                U3 = single_unitary(candidate.row(3))
                HD_eval = (HD.subs([(a0, U0[0,0]), (a1, U0[1,0]), (b0, U1[0,0]), (b1, U1[1,0]), (c0, U2[0,0]), (c1, U2[1,0]), (d0, U3[0,0]), (d1, U3[1,0])])).evalf()
                if abs(HD_eval) > current_max :
                    current_max = abs(HD_eval)
                    param = Matrix(candidate)
                    found_better = True
                    print('Break ! Found new maximum = %.15e'%current_max)
                    print("\nparameters = \n")
                    pprint(param)
                    break
            if found_better == False :
                amplitude = amplitude*factor
        if current_max >= max_HD :
            found_param = True
            print("\n**********  Attempt %d successful to make |HD| greater than %e.  **********"%(nb_try,max_HD))
        else :
            print("\n************  Search failed, let's go for another attempt !  ***************\n")
    if found_param == False :
        print("******************* Search failed after %d attempts : end of program  *******************\n"%max_try)
        return
    # second loop : try to make |HD| equal to the theoretical maximum 
    shots = 25
    ampli_limit = 1.e-13
    accuracy = 1.0e-22
    max_HD = 1/2**8/3**9
    print("\n**********  Now we try to make |HD| equal to %e  **********\n"%max_HD)
    while max_HD - current_max > accuracy and amplitude > ampli_limit :
        print("\n\nStarting a new loop of %d shots to increase |HD|."%shots)
        print('\namplitude = %e (limit = %e)\n'%(amplitude, ampli_limit))
        found_better = False
        for i in range(0, shots) :
            print('shot %d/%d\n'%(i+1, shots))
            direction = Matrix(4, 3, lambda i,j : (2*random()-1)*pi)
            for j in range(0,4) :
                direction[j,2] = 0
            candidate = param + amplitude * direction
            U0 = single_unitary(candidate.row(0))
            U1 = single_unitary(candidate.row(1))
            U2 = single_unitary(candidate.row(2))
            U3 = single_unitary(candidate.row(3))
            HD_eval = (HD.subs([(a0, U0[0,0]), (a1, U0[1,0]), (b0, U1[0,0]), (b1, U1[1,0]), (c0, U2[0,0]), (c1, U2[1,0]), (d0, U3[0,0]), (d1, U3[1,0])])).evalf()
            if abs(HD_eval) > current_max :
                current_max = abs(HD_eval)
                param = Matrix(candidate)
                found_better = True
                print('\nBreak ! Found new maximum = %.15e'%current_max)
                print("\nparameters = \n")
                pprint(param)
                break
        if found_better == False :
            amplitude = amplitude * factor
    if max_HD - current_max <= accuracy :
        print("\n**********  Search successful *********\n")
        param = parameters_to_minus_pi_pi(param)
        print('|HD| = %.15e'%current_max)
        print("\nparameters = \n")
        pprint(param)
    else :
        print("*********  Search failed : end of program  *********\n")

# search for a matrix of parameters param and a phase phase such that :
# |state2> = exp(i*phase) * U(param) * |state1>
# try to make the difference |state2> - exp(i*phase) * U(param) * |state1> equal to the null vector
def search_LU_from_state1_to_state2(state1, state2) :
    global fac_U
    state_to_minimize = state2 - fac_U * state1
    # the random walk tries to minimize up to null vector the state state_to_minimize 
    shots = 100
    factor = 0.5
    ampli_limit = 1.e-8
    zero_limit = 1.e-5
    amplitude = 1.
    phase = (2*random()-1)*pi
    param = Matrix(4, 3, lambda i,j : (2*random()-1)*pi)
    state = (state_to_minimize.subs([(alpha_0, param[0,0]), (beta_0, param[0,1]), (gamma_0,param[0,2]),
                                     (alpha_1, param[1,0]), (beta_1, param[1,1]), (gamma_1,param[1,2]),
                                     (alpha_2, param[2,0]), (beta_2, param[2,1]), (gamma_2,param[2,2]),
                                     (alpha_3, param[3,0]), (beta_3, param[3,1]), (gamma_3,param[3,2]), (phi, phase)])).evalf()
    # state_to_minimize is null iff the sum of the absolute value of its coordinates is zero
    sum_coord = 0
    for i in range(0,16) :
        sum_coord += abs(state[i])
    min = sum_coord
    print('\ninitial minimum = %.15e'%min)
    print('\ninitial parameters = \n')
    pprint(param)
    while min > zero_limit and amplitude > ampli_limit :
        print("\n\nStarting a new loop of %d shots to decrease the minimum."%shots)
        print('\namplitude = %e (limit = %e)\n'%(amplitude, ampli_limit))
        found_better = False
        for i in range(0, shots) :
            print('shot %d/%d\n'%(i+1, shots))
            direction = Matrix(4, 3, lambda i,j : (2*random()-1)*pi)
            param_candidate = param + amplitude * direction
            phase_candidate = phase + amplitude * (2*random()-1)*pi
            state = (state_to_minimize.subs([(alpha_0, param_candidate[0,0]), (beta_0, param_candidate[0,1]), (gamma_0,param_candidate[0,2]),
                                             (alpha_1, param_candidate[1,0]), (beta_1, param_candidate[1,1]), (gamma_1,param_candidate[1,2]),
                                             (alpha_2, param_candidate[2,0]), (beta_2, param_candidate[2,1]), (gamma_2,param_candidate[2,2]),
                                             (alpha_3, param_candidate[3,0]), (beta_3, param_candidate[3,1]), (gamma_3,param_candidate[3,2]), (phi, phase_candidate)])).evalf()
            sum_coord = 0
            for i in range(0,16) :
                sum_coord += abs(state[i])
            if sum_coord < min :
                min = sum_coord
                param = Matrix(param_candidate)
                phase = phase_candidate
                found_better = True
                print('Break ! Found new minimum = %.15e\n'%min)
                print('parameters = \n')
                pprint(param)
                print('\nphase = %e\n'%phase)
                break
        if found_better == False :
            amplitude = amplitude * factor
    if min < limit_zero :
        print("\n************  search successful !  ***************\n")
    else :
        print("\n**********************  search failed  ************************* \n")
    print('minimum found = %e\n'%min)
    print('parameters = \n')
    pprint(param)
    print('\nphase = %e\n'%phase)
    
########################################




####### VERIFICATIONS ##########

def check_psi_max_is_MHS() :
    global psi_max, psi_max_prime
    # computation of the states defined by Identities (24) and (25)
    w1 = 1/sqrt(8) * Matrix([0,1,0,I,0,1,0,-I,1,0,I,0,1,0,-I,0])
    w2 = Rational(1,2) * Matrix([-1,0,0,0,0,0,-I,0,0,0,0,-I,0,1,0,0])
    w3 = Rational(1,2) * Matrix([0,0,1,0,I,0,0,0,0,-I,0,0,0,0,0,1])
    state_eq24 = sqrt(3)/3 * w1 + (3 + sqrt(3))/6 * exp(I*pi/4) * w2 + (3 - sqrt(3))/6 * exp(I*pi/4) * w3
    state_eq25 = sqrt(3)/3 * w1 + (3 + sqrt(3))/6 * exp(-I*pi/4) * w2 + (3 - sqrt(3))/6 * exp(I*3*pi/4) * w3
    print("\n******** Checking Identity (24) ********\n")
    print("\nDefinition of the state |psi_max> from Identity (21) :\n")
    print("|psi_max> = ...\n")
    input("Press ENTER to continue.")
    pprint(psi_max)
    print("\nThe form of |psi_max> resulting from Identity (24) is : \n")
    input("Press ENTER to continue.\n")
    pprint(state_eq24)
    print("\n\nTo check Identity (24), one computes the difference between the two vectors above and find the null vector :\n")
    input("Press ENTER to continue.\n")
    pprint(simplify(re(psi_max - state_eq24))+ I* simplify(im(psi_max - state_eq24)))
    print("\n******** Checking Identity (25) ********\n")
    print("\nDefinition of the state |psi_max'> from Identity (22) :\n")
    print("|psi_max'> = ...\n")
    input("Press ENTER to continue.\n")
    pprint(psi_max_prime)
    print("\nThe form of |psi_max'> resulting from Identity (25) is : \n")
    input("Press ENTER to continue.\n")
    pprint(state_eq25)
    print("\n\nTo check Identity (25), one computes the difference between the two vectors above and find the null vector :\n")
    input("Press ENTER to continue.\n")
    pprint(simplify(re(psi_max_prime - state_eq25))+I*simplify(im(psi_max_prime - state_eq25)))
    print("\n******** Checking Identity (23) ********\n")
    print("\nComputing the value of the hyperdeterminant for the state |psi_max> (takes about five minutes...)")
    pprint(simplify(hyper_det(psi_max)))
    print("\nComputing the value of the hyperdeterminant for the state |psi_max'> (takes about five minutes...)")
    pprint(simplify(hyper_det(psi_max_prime)))

def check_psi_to_psi_prime() :
    global psi_max, psi_max_prime
    print("\n******** This function checks Equality (29) in Remark 2. ********\n")
    print("This equality gives the LU operators that generate |psi_max'> from |psi_max>.\n")
    print("The equality is true iff the difference between its two members is the null vector.\n")
    input("Press ENTER to continue.\n")
    print("Computing the difference ... \n")
    right_state = exp(-I*7*pi/12) * simplify(TensorProduct(P*H*P**(-1), P*X, H, P*H*P**(-1)) * psi_max)
    pprint(simplify(re(psi_max_prime - right_state)) + I * simplify(im(psi_max_prime - right_state)))

def check_distinct_cosets() :
    print("\n******** Proof of proposition 9 *******")
    print("\nConstruction of the list of all matrices A^(i,j)_k (24 matrices) : ")
    A_list = []
    for i in range(0,4) :
        for j in range(0,4) :
            if j != i :
                for k in range(0,4) :
                    if k !=i and k !=j :
                        for l in range(0,4) :
                            if l !=i and l !=j and l !=k :
                                A_list.append(((i,j,k),A_matrix(i,j,k,l)))
    input("\nPress ENTER to continue.")
    for (tuple,mat) in A_list :
        print("A^(%d,%d)_%d="%(tuple[0],tuple[1],tuple[2]))
        pprint(mat)
        print('\n')
    print("\nIn the above list, we search for the matrices which are a line permutation of another matrix of the list :")
    input("\nPress ENTER to continue.")
    for p in range(0,len(A_list)) :
        for q in perm :
            mat = perm[q]*A_list[p][1]
            for r in range(p+1,len(A_list)) :
                if mat == A_list[r][1] :
                    print('\n\n')
                    print("A^(%d,%d)_%d = %s A^(%d,%d)_%d"%(A_list[r][0][0],A_list[r][0][1],A_list[r][0][2],q,A_list[p][0][0],A_list[p][0][1],A_list[p][0][2]))
    print("\n\nWe see that A^(i,j)_k is a line permutation of A^(p,q)_r iff (i,j)=(p,q).\n")

def check_LU_operators() :
    print("\n******** This function  checks the three equalities in Figure 2 ********\n")
    global psi_max, L, Phi5, M2222
    print("\n******** Checking the first equality in Figure 2 ********\n")
    print("This equality is true iff the difference between the state |L> and the state exp(-I*11pi/12) * U(P_psi_to_L) * |psi_max> is the null vector.\n")
    input("Press ENTER to continue.")
    print("Computing this difference ... \n")
    param_psi_to_L = Matrix([[pi, pi/2, -pi/2],
                             [pi/2, -pi/2, pi],
                             [0, pi, pi],
                             [0, -pi/2, pi/2]])
    phase_psi_to_L = -11*pi/12
    check_L = L - factorized_unitary(param_psi_to_L, phase_psi_to_L) * psi_max
    pprint(simplify(re(check_L))+ I * simplify(im(check_L)))
    print("\n******** Checking the second equality in Figure 2 ********\n")
    print("This equality is true iff the difference between the state |Phi5> and the state exp(-I*7pi/12) * U(P_psi_to_Phi5) * |psi_max> is the null vector.\n")
    input("Press ENTER to continue.")
    print("Computing this difference ... (takes about one minute...)\n")
    theta = acos(sqrt(3)/3)
    param_psi_to_Phi5 = Matrix([[-pi/3, theta - pi, -3*pi/4],
                             [pi/3, theta, 3*pi/4],
                             [pi, theta, 3*pi/4],
                             [2*pi/3, pi - theta, pi/4]])
    phase_psi_to_Phi5 = -7*pi/12
    check_Phi5 = Phi5 - factorized_unitary(param_psi_to_Phi5, phase_psi_to_Phi5) * psi_max
    pprint(expand(simplify(re(check_Phi5)))+ I * expand(simplify(im(check_Phi5))))
    print("\n******** Checking the third equality in Figure 2 ********\n")
    print("This equality is true iff the difference between the state |M2222> and the state exp(I*5pi/12) * U(P_psi_to_M2222) * |psi_max> is the null vector.\n")
    input("Press ENTER to continue.")
    print("Computing this difference ... (takes about 30 sec)\n")
    param_psi_to_M2222 = Matrix([[pi/2, pi/4, 0],
                             [-pi/2, -pi/4, pi/2],
                             [0, -pi/2, pi/4],
                             [pi/2, 3*pi/4, pi]])
    phase_psi_to_M2222 = 5*pi/12
    check_M2222 = M2222 - factorized_unitary(param_psi_to_M2222, phase_psi_to_M2222) * psi_max
    pprint(simplify(re(check_M2222))+ I * simplify(im(check_M2222)))
    
def check_circuits() :
    global cnot_prod_M01_2, zero_state, L, Phi5, M2222
    print("\nThis function checks that the three circuits in Figure 3 generate the states |L>, |Phi_5> and |M_2222>.\n")
    theta = acos(sqrt(3)/3)
    print("\nChecking first circuit : the difference between the state generated by this circuit and  exp(-i*7pi/12)|L> must be null.\n")
    input("Press ENTER to continue.")
    print("\nComputing this difference ... \n")
    layer1 = TensorProduct(H, H, R_y(theta), R_y(theta))
    layer2 = TensorProduct(P, P, T, T)
    layer3 = cnot_prod_M01_2
    layer4 = TensorProduct(P, Z, X, P)
    layer5 = TensorProduct(H, H, Id2, H)
    layer6 = TensorProduct(Z, P**(-1), Id2, Z)
    state_circ1 = simplify(layer6 * layer5 * layer4 * layer3 * layer2 * layer1 * zero_state)
    pprint(expand(re(state_circ1 - exp(-I*7*pi/12)*L)) + I * expand(im(state_circ1 - exp(-I*7*pi/12)*L)))
    print("\nChecking second circuit : the difference between the state generated by this circuit and  exp(-i*pi/6)|Phi5> must be null.\n")
    input("Press ENTER to continue.")
    print("\nComputing this difference ... \n")
    layer4 = TensorProduct(T, T, T, T)
    layer5 = TensorProduct(X, P, P, Y)
    layer6 = TensorProduct(R_y(theta), R_y(theta), R_y(theta), R_y(-theta))
    layer7 = TensorProduct(R_z(-pi/3), R_z(pi/3), Z, R_z(2*pi/3))
    state_circ2 = simplify(layer7 * layer6 * layer5 * layer4 * layer3 * layer2 * layer1 * zero_state)
    pprint(expand(re(state_circ2 - exp(-I*pi/6)*Phi5)) + I * expand(im(state_circ2 - exp(-I*pi/6)*Phi5)))
    print("\nChecking third circuit : the difference between the state generated by this circuit and  exp(-i*13pi/24)|M2222> must be null.\n")
    input("Press ENTER to continue.")
    print("\nComputing this difference ... (takes about one minute) \n")
    layer4 = TensorProduct(Id2, P, T, Z)
    layer5 = TensorProduct(R_y(pi/4), R_y(-pi/4), H, R_y(3*pi/4))
    layer6 = TensorProduct(P, P**(-1), Z, P)
    state_circ3 = simplify(layer6 * layer5 * layer4 * layer3 * layer2 * layer1 * zero_state)
    # r is an ad hoc coefficient to make simplification easier
    r = 96*sqrt(2+sqrt(2))
    pprint(simplify(expand(r*factor(re(state_circ3 - exp(-I*13*pi/24)*M2222)))) + I * simplify(expand(r*factor((im(state_circ3 - exp(-I*13*pi/24)*M2222))))))

########################################

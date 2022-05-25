from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.special import block_matrix

### Input and output transformations for the linear equivalence
M = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])
N = Matrix(GF(2), [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]])
input_transformation = block_matrix([[M if i == j else N if (i+1)%32 == j else 0 for i in range(32)] for j in range(32)])
output_transformation = input_transformation.inverse()

### Build the alternative linear layer
A = Matrix(GF(2),[[0, 1, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
B = Matrix(GF(2),[[1, 1, 0, 1], [1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
C = Matrix(GF(2),[[0, 0, 0, 0], [1, 1, 0, 1], [1, 1, 1, 1], [0, 0, 0, 0]])
D = Matrix(GF(2),[[1, 1, 0, 0], [0, 0, 0, 0], [1, 1, 0, 1], [0, 0, 0, 0]])
E = Matrix(GF(2),[[1, 1, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
F = Matrix(GF(2),[[1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 0, 0]])
G = Matrix(GF(2),[[1, 1, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 0, 1]])
H = Matrix(GF(2),[[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 0, 1]])
I = Matrix(GF(2),[[1, 1, 0, 1], [1, 0, 0, 1], [0, 0, 0, 0], [1, 1, 0, 1]])
J = Matrix(GF(2),[[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
K = Matrix(GF(2),[[0, 0, 0, 0], [1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
L = Matrix(GF(2),[[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [0, 0, 0, 0]])
M = Matrix(GF(2),[[0, 0, 0, 0], [1, 1, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
N = Matrix(GF(2),[[1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
O = Matrix(GF(2),[[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 0, 1], [0, 0, 0, 0]])

alternative_linear_layer = block_matrix([
	[A, B, C, D, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0],
	[0, 0, 0, F, H, B, C, D, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, F, H, B, C, D, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, H, B, C, D, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, H, B, C, D, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, H, B, C, D, E, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, H, B, C, D, E, 0, 0, 0],
	[E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, H, B, C, D],
	[I, C, D, J, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F],
	[0, 0, F, G, K, C, D, J, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, F, G, K, C, D, J, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, K, C, D, J, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, K, C, D, J, E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, K, C, D, J, E, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, K, C, D, J, E, 0, 0, 0],
	[E, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, K, C, D, J],
	[L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G],
	[0, F, G, 0, L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, F, G, 0, L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, L, D, J, B, M, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, L, D, J, B, M, 0, 0, 0],
	[M, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, L, D, J, B],
	[N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0],
	[F, G, 0, 0, N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, F, G, 0, 0, N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0, N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0, N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0, N, J, B, C, O, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0, N, J, B, C, O, 0, 0, 0],
	[O, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, F, G, 0, 0, N, J, B, C]
])

### Build GIFT-128 permutation-matrix
gift_permutation = [0  ,33 ,66 ,99 ,96 ,1  ,34 ,67 ,64 ,97 ,2  ,35 ,32 ,65 ,98 ,3  ,
					4  ,37 ,70 ,103,100,5  ,38 ,71 ,68 ,101,6  ,39 ,36 ,69 ,102,7  ,
					8  ,41 ,74 ,107,104,9  ,42 ,75 ,72 ,105,10 ,43 ,40 ,73 ,106,11 ,
					12 ,45 ,78 ,111,108,13 ,46 ,79 ,76 ,109,14 ,47 ,44 ,77 ,110,15 ,
					16 ,49 ,82 ,115,112,17 ,50 ,83 ,80 ,113,18 ,51 ,48 ,81 ,114,19 ,
					20 ,53 ,86 ,119,116,21 ,54 ,87 ,84 ,117,22 ,55 ,52 ,85 ,118,23 ,
					24 ,57 ,90 ,123,120,25 ,58 ,91 ,88 ,121,26 ,59 ,56 ,89 ,122,27 ,
					28 ,61 ,94 ,127,124,29 ,62 ,95, 92 ,125,30 ,63 ,60 ,93 ,126,31 ]
permutation = Matrix(GF(2), Permutation([1 + i for i in gift_permutation]).to_matrix())

# Setup round constants as 128 bit vectors
round_const = [vector(GF(2), 128, ZZ(
	1<<127 | ((c & 32) << 18) | ((c & 16) << 15) | ((c & 8) << 12) | ((c & 4) << 9) | ((c & 2) << 6) | ((c & 1) << 3)
).digits(2, padto=128)) for c in [1, 3, 7, 15, 31, 62, 61, 59, 55, 47, 30, 60, 57, 51, 39, 14, 29, 58, 53, 43, 22, 44, 24, 48, 33, 2, 5, 11]]

def substitution_layer_default(x):  # Substitution layer for DEFAULT Layer
	s = copy(x)
	sum0 = x[0::4] + x[3::4]
	sum1 = x[1::4] + x[2::4]
	prod = sum0.pairwise_product(sum1)
	s[0::4] = x[0::4] + x[1::4] + x[2::4]
	s[1::4] = prod + x[0::4] + x[1::4]
	s[2::4] = x[1::4] + x[2::4] + x[3::4]
	s[3::4] = prod + x[2::4] + x[3::4]
	return s
	
def substitution_layer_core(x):  # Substitution layer for DEFAULT Core
	s = copy(x)
	one = vector(GF(2), 32, ZZ(0xffffffff).digits(2, padto=32))
	s[0::4] = x[0::4].pairwise_product(x[1::4] + x[2::4]) + x[1::4] + x[3::4] + one
	s[1::4] = x[0::4].pairwise_product(x[2::4]) + x[1::4] + x[2::4] + x[3::4]
	s[2::4] = x[0::4].pairwise_product(x[3::4]) + x[1::4] + x[2::4]
	s[3::4] = x[0::4].pairwise_product(x[3::4]).pairwise_product(x[1::4] + one) + x[2::4].pairwise_product(x[1::4] + x[3::4]) + x[0::4] + x[3::4]
	return s

def round_default(x, k):  # Round function for DEFAULT Layer
	x = substitution_layer_default(x)
	x = permutation * x
	x += k
	return x
	
def round_default_linear_equivalent(x, k):  # Round function for DEFAULT Layer of the linear equivalent version
	x = substitution_layer_default(x)
	x = alternative_linear_layer * x
	x += input_transformation * k
	return x
	
def round_core(x, k):  # Round function for DEFAULT Core
	x = substitution_layer_core(x)
	x = permutation * x
	x += k
	return x

def key_schedule(key, rounds):  # Key schedule (for both versions)
	keys = [copy(key)]
	c = vector(GF(2), 128, ZZ(1<<127).digits(2, padto=128))
	for i in range(3):
		keys.append(round_default(round_default(round_default(round_default(keys[i], c), c), c), c))
	ks = []
	for i in range(rounds):
		ks.append(keys[i % 4] + round_const[i])
	return ks
	
def encrypt(x, key):  # Encrypt using the standard representation
	ks = key_schedule(key, 28)
	for k in ks:
		x = round_default(x, k)
	ks = key_schedule(key, 24)
	for k in ks:
		x = round_core(x, k)
	ks = key_schedule(key, 28)
	for k in ks:
		x = round_default(x, k)
	return x
	
def encrypt_linear_equivalent(x, key):  # Encrypt using the linear equivalent version of the DEFAULT Layer
	# Prepare input according to linear equivalence
	x = input_transformation * x
	ks = key_schedule(key, 28)
	for k in ks:
		x = round_default_linear_equivalent(x, k)
	# Prepare output according to linear equivalence
	x = output_transformation * x
	ks = key_schedule(key, 24)
	for k in ks:
		x = round_core(x, k)
	# Prepare input according to linear equivalence
	x = input_transformation * x
	ks = key_schedule(key, 28)
	for k in ks:
		x = round_default_linear_equivalent(x, k)
	# Prepare output according to linear equivalence
	x = output_transformation * x
	return x
	
def check_test_vectors():
	# Check the test vectors provided by the designers for both versions
	for enc_func in [encrypt, encrypt_linear_equivalent]:
		pts = [0x00000000000000000000000000000000, 0x33333333333333333333333333333333, 0x55555555555555555555555555555555, 0xe1e51e2e08f8588d6fb85911b25a1829]
		mks = [0x00000000000000000000000000000000, 0x33333333333333333333333333333333, 0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa, 0x974c0adaa33900495909bea963df0a19]
		cts = [0x93faff138c527a052e5c996278280244, 0x68902d38bed0d8a19c420cfc3c0d3d9a, 0xb601610542b82ae8432c1117875b16be, 0xf9194b9928ff08c768398afaa59bd0f3]
		for i in range(len(pts)):
			pt = vector(GF(2), ZZ(pts[i]).digits(2, padto=128))
			mk = vector(GF(2), ZZ(mks[i]).digits(2, padto=128))
			ct = enc_func(pt, mk)
			print(ZZ(list(ct), 2).hex(), ct == vector(GF(2), ZZ(cts[i]).digits(2, padto=128)))
	
check_test_vectors()

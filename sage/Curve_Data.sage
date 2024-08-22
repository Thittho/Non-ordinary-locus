# Function to compute Cartier matrix
def CartierMatrix(f, p, r):
    Fp = GF(p**r)
    g = (f.degree()-1) / 2
    h = f**((p-1) / 2)
    M = matrix(Fp, g, g)

    for i in range(1, g+1):
        for j in range(1, g+1):
            if ((p*i) - j) >= len(h.list()):
                M[i-1,j-1] = 0
            else:
                M[i-1,j-1] = h.list()[(p*i)-j]
   
    return M

# Function for checking if polynomial is separable
def IsSeparable(f):
    f_prime = f.derivative()
    return f.gcd(f_prime) == 1

# Function to compute data for all curves over GF(p^r)
def CurveDataFull(p, r):
    Fp = GF(p**r)
    R.<x> = PolynomialRing(Fp)

    Poly_Counter = 0
    p_ranks = []

    for Tuple in Tuples(Fp, 6):
        f = R(x**9 + x**7 + x**6 + sum(Tuple[i]*x**i for i in range(6)))
        if IsSeparable(f):
            p_ranks.append((CartierMatrix(f, p, r)**4).rank())
            Poly_Counter += 1

    res = [0, 0, 0, 0, 0]
    for i in p_ranks:
        res[i] += 1
    
    return p, r, Poly_Counter, res[0], res[1], res[2], res[3], res[4]

# Function to compute data for sample of curves over GF(p^r)
def CurveDataSample(p, r, s, *args):
    Fp = GF(p**r)
    R.<x> = PolynomialRing(Fp)

    Poly_Counter = 0
    p_ranks = []

    if len(args) == 0:
        Tuple_List = [ [Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element()] for i in range(5*s / 4)]
    elif len(args) == 1:
        Tuple_List = [ [args[0], Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element()] for i in range(5*s / 4)]
    elif len(args) == 2:
        Tuple_List = [ [args[0], args[1], Fp.random_element(), Fp.random_element(), Fp.random_element(), Fp.random_element()] for i in range(5*s / 4)]

    else:
        return print("Too many arguments specified!")

    L = list(Set(Tuple_List))

    for Tuple in L:
        if Poly_Counter < s:
            f = R(x**9 + x**7 + x**6 + sum(Tuple[i]*x**i for i in range(6)))
            if IsSeparable(f):
                p_ranks.append((CartierMatrix(f, p, r)**4).rank())
                Poly_Counter += 1

    if len(p_ranks) < s:
        print("NOTE: Only picked up", len(p_ranks), "curves instead of", s)

    res = [0, 0, 0, 0, 0]
    for i in p_ranks:
        res[i] += 1
    
    return p, r, Poly_Counter, res[0], res[1], res[2], res[3], res[4]
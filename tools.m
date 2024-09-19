intrinsic NotRoot(Fp::FldFin, s::RngInt) -> FldFinElt
    {Given a field Fp and an integer s, returns an element of Fp which is not a s-th power. Return an error if there is no such element.}
    Fpx<x> := PolynomialRing(Fp);
    if #Roots(x^s-1) eq 1 then
        return "Error: every element of Fp is a s-th power.";
    end if;
    a := Random(Fp);
    while IsPower(a, s) do
        a := Random(Fp);
    end while;
    return a;
end intrinsic;

intrinsic CartierMatrix(f::RngUPolElt, p::RngIntElt : r := 1) -> .
    {Given a univariate polynomial f with coefficients in F_p^r (or coercible in F_p^r), return its Cartier-Manin matrix.}
    
    if Type(BaseRing(Parent(f))) ne FldFin then
        assert IsCoercible(PolynomialRing(GF(p^r)), f);
        f := PolynomialRing(GF(p^r))!f;
    end if;

    g := (Degree(f)-1) div 2;
    h := f^((p-1) div 2);
    return Matrix([[p*i-j ge 0 select Coefficient(h, p*i-j) else 0 : i in [1..g]] : j in [1..g]]);
end intrinsic;

intrinsic StableCartierManin(M::Any, p::RngIntElt : r := 1) -> .
    {}
    if r eq 1 then 
        return M^(Ncols(M));
    end if;
    return &*([M] cat [Matrix([[M[i,j]^(p^k) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]) : k in [1..Ncols(M)-1]]);
end intrinsic;

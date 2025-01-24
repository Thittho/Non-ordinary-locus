intrinsic HasseWittMatrix(f::RngUPolElt, p::RngIntElt : r := 1) -> .
    {Given a univariate polynomial f with coefficients in F_p^r (or coercible in F_p^r), return its Cartier-Manin matrix.}
    if Type(BaseRing(Parent(f))) ne FldFin then
        assert IsCoercible(PolynomialRing(GF(p^r)), f);
        f := PolynomialRing(GF(p^r))!f;
    end if;
    g := (Degree(f)-1) div 2;
    h := f^((p-1) div 2);
    return Matrix([[p*j-i ge 0 select Coefficient(h, p*j-i) else 0 : j in [1..g]] : i in [1..g]]);
end intrinsic;

intrinsic StableHasseWitt(M::Any, p::RngIntElt : r := 1) -> .
    {}
    if r eq 1 then 
        return M^(Ncols(M));
    end if;
    return &*([M] cat [Matrix([[M[i,j]^(p^k) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]) : k in [1..Ncols(M)-1]]);
end intrinsic;

intrinsic IsSupersingular(f::Any, p::RngIntElt, r::RngIntElt) -> Bool
    {Check whether the hyperelliptic curve defined by f over F_p^r is supersingular or not.}
    Lpol := Coefficients(LPolynomial(HyperellipticCurve(f)));
    if &and[Valuation(Lpol[i], p) ge r*(i-1)/2 : i in [1..#Lpol]] then
        return true;
    end if;
    return false;
end intrinsic;

intrinsic IsomorphicNormalRepresentative(f::RngUPolElt) -> RngUPolElt
    {Returns the "smallest" polynomial P (in the Magma sense) in our family which defines a hyperelliptic curve y^2 = P(x) which is isomorphic over F_q to y^2 = f(x).}
    R<x> := Parent(f);
    d := Degree(f);
    if d mod Characteristic(BaseRing(R)) eq 0 then
        printf "Problem of characteristic.\n";
        return f;
    end if;
    roots := Roots(f);
    res := f;
    for root in roots do
        if root[2] eq 1 then
            r := root[1];
            f1 := R!((-x+r)^(d+1)*Evaluate(f, (-r*x+(r^2-1))/(-x+r)));
            if Coefficient(f1, d) ne 0 then
                f1 /:= Coefficient(f1, d);
            else 
                "Problem, this coefficient shouldn't vanish.";break;
            end if;
            f1 := Evaluate(f1, x-Coefficient(f1, d-1)/d);
            if Coefficient(f1, d-2)*Coefficient(f1, d-3) ne 0 then
                u := Coefficient(f1, d-2)/Coefficient(f1, d-3);
                res := Min([res, u^d*Evaluate(f1, x/u)]);
            end if;
        end if;
    end for;
    return res;
end intrinsic;

intrinsic CurveComplete(g::RngIntElt, p::RngIntElt, r::RngIntElt) -> SeqEnum
    {Returns a list of all the curves of genus g in our family over F_p^r.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating curves...";
    if g eq 2 then
        time PolyList := [F![a1,a2,a3,a3,0,1] : a1,a2,a3 in Fp | a3 ne 0];
    elif g eq 3 then
        time PolyList := [F![a1,a2,a3,a4,a5,a5,0,1] : a1,a2,a3,a4,a5 in Fp | a5 ne 0];
    elif g eq 4 then
        time PolyList := [F![a1,a2,a3,a4,a5,a6,a7,a7,0,1] : a1,a2,a3,a4,a5,a6,a7 in Fp | a7 ne 0];
    elif g eq 5 then
        time PolyList := [F![a1,a2,a3,a4,a5,a6,a7,a8,a9,a9,0,1] : a1,a2,a3,a4,a5,a6,a7,a8,a9 in Fp | a9 ne 0];
    else 
        "Not implemented, but it is easy to do";
        return [];
    end if;
    return PolyList;    
end intrinsic;

intrinsic CurveSamples(g::RngIntElt, p::RngIntElt, r::RngIntElt, s::RngIntElt) -> SeqEnum
    {Returns a list of around s random curves of genus g over F_p^r in our family.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;

    "Creating curves...";
    time TupleSet := Set([F!([Random(Fp) : j in [1..2*g-2]] cat [rand,rand,0,1]) where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);

    if #TupleSet lt s then
        printf "Only picked up %o curves instead of %o\n", #TupleSet, s;
    end if;

    return Setseq(TupleSet);
end intrinsic;


intrinsic CurveDataForP(g::RngIntElt, p::RngIntElt, r::RngIntElt : s := 0) -> SeqEnum
    {Computes the p-rank of s elements of the family of genus g we consider, with coefficients in F_p^r. If s is 0, returns the data for the whole family.}

    if (s eq 0) or (s ge (p^r)^(2*g-1)) then
        PolyList := CurveComplete(g, p, r);
    else
        PolyList := CurveSamples(g, p, r, Floor(1.1*s));
    end if;

    "Checking which curves are smooth...";
    time PolyList := [f : f in PolyList | Discriminant(f) ne 0];

    "Removing isomorphic curves...";
    time PolyList := Setseq(Set([IsomorphicNormalRepresentative(f) : f in PolyList]));    
    
    if #PolyList gt s and s gt 0 then
        PolyList := PolyList[1..s];
    elif s gt 0 then
        printf "Only picked up %o curves instead of %o\n", #PolyList, s;
    end if;

    "Computing p-ranks...";
    time p_ranks := [Rank(StableHasseWitt(HasseWittMatrix(f, p : r := r), p : r := r)) : f in PolyList];

    res := [0 : i in [1..g+1]];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

/*
intrinsic NotRoot(Fp::FldFin, s::RngIntElt) -> .
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
*/
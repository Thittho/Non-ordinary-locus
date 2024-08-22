intrinsic CartierMatrix(f::RngUPolElt, p::RngIntElt : r := 1) -> ModMatFldElt
    {Given a univariate polynomial f with coefficients in F_p^r (or coercible in F_p^r), return its Cartier-Manin matrix.}
    
    assert IsCoercible(PolynomialRing(GF(p^r)), f);

    f := PolynomialRing(GF(p^r))!f;
    g := (Degree(f)-1) div 2;
    h := f^((p-1) div 2);
    return Matrix([[p*i-j ge 0 select Coefficient(h, p*i-j) else 0 : i in [1..g]] : j in [1..g]]);
end intrinsic;


intrinsic CurveDataForP(p::RngIntElt, r::RngIntElt, s::RngIntElt : i0 := "not def") -> SeqEnum
    {Computes the p-rank of s elements of the family we consider, with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if Type(i0) eq MonStgElt then
        TupleSet := Set([<Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp)> : i in [1..Floor(5*s/4)]]);
    else
        TupleSet := Set([<Fp!i0,Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp)> : i in [1..Floor(5*s/4)]]);
    end if;

    TupleList := Setseq(TupleSet);

    PolyList10 := [x^9+x^7+x^6+(&+[t[i+1]*x^i : i in [0..5]]) : t in TupleList];

    time PolyList10_2 := [f : f in PolyList10 | Discriminant(f) ne 0];

    if #PolyList10_2 gt s then
        PolyList10_2 := PolyList10_2[1..s];
    else 
        "Only picked up", #PolyList10_2, "curves instead of", s;
    end if;

    p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList10_2];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPComplete(p::RngIntElt, r::RngIntElt : i0 := "not def") -> SeqEnum
    {Computes the p-rank of every member of the family with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    if Type(i0) eq MonStgElt then    
        Fp := GF(p^r);
        F<x> := PolynomialRing(Fp);
        TupleList := [<a1,a2,a3,a4,a5,a6> : a1,a2,a3,a4,a5,a6 in Fp];
    else
        if Type(Parent(i0)) eq FldFin then
            Fp := Parent(i0);
        else 
            Fp := GF(p^r);
        end if;
        F<x> := PolynomialRing(Fp);
        TupleList := [<Fp!i0,a2,a3,a4,a5,a6> : a2,a3,a4,a5,a6 in Fp];
    end if;

    PolyList10 := [x^9+x^7+x^6+(&+[t[i+1]*x^i : i in [0..5]]) : t in TupleList];

    time PolyList10_2 := [f : f in PolyList10 | Discriminant(f) ne 0];

    p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList10_2];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


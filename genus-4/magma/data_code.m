intrinsic CartierMatrix(f::RngUPolElt, p::RngIntElt : r := 1) -> ModMatFldElt
    {Given a univariate polynomial f with coefficients in F_p^r (or coercible in F_p^r), return its Cartier-Manin matrix.}
    
    if Type(BaseRing(Parent(f))) ne FldFin then
        assert IsCoercible(PolynomialRing(GF(p^r)), f);
        f := PolynomialRing(GF(p^r))!f;
    end if;

    g := (Degree(f)-1) div 2;
    h := f^((p-1) div 2);
    return Matrix([[p*i-j ge 0 select Coefficient(h, p*i-j) else 0 : i in [1..g]] : j in [1..g]]);
end intrinsic;

function StableCartierManin(M, p : r := 1)
    if r eq 1 then 
        return M^(Ncols(M));
    end if;
    return &*([M] cat [Matrix([[M[i,j]^(p^k) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]) : k in [1..Ncols(M)-1]]);
end function;

intrinsic CurveDataForPOld(p::RngIntElt, r::RngIntElt, s::RngIntElt : i0 := "not def") -> SeqEnum
    {Computes the p-rank of s elements of the family we consider, with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating curves...";
    if Type(i0) eq MonStgElt then
        time TupleSet := Set([<Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp)> : i in [1..Floor(5*s/4)]]);
    else
        time TupleSet := Set([<Fp!i0,Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp)> : i in [1..Floor(5*s/4)]]);
    end if;

    TupleList := Setseq(TupleSet);

    time PolyList10 := [x^9+x^7+x^6+(&+[t[i+1]*x^i : i in [0..5]]) : t in TupleList];

    "Checking which are smooth...";
    time PolyList10_2 := [f : f in PolyList10 | Discriminant(f) ne 0];
    
    if #PolyList10_2 gt s then
        PolyList10_2 := PolyList10_2[1..s];
    else 
        "Only picked up", #PolyList10_2, "curves instead of", s;
    end if;
    
    "Computing p-ranks...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList10_2];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPCompleteOld(p::RngIntElt, r::RngIntElt : i0 := "not def") -> SeqEnum
    {Computes the p-rank of every member of the family with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    "Creating curves...";
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

    time PolyList10 := [x^9+x^7+x^6+(&+[t[i+1]*x^i : i in [0..5]]) : t in TupleList];

    "Checking which curves are smooth...";
    time PolyList10_2 := [f : f in PolyList10 | Discriminant(f) ne 0];

    "Computing p-ranks...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList10_2];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForP(p::RngIntElt, r::RngIntElt, s::RngIntElt : i0 := "not def") -> SeqEnum
    {Computes the p-rank of s elements of the family we consider, with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;

    "Creating list of curves...";
    if Type(i0) eq MonStgElt  then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
         time TupleSet := Set([F![Fp!i0,Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    end if;

    PolyList := Setseq(TupleSet);
    
    "Checking which curves are smooth...";
    time PolyList := [f : f in PolyList | Discriminant(f) ne 0];
    
    if #PolyList gt s then
        PolyList := PolyList[1..s];
    else 
        "Only picked up", #PolyList, "curves instead of", s;
    end if;

    "Computing p-ranks...";
    //time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList];
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList];

    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPComplete(p::RngIntElt, r::RngIntElt : i0 := "not def", j0 := "not def") -> SeqEnum
    {Computes the p-rank of every member of the family with coefficients in F_p^r. The parameter i0 can be defined to be an integer or element of F_p^r, and in that case fixes the value of the constant coefficient of the polynomials considered to be i0.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating all curves...";
    if Type(i0) eq MonStgElt and Type(j0) eq MonStgElt then    
        time PolyList := [F![a1,a2,a3,a4,a5,a6,a7,a7,0,1] : a1,a2,a3,a4,a5,a6,a7 in Fp | a7 ne 0];
    elif Type(i0) ne MonStgElt and Type(j0) eq MonStgElt then
        time PolyList := [F![Fp!i0,a2,a3,a4,a5,a6,a7,a7,0,1] : a2,a3,a4,a5,a6,a7 in Fp | a7 ne 0];
    elif Type(j0) ne MonStgElt and Type(i0) eq MonStgElt then
        time PolyList := [F![a2,Fp!j0,a3,a4,a5,a6,a7,a7,0,1] : a2,a3,a4,a5,a6,a7 in Fp | a7 ne 0];
    else
        time PolyList := [F![Fp!i0,Fp!j0,a3,a4,a5,a6,a7,a7,0,1] : a3,a4,a5,a6,a7 in Fp | a7 ne 0];
    end if;

    "Computing p-ranks of smooth curves...";
    //time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList | Discriminant(f) ne 0];

    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

/*
F<x> := PolynomialRing(GF(5));
Fp<alpha> := ext<GF(5) | x^3+3*x+3>;
Fpx<y> := PolynomialRing(Fp);
f := y^5 + y^4 + alpha^92*y^3 + alpha^18*y^2 + alpha^56*y;
Attach("~/github/Non-ordinary-locus/magma/data_code.m");
M := Transpose(CartierMatrix(f, 5));
M1 := Matrix([[M[i,j]^5 : j in [1..2]] : i in [1..2]]);
p := 5;
r := 3;
StableCartierManin(Transpose(CartierMatrix(f, p : r := r)), p : r := r);
M*M1;
*/
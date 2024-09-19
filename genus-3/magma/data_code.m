

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
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
         time TupleSet := Set([F![Fp!i0,Random(Fp),Random(Fp),Random(Fp),rand,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
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

    res := [0,0,0,0];
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
        time PolyList := [F![a1,a2,a3,a4,a5,a5,0,1] : a1,a2,a3,a4,a5 in Fp | a5 ne 0];
    elif Type(i0) ne MonStgElt and Type(j0) eq MonStgElt then
        time PolyList := [F![Fp!i0,a2,a3,a4,a5,a5,0,1] : a2,a3,a4,a5 in Fp | a5 ne 0];
    elif Type(j0) ne MonStgElt and Type(i0) eq MonStgElt then
        time PolyList := [F![a2,Fp!j0,a3,a4,a5,a5,0,1] : a2,a3,a4,a5 in Fp | a5 ne 0];
    else
        time PolyList := [F![Fp!i0,Fp!j0,a3,a4,a5,a5,0,1] : a3,a4,a5 in Fp | a5 ne 0];
    end if;

    "Computing p-ranks of smooth curves...";
    //time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList | Discriminant(f) ne 0];

    res := [0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPCompleteDim4(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 4-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1 or 2.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        alpha := NotRoot(Fp, 2);
        time PolyList1 := Set([Min([f, -Evaluate(f, -x)]) where f := F![a1,a2,a3,a4,0,a4,0,1] : a1,a2,a3,a4 in Fp | a4 ne 0]);
        time PolyList2 := Set([Min([f, -Evaluate(f, -x)]) where f := F![a1,a2,a3,alpha*a4,0,a4,0,1] : a1,a2,a3,a4 in Fp | a4 ne 0 and (a1 ne 0 or a3 ne 0)]); // otherwise we overcount the curves with gcd 2
        PolyList := Setseq(PolyList1 join PolyList2);
    elif fam eq 2 then
        time PolyList := [F![a1,a2,a3,a4,a4,0,0,1] : a1,a2,a3,a4 in Fp | a4 ne 0];
    else
        "fam parameter should be 1 or 2";
    end if;

    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPCompleteDim3(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 3-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2 or 3.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        if #Roots(x^3-1) eq 1 then // bijection
            time PolyList := [F![a1,a2,a3,0,0,a3,0,1] : a1,a2,a3 in Fp | a3 ne 0];
        else
            roots := [s[1] : s in Roots(x^3-1)];
            assert #roots eq 3;
            alpha := NotRoot(Fp, 3);
            time PolyList1 := Set([Min([roots[i]^7*Evaluate(f, x/roots[i]) : i in [1..#roots]]) where f := F![a1,a2,a3,0,0,a3,0,1] : a1,a2,a3 in Fp | a3 ne 0]);
            time PolyList2 := Set([Min([roots[i]^7*Evaluate(f, x/roots[i]) : i in [1..#roots]]) where f := F![a1,a2,alpha*a3,0,0,a3,0,1] : a1,a2,a3 in Fp | a3 ne 0]);
            time PolyList3 := Set([Min([roots[i]^7*Evaluate(f, x/roots[i]) : i in [1..#roots]]) where f := F![a1,a2,alpha^2*a3,0,0,a3,0,1] : a1,a2,a3 in Fp | a3 ne 0]);
            PolyList := Setseq(PolyList1 join PolyList2 join PolyList3);
        end if;
    elif fam eq 2 then
        alpha := NotRoot(Fp, 2);
        time PolyList1 := Set([Min([f, -Evaluate(f, -x)]) where f := F![a1,a2,a3,0,a3,0,0,1] : a1,a2,a3 in Fp | a3 ne 0]);
        time PolyList2 := Set([Min([f, -Evaluate(f, -x)]) where f := F![a1,a2,alpha*a3,0,a3,0,0,1] : a1,a2,a3 in Fp | a3 ne 0]);
        PolyList := Setseq(PolyList1 join PolyList2);
    elif fam eq 3 then
        time PolyList := [F![a1,a2,a3,a3,0,0,0,1] : a1,a2,a3 in Fp | a3 ne 0];
    else
        "fam parameter should be 1 or 2";
    end if;

    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

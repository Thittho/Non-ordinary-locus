function RemoveDuplicates(TupleSet0, m : pain := 0)
    //#TupleSet0;
    F<x> := Universe(TupleSet0);
    
    roots := [r[1] : r in Roots(F!(x^m-1))];
    if #roots eq 1 then
        #TupleSet0;
        return TupleSet0; 
    end if;

    roots;
    TupleSetSym := Set(&cat[[r^9*Evaluate(f, 1/r*x) : r in roots | r^9*Evaluate(f, 1/r*x) ne f] : f in TupleSet0]);
    //#TupleSetSym;
    TupleSet := TupleSet0 diff TupleSetSym;
    res := TupleSet0 meet TupleSetSym;
    res2 := res;
    for f in res do
        if f in res2 then
            Include(~TupleSet, f);
            for r in roots do
                Exclude(~res2, r^9*Evaluate(f, 1/r*x));
            end for;
        end if;
    end for;

    if m eq 4 and #roots eq 2 and pain eq 1 then // there is a specific case where we overcount curves
        TupleSetSym := {f : f in TupleSet | -Evaluate(f, -x) eq f};
        for f in TupleSetSym do
            if f in TupleSet then
                Exclude(~TupleSet, x^9-Coefficient(f,7)-Coefficient(f,3)+Coefficient(f,1));
            end if;
        end for;
    end if;
    //#TupleSet;

    return TupleSet;
end function;

intrinsic CurveDataForPDim6(p::RngIntElt, r::RngIntElt, s::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 6-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1 or 2.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;

    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,0,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),alpha*rand,0,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]); 
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 2 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,rand,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
        "fam parameter should be 1 or 2";
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
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

intrinsic CurveDataForPDim5(p::RngIntElt, r::RngIntElt, s::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 5-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2 or 3.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;
    
    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,0,0,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 2 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,0,rand,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 3 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),Random(Fp),rand,rand,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
        "fam parameter should be 1, 2 or 3";
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
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

intrinsic CurveDataForPDim4(p::RngIntElt, r::RngIntElt, s::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 4-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2, 3 or 4.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;
    
    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),rand,0,0,0,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 4 : pain := 1);
    elif fam eq 2 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),rand,0,0,rand,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]); 
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 3 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),rand,0,rand,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 4 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),Random(Fp),rand,rand,0,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
        "fam parameter should be 1, 2, 3 or 4";
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
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPDim3(p::RngIntElt, r::RngIntElt, s::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 4-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2, 3, 4 or 5.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    if r eq 1 then 
        Fpmult := [1..p-1];
    else
        Fpmult := [(Fp.1)^j : j in [0..#MultiplicativeGroup(Fp)-1]];
    end if;
    
    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),rand,0,0,0,0,rand,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 5);
    elif fam eq 2 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),rand,0,0,0,rand,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 4);
    elif fam eq 3 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),rand,0,0,rand,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 4 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),rand,0,rand,0,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 5 then
        time TupleSet := Set([F![Random(Fp),Random(Fp),rand,rand,0,0,0,0,0,1] where rand := Random(Fpmult) : i in [1..Floor(5*s/4)]]);
    else
        "fam parameter should be 1, 2, 3, 4 or 5";
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
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


















/////////////////////////////////////////////////////////////////////////////

intrinsic CurveDataForPCompleteDim6(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 6-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1 or 2.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        alpha := NotRoot(Fp, 2);
        time PolyList1 := Set([Min([f, -Evaluate(f, -x)]) where f :=F![a1,a2,a3,a4,a5,a6,0,a6,0,1] : a1,a2,a3,a4,a5,a6 in Fp | a6 ne 0]);
        time PolyList2 := Set([Min([f, -Evaluate(f, -x)]) where f :=F![a1,a2,a3,a4,a5,alpha*a6,0,a6,0,1] : a1,a2,a3,a4,a5,a6 in Fp | a6 ne 0 and (a1 ne 0 or a3 ne 0 or a5 ne 0)]);
        PolyList := Setseq(PolyList1 join PolyList2);
    elif fam eq 2 then
        time PolyList := [F![a1,a2,a3,a4,a5,a6,a6,0,0,1] : a1,a2,a3,a4,a5,a6 in Fp | a6 ne 0];
    else
        "fam parameter should be 1 or 2";
    end if;
    
    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

intrinsic CurveDataForPCompleteDim5(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 5-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2 or 3.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![a2,a3,a4,a5,a6,0,0,a6,0,1] : a2,a3,a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 2 then
        time TupleSet := Set([F![a2,a3,a4,a5,a6,0,a6,0,0,1] : a2,a3,a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 3 then
        time TupleSet := Set([F![a2,a3,a4,a5,a6,a6,0,0,0,1] : a2,a3,a4,a5,a6 in Fp | a6 ne 0]);
    else
        "fam parameter should be 1, 2 or 3";
    end if;

    PolyList := Setseq(TupleSet);
    
    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;

intrinsic CurveDataForPCompleteDim4(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 4-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2, 3 or 4.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![a3,a4,a5,a6,0,0,0,a6,0,1] : a3,a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 4);
    elif fam eq 2 then
        time TupleSet := Set([F![a3,a4,a5,a6,0,0,a6,0,0,1] : a3,a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 3 then
        time TupleSet := Set([F![a3,a4,a5,a6,0,a6,0,0,0,1] : a3,a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 4 then
        time TupleSet := Set([F![a3,a4,a5,a6,a6,0,0,0,0,1] : a3,a4,a5,a6 in Fp | a6 ne 0]);
    else
        "fam parameter should be 1, 2, 3 or 4";
    end if;

    PolyList := Setseq(TupleSet);
    
    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;


intrinsic CurveDataForPCompleteDim3(p::RngIntElt, r::RngIntElt : fam := 1) -> SeqEnum
    {Computes the p-rank of s elements of the 4-dim families we consider, with coefficients in F_p^r. The parameter fam is related to the family considered, it must be 1, 2, 3, 4 or 5.}
    
    Fp := GF(p^r);
    F<x> := PolynomialRing(Fp);

    "Creating list of curves...";

    if fam eq 1 then
        time TupleSet := Set([F![a4,a5,a6,0,0,0,0,a6,0,1] : a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 5);
    elif fam eq 2 then
        time TupleSet := Set([F![a4,a5,a6,0,0,0,a6,0,0,1] : a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 4);
    elif fam eq 3 then
        time TupleSet := Set([F![a4,a5,a6,0,0,a6,0,0,0,1] : a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 3);
    elif fam eq 4 then
        time TupleSet := Set([F![a4,a5,a6,0,a6,0,0,0,0,1] : a4,a5,a6 in Fp | a6 ne 0]);
        TupleSet := RemoveDuplicates(TupleSet, 2);
    elif fam eq 5 then
        time TupleSet := Set([F![a4,a5,a6,a6,0,0,0,0,0,1] : a4,a5,a6 in Fp | a6 ne 0]);
    else
        "fam parameter should be 1, 2, 3, 4 or 5";
    end if;

    PolyList := Setseq(TupleSet);
    
    "Computing p-ranks of smooth curves...";
    time p_ranks := [Rank(CartierMatrix(f, p : r := r)^4) : f in PolyList | Discriminant(f) ne 0];
    
    res := [0,0,0,0,0];
    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end intrinsic;
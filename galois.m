// single-galois-type.magma
// Version 2.0
// November 2024

// Thomas Bouchet, Erik Davis, Steven Groen, Zachary Porat, Benjamin York

/*
================================================================================
// Loading Everett Howe's genus-2 and genus-3 code.
================================================================================
*/

load "~/github/hyperelliptic/Hyperelliptic3.magma";
Attach("tools.m");

/*
================================================================================
// Finding hyperelliptic curves with Galois type (1,1,...,1)
================================================================================
*/

function sample_eightpoints(K, sample)
    // York: Modified original method to only find a number of hyperelliptic curves
    // equal to variable sample
    //
    // The following is a modification of Howe's documentation for the new setting
    //
    // Howe: Find PGL2(K) orbit reps for eight points of P^1(K). 
    //
    // What is our normal form?
    //   Given eight points, put one at infinity, one at 0, and ask that the sum
    //   of the rest be 1. This may not be possible for the given choice of points
    //   to put at oo and 0... but it is possible for *some* choice, with one
    //   exception: When K = GF(7) and the eight points are the elements 
    //   of P^1(K). We treat that one case separately.
    // 
    // So: Loop a1 < a2 < ... < a5 and set a6 = 1 - a1 - a2 - a3 - a4 - a5.
    // Demand that a5 < a6.
    //
    // In turn, shift each a_i to 0 and readjust scaling. Demand that our original
    // a1 through a6 be smaller than these new values.
    //
    // In turn, move each a_i to infinity, and repeat same loop as above, 
    // and again demand that our original a1 through a6 are the smallest.
    //
    // In Howe's original code, the automorphism group of the hyperelliptic curve was tracked
    // We have modified the code to eliminate curves with non-trivial automorphisms from the data set
    //

    q := #K;
    if q lt 7 then return []; end if;
    sortKstar := Sort([a : a in K | a ne 0]);
    R<x> := PolynomialRing(K);
    polys := [];
    counter := 0;
  
    
    for i1 in [1..q-6] do
        a1 := sortKstar[i1];
        for i2 in [i1+1..q-5] do
            a2 := sortKstar[i2];
            for i3 in [i2+1..q-4] do
                a3 := sortKstar[i3];
                for i4 in [i3+1..q-3] do
                    a4 := sortKstar[i4];
                    for i5 in [i4+1..q-2] do
                        a5 := sortKstar[i5];
                        a6 := 1 - a1 - a2 - a3 - a4 - a5;
                        if a5 ge a6 then continue i5; end if;
            
                        basecase := Sort([a1,a2,a3,a4,a5,a6]);
            
                        for i in [1..6] do
                            if 1 ne 7*basecase[i] then
                                s := 1/(1-7*basecase[i]);
                                newcase := Sort([-basecase[i]*s] cat [(basecase[j]-basecase[i])*s : j in [1..6] | j ne i]);
                                if newcase le basecase then continue i5; end if;
                            end if;
                        end for;
            
                        // Put 0 at oo.
            
                        newbase := [1/basecase[i] : i in [1..6]];
                
                        sum := &+newbase;
                        if sum ne 0 then
                            s := 1/sum;
                            newcase := Sort([a*s : a in newbase]);
                            if newcase le basecase then continue i5; end if;
                        end if;
                
                        for j in [1..6] do
                            if sum ne 7*newbase[j] then
                                s := 1/(sum-7*newbase[j]);
                                newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..6] | k ne j]);
                                if newcase le basecase then continue i5; end if;
                            end if;
                        end for;

                        for i in [1..6] do
                        // send ai to oo and oo to 0... x --> 1/(x-ai)
                            newbase := [-1/basecase[i]] cat [1/(basecase[j]-basecase[i]) : j in [1..6] | j ne i];
                            
                            sum := &+newbase;
                            if sum ne 0 then
                                s := 1/sum;
                                newcase := Sort([a*s : a in newbase]);
                                if newcase le basecase then continue i5; end if;
                            end if;
                            
                            for j in [1..6] do
                                if sum ne 7*newbase[j] then
                                    s := 1/(sum-7*newbase[j]);
                                    newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..6] | k ne j]);
                                    if newcase le basecase then continue i5; end if;
                                end if;
                            end for;
                        end for;
                        
                        // We got one!
                        
                        f := x*(x-a1)*(x-a2)*(x-a3)*(x-a4)*(x-a5)*(x-a6);
                        polys cat:= [f];

                        // if number of curves found is greater than sample, end loop
                        counter +:= 1;
                        if counter ge sample then
                            break i1;
                        end if;
                
                    end for;
                end for;
            end for;
        end for;
    end for;
    
    return polys;
end function;




function sample_tenpoints(K, sample)
    // York: Modified original method to only find a number of hyperelliptic curves
    // equal to variable sample
    // 
    // The following is a modification of Howe's documentation for the new setting
    //
    // Find PGL2(K) orbit reps for ten points of P^1(K).
    //
    // What is our normal form?
    //   Given eight points, put one at infinity, one at 0, and ask that the sum
    //   of the rest be 1. This may not be possible for the given choice of points
    //   to put at oo and 0, but it is possible for *some* choice, 
    //   with one exception that we exclude at the start
    // 
    // So: Loop a1 < a2 < ... < a6 < a7 and set a8 = 1 - a1 - a2 - a3 - a4 - a5 - a6 - a7
    // Demand that a7 < a8.
    //
    // In turn, shift each a_i to 0 and readjust scaling. Demand that our original
    // a1 through a8 be smaller than these new values.
    //
    // In turn, move each a_i to infinity, and repeat same loop as above, 
    // and again demand that our original a1 through a8 are the smallest.
    //
    // In Howe's original code, the automorphism group of the hyperelliptic curve was tracked
    // We have modified the code to eliminate curves with non-trivial automorphisms from the data set
    //

    q := #K;
    if q lt 9 then return []; end if;
    sortKstar := Sort([a : a in K | a ne 0]);
    R<x>:=PolynomialRing(K);
    polys := [];
    counter := 0;
  
    
    for i1 in [1..q-8] do
        a1 := sortKstar[i1];
        for i2 in [i1+1..q-7] do
            a2 := sortKstar[i2];
            for i3 in [i2+1..q-6] do
                a3 := sortKstar[i3];
                for i4 in [i3+1..q-5] do
                    a4 := sortKstar[i4];
                    for i5 in [i4+1..q-4] do
                        a5 := sortKstar[i5];
                        for i6 in [i5+1..q-3] do
                            a6 := sortKstar[i6];
                            for i7 in [i6+1..q-2] do
                                a7 := sortKstar[i7];
                                a8 := 1 - a1 - a2 - a3 - a4 - a5 - a6 - a7;
                                if a7 ge a8 then continue i7; end if;
                                  
                                basecase := Sort([a1,a2,a3,a4,a5,a6,a7,a8]);
                    
                                for i in [1..8] do
                                    if 1 ne 9*basecase[i] then
                                        s := 1/(1-9*basecase[i]);
                                        newcase := Sort([-basecase[i]*s] cat [(basecase[j]-basecase[i])*s : j in [1..8] | j ne i]);
                                        if newcase le basecase then continue i7; end if;
                                    end if;
                                end for;
                    
                                // Put 0 at oo.
                    
                                newbase := [1/basecase[i] : i in [1..8]];
                        
                                sum := &+newbase;
                                if sum ne 0 then
                                    s := 1/sum;
                                    newcase := Sort([a*s : a in newbase]);
                                    if newcase le basecase then continue i7; end if;
                                end if;
                        
                                for j in [1..8] do
                                    if sum ne 9*newbase[j] then
                                        s := 1/(sum-9*newbase[j]);
                                        newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..8] | k ne j]);
                                        if newcase le basecase then continue i7; end if;
                                    end if;
                                end for;

                                for i in [1..8] do
                                // send ai to oo and oo to 0... x --> 1/(x-ai)
                                    newbase := [-1/basecase[i]] cat [1/(basecase[j]-basecase[i]) : j in [1..8] | j ne i];
                                    
                                    sum := &+newbase;
                                    if sum ne 0 then
                                        s := 1/sum;
                                        newcase := Sort([a*s : a in newbase]);
                                        if newcase le basecase then continue i7; end if;
                                    end if;
                                    
                                    for j in [1..8] do
                                        if sum ne 9*newbase[j] then
                                            s := 1/(sum-9*newbase[j]);
                                            newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..8] | k ne j]);
                                            if newcase le basecase then continue i7; end if;
                                        end if;
                                    end for;
                                end for;
                                
                                // We got one!
                                
                                f := x*(x-a1)*(x-a2)*(x-a3)*(x-a4)*(x-a5)*(x-a6)*(x-a7)*(x-a8);
                                polys cat:= [f];

                                // if number of curves found is greater than sample, end loop
                                counter +:= 1;
                                if counter ge sample then
                                    break i1;
                                end if;
                            end for;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return polys;
end function;


function sample_twelvepoints(K, sample)
    // York: Modified original function to only find a number of hyperelliptic curves
    // equal to variable sample
    // 
    // The following is a modification of Howe's documentation for the new setting
    //
    // Find PGL2(K) orbit reps for twelve points of P^1(K).
    //
    // What is our normal form?
    //   Given twelve points, put one at infinity, one at 0, and ask that the sum
    //   of the rest be 1. This may not be possible for the given choice of points
    //   to put at oo and 0, but it is possible for *some* choice, with one exception
    //   that we eliminate at the start
    // 
    // So: Loop a1 < a2 < ... < a8 < a9 and set a10 = 1 - a1 - a2 - a3 - a4 - a5 - a6 - a7 - a8 - a9.
    // Demand that a9 < a10.
    //
    // In turn, shift each a_i to 0 and readjust scaling. Demand that our original
    // a1 through a10 be smaller than these new values.
    //
    // In turn, move each a_i to infinity, and repeat same loop as above, 
    // and again demand that our original a1 through a10 are the smallest.
    //
    // In Howe's original code, the automorphism group of the hyperelliptic curve was tracked
    // We have modified the code to eliminate curves with non-trivial automorphisms from the data set
    //

    q := #K;
    if q lt 11 then return []; end if;
    sortKstar := Sort([a : a in K | a ne 0]);
    R<x>:=PolynomialRing(K);
    polys := [];
    counter := 0;
  
    for i1 in [1..q-10] do
        a1 := sortKstar[i1];
        for i2 in [i1+1..q-9] do
            a2 := sortKstar[i2];
            for i3 in [i2+1..q-8] do
                a3 := sortKstar[i3];
                for i4 in [i3+1..q-7] do
                    a4 := sortKstar[i4];
                    for i5 in [i4+1..q-6] do
                        a5 := sortKstar[i5];
                        for i6 in [i5+1..q-5] do
                            a6 := sortKstar[i6];
                            for i7 in [i6+1..q-4] do
                                a7 := sortKstar[i7];
                                for i8 in [i7+1..q-3] do
                                    a8 := sortKstar[i8];
                                    for i9 in [i8+1..q-2] do
                                        a9 := sortKstar[i9];
                                        a10 := 1 - a1 - a2 - a3 - a4 - a5 - a6 - a7 - a8 - a9;
                                        if a9 ge a10 then continue i9; end if;
                            
                                        basecase := Sort([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10]);
                            
                                        for i in [1..10] do
                                            if 1 ne 11*basecase[i] then
                                                s := 1/(1-11*basecase[i]);
                                                newcase := Sort([-basecase[i]*s] cat [(basecase[j]-basecase[i])*s : j in [1..10] | j ne i]);
                                                if newcase le basecase then continue i9; end if;
                                            end if;
                                        end for;
                            
                                        // Put 0 at oo.
                            
                                        newbase := [1/basecase[i] : i in [1..10]];
                                
                                        sum := &+newbase;
                                        if sum ne 0 then
                                            s := 1/sum;
                                            newcase := Sort([a*s : a in newbase]);
                                            if newcase le basecase then continue i9; end if;
                                        end if;
                                
                                        for j in [1..10] do
                                            if sum ne 11*newbase[j] then
                                                s := 1/(sum-11*newbase[j]);
                                                newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..10] | k ne j]);
                                                if newcase le basecase then continue i9; end if;
                                            end if;
                                        end for;

                                        for i in [1..10] do
                                        // send ai to oo and oo to 0... x --> 1/(x-ai)
                                            newbase := [-1/basecase[i]] cat [1/(basecase[j]-basecase[i]) : j in [1..10] | j ne i];
                                            
                                            sum := &+newbase;
                                            if sum ne 0 then
                                                s := 1/sum;
                                                newcase := Sort([a*s : a in newbase]);
                                                if newcase le basecase then continue i9; end if;
                                            end if;
                                            
                                            for j in [1..10] do
                                                if sum ne 11*newbase[j] then
                                                    s := 1/(sum-11*newbase[j]);
                                                    newcase := Sort([-newbase[j]*s] cat [(newbase[k]-newbase[j])*s : k in [1..10] | k ne j]);
                                                    if newcase le basecase then continue i9; end if;
                                                end if;
                                            end for;
                                        end for;
                                        
                                        // We got one!
                                        
                                        f := x*(x-a1)*(x-a2)*(x-a3)*(x-a4)*(x-a5)*(x-a6)*(x-a7)*(x-a8)*(x-a9)*(x-a10);
                                        polys cat:= [f];

                                        // if number of curves found is greater than sample, end loop
                                        counter +:= 1;
                                        if counter ge sample then
                                            break i1;
                                        end if;
                                    end for;
                                end for;
                            end for;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return polys;
end function;

/*
================================================================================
// Computing the Cartier-Manin Matrix
================================================================================
*/

function CartierMatrix(f, p : r := 1)
    if Type(BaseRing(Parent(f))) ne FldFin then
        assert IsCoercible(PolynomialRing(GF(p^r)), f);
        f := PolynomialRing(GF(p^r))!f;
    end if;
    
    g := (Degree(f)-1) div 2;
    h := f^((p-1) div 2);
    return Matrix([[p*i-j ge 0 select Coefficient(h, p*i-j) else 0 : i in [1..g]] : j in [1..g]]);
end function;

function StableCartierManin(M, p : r := 1)
    if r eq 1 then 
        return M^(Ncols(M));
    end if;
    return &*([M] cat [Matrix([[M[i,j]^(p^k) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]) : k in [1..Ncols(M)-1]]);
end function;


/*
================================================================================
// Code modified from Everett Howe finding hyperelliptic curves with Galois type (1,1,...,1)
================================================================================
*/

function CurveDataForPGalois(g, p, r, s)
    // Given a genus g = 3, 4, 5
    // we compute s hyperelliptic curves of genus g defined over F_{p^r}
    // and then compute their p-ranks via the Cartier-Manin matrix
    // we report this data as an ordered list of form 
    // [p-rank 0, p-rank 1, ...]

    Fq := GF(p^r);
    F<x> := PolynomialRing(Fq);

    "Creating list of curves...";
    
    if g eq 3 then
        time PolyList := sample_eightpoints(Fq, s);
    elif g eq 4 then
        time PolyList := sample_tenpoints(Fq, s);
    elif g eq 5 then
        time PolyList := sample_twelvepoints(Fq, s);
    else
        PolyList := [];
        "Genus is unsupported";
    end if;
    
    if #PolyList lt s then 
        "Only picked up", #PolyList, "curves instead of", s;
    end if;

    "Computing p-ranks...";
    time p_ranks := [Rank(StableCartierManin(CartierMatrix(f, p : r := r), p : r := r)) : f in PolyList];

    res := [0]; //make list of zeroes of length g+1
    for k in [1..g] do
        res cat:= [0];
    end for;

    for i in p_ranks do
        res[i+1] +:= 1;
    end for;

    return res;
end function;


/*

//Example testing hyperelliptic curve generating functions for genus g = 3, 4, 5
for g in [3,4,5] do
    for p in PrimesInInterval(12,30) do
        K := GF(p^2);
        s := 20;
        if g eq 3 then
            time sample_eightpoints(K,s);
        elif g eq 4 then
            time sample_tenpoints(K,s);
        elif g eq 5 then
            time sample_twelvepoints(K,s);
        end if;
    end for;
end for;

*/

/*

//Example showing computation of p-ranks for genus g = 3, 4, 5
for g in [3,4,5] do
    for p in PrimesInInterval(12,30) do
        r := 2;
        s := 10^4;
        time p_ranks := CurveDataForP(g, p, r, s); p_ranks;

        non_ord := &+p_ranks[1..g];

        tot := &+p_ranks;

        RealField(20)!(p^r * non_ord/tot);
    end for;
end for;

*/

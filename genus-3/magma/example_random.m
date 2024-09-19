SetSeed(1); //to have consistent and reproducible data
Attach("~/github/Non-ordinary-locus/tools.m");
Attach("~/github/Non-ordinary-locus/genus-3/magma/data_code.m");

s := 10000000;
r := 2;
for p in PrimesInInterval(3,3) do
    "p =", p;
    time p_r := CurveDataForP(p, r, s);
    "p-rank data:", p_r;
    "number of components non-ordinary locus:",RealField(5)!((&+p_r[1..#p_r-1])*(p^r)/&+p_r);
    Write("~/g3-p_rank_fp2_new.m", "[*" cat Sprint(p) cat ", " cat Sprint(r) cat ", " cat Sprint(&+p_r) cat ", " cat Sprint(p_r) cat "*],");
end for;

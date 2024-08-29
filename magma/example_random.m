SetSeed(1); //to have consistent and reproducible data
Attach("~/github/Non-ordinary-locus/magma/data_code.m");

s := 50000000;
r := 2;
for p in PrimesInInterval(499,499) do
    "p =", p;
    time p_r := CurveDataForP(p, r, s);
    "p-rank data:", p_r;
    "number of components non-ordinary locus:",RealField(5)!((&+p_r[1..4])*(p^r)/&+p_r);
    Write("~/p_rank_fp2_new.m", "[*" cat Sprint(p) cat ", " cat Sprint(p_r) cat "*],");
end for;

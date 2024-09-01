SetSeed(1); //to have consistent and reproducible data
Attach("~/github/Non-ordinary-locus/magma/data_code.m");

s := 100000000;
r := 1;
for p in PrimesInInterval(450,500) do
    "p =", p;
    time p_r := CurveDataForP(p, r, s);
    "p-rank data:", p_r;
    "number of components non-ordinary locus:",RealField(5)!((&+p_r[1..4])*(p^r)/&+p_r);
    Write("~/p_rank_fp_new.m", "[*" cat Sprint(p) cat ", " cat Sprint(p_r) cat "*],");
end for;

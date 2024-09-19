SetSeed(1); //to have consistent and reproducible data
Attach("~/github/Non-ordinary-locus/genus-4/magma/data_code.m");

s := 100000000;
r := 2;
for p in PrimesInInterval(280, 300) do
    "p =", p;
    time p_r := CurveDataForP(p, r, s);
    "p-rank data:", p_r;
    "number of components non-ordinary locus:",RealField(5)!((&+p_r[1..4])*(p^r)/&+p_r);
    Write("~/p_rank_fp2_new.m", "[*" cat Sprint(p) cat ", " cat Sprint(r) cat ", " cat Sprint(&+p_r) cat ", " cat Sprint(p_r) cat "*],");
end for;

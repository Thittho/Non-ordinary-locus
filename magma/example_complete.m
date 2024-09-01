Attach("~/github/Non-ordinary-locus/magma/data_code.m");

for p in PrimesInInterval(31,31) do
    res := [0,0,0,0,0];

    for i0 in [30..30] do
        "i0 := ", i0;
        p_r := CurveDataForPComplete(p, 1 : i0 := i0);
        "p-rank data:", p_r;
        res := [res[i]+p_r[i] : i in [1..5]];
    end for;
    res;
    Write("~/p_rank_fp_complete_new.m", "[*" cat Sprint(p) cat ", 1, " cat Sprint(res) cat "*],");
end for;

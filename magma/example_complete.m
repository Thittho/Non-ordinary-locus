Attach("~/github/Non-ordinary-locus/magma/data_code.m");

p := 47;

/*for i0 in GF(p^2) do
    "i0 := ", i0;
    time p_r := CurveDataForPComplete(p, 2 : i0 := i0);
    "p-rank data:", p_r;
    Write("~/p_rank_fp2_complete.m", "[*" cat Sprint(p) cat ", " cat Sprint(p_r) cat "*],");
end for;
*/

for i0 in [43..46] do
    "i0 := ", i0;
    time p_r := CurveDataForPComplete(p, 1 : i0 := i0);
    "p-rank data:", p_r;
    Write("~/p_rank_fp_complete.m", "[*" cat Sprint(p) cat ", " cat Sprint(p_r) cat "*],");
end for;

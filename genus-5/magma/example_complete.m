Attach("~/github/Non-ordinary-locus/genus-5/magma/data_code.m");

r := 1;
for p in PrimesInInterval(31,31) do
    res := [0,0,0,0,0];

    for i0 in [15..19] do
        "i0 := ", i0;
        p_r := CurveDataForPComplete(p, r : i0 := i0);
        "p-rank data:", p_r;
        res := [res[i]+p_r[i] : i in [1..#res]];
    end for;
    res;
    Write("~/p_rank_fp_complete_new.m", "[*" cat Sprint(p) cat ", " cat Sprint(r) cat ", " cat Sprint(&+res) cat ", " cat Sprint(res) cat "*],");
end for;

/*
Attach("~/github/Non-ordinary-locus/magma/smaller_families.m");
res := [];

for p in PrimesInInterval(13,13) do
    pr7 := CurveDataForPComplete(p, 1);
    pr6 := [CurveDataForPCompleteDim6(p, 1 : fam := i) : i in [1..2]];
    pr5 := [CurveDataForPCompleteDim5(p, 1 : fam := i) : i in [1..3]];
    pr4 := [CurveDataForPCompleteDim4(p, 1 : fam := i) : i in [1..4]];
    pr3 := [CurveDataForPCompleteDim3(p, 1 : fam := i) : i in [1..5]];
    pr7, pr6, pr5, pr4, pr3;
    [pr7[i]+pr6[1][i]+pr6[2][i]+pr5[1][i]+pr5[2][i]+pr5[3][i]+pr4[1][i]+pr4[2][i]+pr4[3][i]+pr4[4][i]+pr3[1][i]+pr3[2][i]+pr3[3][i]+pr3[4][i]+pr3[5][i] : i in [1..5]];
    Write("~/p_rank_fp_complete_new.m", Sprint(p) cat ", 1, " cat Sprint([pr7[i]+pr6[1][i]+pr6[2][i]+pr5[1][i]+pr5[2][i]+pr5[3][i]+pr4[1][i]+pr4[2][i]+pr4[3][i]+pr4[4][i]+pr3[1][i]+pr3[2][i]+pr3[3][i]+pr3[4][i]+pr3[5][i] : i in [1..5]]));
    //Append(~res, [pr7[i]+pr6[1][i]+pr6[2][i]+pr5[1][i]+pr5[2][i]+pr5[3][i]+pr4[1][i]+pr4[2][i]+pr4[3][i]+pr4[4][i]+pr3[1][i]+pr3[2][i]+pr3[3][i]+pr3[4][i]+pr3[5][i] : i in [1..5]]);
end for;
res;
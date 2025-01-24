load "galois.m";

g := 4; // g is the genus
r := 3; // r is the power of p
s := 100000000; // s is the number of samples

doc := "~/github/Non-ordinary-locus/data/galois_family/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(5,5) do
    printf "p = %o\n", p;
    time p_ranks := CurveDataForPGalois(g, p, r, s);
    printf "p-rank data = %o\n", p_ranks;
    non_ord := &+p_ranks[1..g];
    tot := &+p_ranks;
    if tot ne 0 then
        RealField(20)!(p^r * non_ord/tot);
        //Write(doc, [p, r, tot] cat p_ranks);
    end if;
end for;

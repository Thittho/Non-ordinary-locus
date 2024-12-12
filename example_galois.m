load "galois.m";

g := 5;
r := 2;
s := 10000000;

doc := "~/github/Non-ordinary-locus/data/galois_family/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(3,100) do
    time p_ranks := CurveDataForPGalois(g, p, r, s);
    "p-rank data:", p_ranks;
    non_ord := &+p_ranks[1..g];
    tot := &+p_ranks;
    if tot ne 0 then
        RealField(20)!(p^r * non_ord/tot);
    end if;
    Write(doc, [p, r, tot] cat p_ranks);
end for;

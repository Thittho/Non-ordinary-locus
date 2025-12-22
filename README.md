# Irreducibility of p-Rank Strata

In this repository, you will find all of the code and data related to the paper "Heuristics for (ir)reducibility of p-rank strata of the moduli space of hyperelliptic curves," which is available at [arXiv:2506.06457](https://arxiv.org/abs/2506.06457).  Information about usage of the code can be found below:

## Family Method

## Galois Type Method

The following information concerns the so-called Galois Type Method, described in Section 3.2 of the paper.  


```cpp

load "galois.m";

g := 5; // g is the genus
r := 3; // r is the power of p
s := 100000000; // s is the number of samples

doc := "~/github/Non-ordinary-locus/data/galois_family/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(97,97) do
    printf "p = %o\n", p;
    time p_ranks := CurveDataForPGalois(g, p, r, s);
    printf "p-rank data = %o\n", p_ranks;
    non_ord := &+p_ranks[1..g];
    tot := &+p_ranks;
    if tot ne 0 then
        "p = ", p;
        RealField(20)!(p^r * non_ord/tot);
        Write(doc, [p, r, tot] cat p_ranks);
    end if;
end for;
```

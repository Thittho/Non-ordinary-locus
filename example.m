Attach("~/github/Non-ordinary-locus/tools.m");

r := 3;
g := 5;
s := 100000000; // put s to 0 if you want the whole family
doc := "~/github/Non-ordinary-locus/data/7_dim_family/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(67,67) do
  if (2*g+1) mod p ne 0 then
    p;
    res := CurveDataForP(g, p, r : s := s);
    "p-rank data:", res;
    Write(doc, [p, r, &+res] cat res);
  end if;
end for;
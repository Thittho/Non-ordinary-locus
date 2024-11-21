Attach("~/github/Non-ordinary-locus/tools.m");

r := 1;
g := 3;
s := 10000000; // put s to 0 if you want the whole family
doc := "~/github/Non-ordinary-locus/data/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(900,1000) do
  if (2*g+1) mod p ne 0 then
    p;
    res := CurveDataForP(g, p, r : s := s);
    "p-rank data:", res;
    Write(doc, [p, r, &+res] cat res);
  end if;
end for;

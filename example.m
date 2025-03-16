Attach("~/github/Non-ordinary-locus/tools.m");

r := 2;
g := 20;
s := 1000000; // put s to 0 if you want the whole family
doc := "~/p_rank_fp" cat Sprint(r) cat "g" cat Sprint(g) cat ".m";

for p in PrimesInInterval(3,100) do
  if (2*g+1) mod p ne 0 then
    p;
    res := CurveDataForP(g, p, r : s := s);
    "p-rank data:", res;
    Write(doc, [p, r, &+res] cat res);
  end if;
end for;

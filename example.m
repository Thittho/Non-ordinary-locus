Attach("~/github/Non-ordinary-locus/tools.m");

r := 3;
g := 5;
s := 100000000; // put s to 0 if you want the whole family
doc := "~/github/Non-ordinary-locus/data/genus-" cat Sprint(g) cat "/p_rank_fp" cat Sprint(r) cat ".m";

for p in PrimesInInterval(70,80) do
  if (2*g+1) mod p ne 0 then
    p;
    res := CurveDataForP(g, p, r : s := s);
    "p-rank data:", res;
    Write(doc, [p, r, &+res] cat res);
  end if;
end for;


P1P1 := ProductProjectiveSpace(QQ, [1,1]);
R<x,y,u,v> := CoordinateRing(P1P1);
f := x^3*u^3+2*x^2*y*u*v^2+y^3*u^2*v+3*x^2*y*v^3;
C := Curve(P1P1, f);
BasisOfHolomorphicDifferentials(C);


R<x,y> := PolynomialRing(Rationals(), 2);
C := Curve(Spec(R), y^5-x*(x-1)*(x-2));
E, Q := Explode(Equations(Image(CanonicalEmbedding(C))));
S<x,y,z,t> := PolynomialRing(Rationals(), 4);
E := S!E;
Q := S!Q;
Q;
E;

lambda := 2;
Q0 := x*z-y^2;
E0 := t^3-1/3*(lambda^2-lambda+1)*t*z^2-1/27*(1+lambda)*(lambda-2)*(lambda+1/2)*z^3-x^2*y;
InvariantsGenus4Curves(Q0,E0);

T<s, t, w> := PolynomialRing(BaseRing(Parent(E0)), [1,1,2]);
f_weighted := Evaluate(E, [w, s^2, s*t, t^2]);
		vprint Genus4 : "Computing normal form of the sextic...";
		require MonomialCoefficient(f_weighted, w^3) ne 0: "The curve is not smooth";

		// we put the curve in normal form
		alpha := MonomialCoefficient(f_weighted, w^3);
		//"alpha", alpha;
		f_weighted /:= alpha;
		f_weighted := Evaluate(f_weighted, [s, t, w-ExactQuotient(Terms(f_weighted, w)[3], 3*w^2)]);

L<x,y,l> := PolynomialRing(Rationals(), 3);
f := y^5-x*(x-1)*(x-l);
 
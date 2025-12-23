# Irreducibility of p-Rank Strata

In this repository, you will find all of the code and data related to the paper "Heuristics for (ir)reducibility of p-rank strata of the moduli space of hyperelliptic curves," which is available at [arXiv:2506.06457](https://arxiv.org/abs/2506.06457).  Information about usage of the code can be found below.  The first two sections detail how to build a sample of hyperelliptic curves, and the third section describes how to compute $p$-ranks of these curves.  

## Family Method

Additional documentation coming shortly...

## Galois Type Method

The following information concerns the **Galois Type Method**, described in Section 3.2 of the paper.  The program [galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/galois.m) loads one primary command, `CurveDataForPGalois`, which provides a method for building a sample of hyperellptic curves over $\mathbb{F}_q$.

Let $p$ be an odd prime, $r \geq 1$ an integer, and $q = p^r$.  Then, for $g \in \{3, 4, 5\}$, a prime power $q > 2g + 1$, and an integer $s \geq 1$, the command `CurveDataForPGalois(g, p, r, s)` returns a list $L$ of size $s$ of polynomials $f(x) \in \mathbb{F}_q[x]$ that split completely over $\mathbb{F}_q$, and which represent hyperelliptic curves $y^2 = f(x)$ of genus $g$ that are unique up to geometric isomorphism (see Algorithm 3.7).

### Input
* `g` - the genus (an integer 3, 4, or 5)
* `p` - a rational prime
* `r` - the power of p
* `s` - the number of curves desired in the sample

### Output
* `L` - a list of polynomials $f(x) \in \mathbb{F}_q[x]$ that split completely over $\mathbb{F}_q$

### Example
A small example is given below.  A longer example can be found in the file [example_galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/example_galois.m).
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

**Note**: In order to load [galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/galois.m) successfully, our helper program [tools.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/tools.m), as well as Everett Howe's [Hyperelliptic3.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic3.magma) and [Hyperelliptic2.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic2.magma) must be accessible. 

## Computing p-Ranks

Additional documentation coming shortly...

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
A small example is given below.  This example was run on a ThinkPad X1 Carbon (6th Gen) equipped with an [IntelCore i5-8350U Processor](https://www.intel.com/content/www/us/en/products/sku/124969/intel-core-i58350u-processor-6m-cache-up-to-3-60-ghz/specifications.html) and 16GB of RAM.

A longer example can be found in the file [example_galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/example_galois.m).
```cpp
> load "galois.m";
Loading "galois.m"
Loading "Hyperelliptic3.magma"
Loading "Hyperelliptic2.magma"
>
> CurveDataForPGalois(4, 13, 2, 1000000);
Creating list of curves...
Time: 587.430
Computing p-ranks...
Time: 60.150
[ 0, 0, 48, 6063, 993889 ]
```

**Note**: In order to load [galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/galois.m) successfully, our helper program [tools.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/tools.m) as well as Everett Howe's [Hyperelliptic3.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic3.magma) and [Hyperelliptic2.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic2.magma) must be accessible. 

### Data
The data collected using this method can be found in the directory [data/galois_family](https://github.com/Thittho/Non-ordinary-locus/tree/main/data/galois_family).  The data is broken up by genus. 

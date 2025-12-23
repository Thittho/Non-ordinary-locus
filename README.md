# Irreducibility of p-Rank Strata

In this repository, you will find all of the code and data related to the paper "Heuristics for (ir)reducibility of p-rank strata of the moduli space of hyperelliptic curves," which is available at [arXiv:2506.06457](https://arxiv.org/abs/2506.06457).  Information about usage of the code can be found below.  The first two sections detail how to build a sample of hyperelliptic curves, and the third section describes how to compute $p$-ranks of these curves.  

## Family Method

The following information concerns the **Family Method**, described in Section 3.1 of the paper.  The program [tools.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/tools.m) loads several commands.  The primary command, `CurveDataForP`, is used to create a sample of hyperelliptic curves over $\mathbb{F}_q$ in a specific family and compute their $p$-ranks.

Let $p$ be an odd prime, $r \geq 1$ an integer, and $q = p^r$.  Then, the command `CurveDataForP(g, p, r, s)` first builds a list $L$ of size $s$ of polynomials $f(x) \in \mathbb{F}_q[x]$ that are in the family, and which represent hyperelliptic curves $y^2 = f(x)$ of genus $g$ that are unique up to geometric isomorphism, and then computes the $p$-ranks of these curves.  The command returns a list containing the number of curves in the sample of each $p$-rank.  

### Input
* `g` - the genus
* `p` - a rational prime
* `r` - the power of $p$
* `s` - the number of curves desired in the sample; set `s := 0` if you want the whole family

### Output
* `res` - a list of the number of curves in the sample of each $p$-rank
  - the list is organized as [# of $p$-rank 0 curves, # of $p$-rank 1 curves, ...]

### Example
A small example is given below.  This example was run on a ThinkPad X1 Carbon (6th Gen) equipped with an IntelCore i5-8350U processor and 16GB of RAM.

```cpp
> load "tools.m";
Loading "tools.m"
>
> CurveDataForP(4, 13, 2 : s := 1000000);
Creating curves...
Time: 5.500
Checking which curves are smooth...
Time: 2.290
Removing isomorphic curves...
Time: 74.940
Computing p-ranks...
Time: 59.720
[ 0, 0, 35, 5959, 994006 ]
```

A longer example can be found in the file [example.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/example.m).

### Data
The data collected using this method for genus $g = 3, 4, 5$ can be found in the directory [data/7_dim_family](https://github.com/Thittho/Non-ordinary-locus/tree/main/data/7_dim_family).  The data for genus $6 \leq g \leq 20$ can be found in [data/higher_genus_data](https://github.com/Thittho/Non-ordinary-locus/tree/main/data/higher_genus_data).


## Galois Type Method

The following information concerns the **Galois Type Method**, described in Section 3.2 of the paper.  The program [galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/galois.m) loads one primary command, `CurveDataForPGalois`, which provides a method for building a sample of hyperelliptic curves that split completely over $\mathbb{F}_q$ and computing their $p$-ranks.

Let $p$ be an odd prime, $r \geq 1$ an integer, and $q = p^r$.  Then, for $g \in \{3, 4, 5\}$, a prime power $q > 2g + 1$, and an integer $s \geq 1$, the command `CurveDataForPGalois(g, p, r, s)` first builds a list $L$ of size $s$ of polynomials $f(x) \in \mathbb{F}_q[x]$ that split completely over $\mathbb{F}_q$, and which represent hyperelliptic curves $y^2 = f(x)$ of genus $g$ that are unique up to geometric isomorphism (see Algorithm 3.7), and then computes the $p$-ranks of these curves.  The command returns a list containing the number of curves in the sample of each $p$-rank.  

### Input
* `g` - the genus (an integer 3, 4, or 5)
* `p` - a rational prime
* `r` - the power of $p$
* `s` - the number of curves desired in the sample

### Output
* `res` - a list of the number of curves in the sample of each $p$-rank
  - the list is organized as [# of $p$-rank 0 curves, # of $p$-rank 1 curves, ...]

### Exampl
A small example is given below.  This example was run on a ThinkPad X1 Carbon (6th Gen) equipped with an IntelCore i5-8350U processor and 16GB of RAM.

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

A longer example can be found in the file [example_galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/example_galois.m).

**Note**: In order to load [galois.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/galois.m) successfully, our helper program [tools.m](https://github.com/Thittho/Non-ordinary-locus/blob/main/tools.m) as well as Everett Howe's [Hyperelliptic3.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic3.magma) and [Hyperelliptic2.magma](https://github.com/everetthowe/hyperelliptic/blob/main/Hyperelliptic2.magma) must be accessible. 

### Data
The data collected using this method can be found in the directory [data/galois_family](https://github.com/Thittho/Non-ordinary-locus/tree/main/data/galois_family).  The data is broken up by genus. 

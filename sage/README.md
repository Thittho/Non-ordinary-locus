## Sage Code Examples

### Full Family Example

The command `CurveDataFull(p, r)` returns a count of the number of curves with various $p$-ranks for the full family of curves over $\mathbb{F}_q$, where $q = p^r$.

INPUT:
* `p` - a rational prime number
* `r` - a positive, rational integer

```python
sage: attach("PATH/Curve_Data.sage")
sage: CurveDataFull(23, 1)
(23, 1, 141599590, 941, 11198, 271326, 5927070, 135389055)
```
The data is output in the following format:

(p, r, family size, # p-rank 0 curves, # p-rank 1 curves, # p-rank 2 curves, # p-rank 3 curves, # p-rank 4 curves)

### Sample Example

The command `CurveDataSample(p, r, s, t0, t1)` returns a count of the number of curves with various $p$-ranks for a sample of curves over $\mathbb{F}_q$.  Again, $q = p^r$ for $p$ a rational prime, but this command also allows the user to specify the sample size $s$.

Further, this function takes two optional arguments, allowing the user to specify coefficients in the model of the hyperelliptic curve family.  The user can set the constant coefficient $t_0$, or the constant and linear coefficients $t_0$ and $t_1$.  If the optional arguments are not specified, then these coefficients will be randomly chosen from $\mathbb{F}_q$. 

INPUT:
* `p` - a rational prime number
* `r` - a positive, rational integer
* `s` - sample size
* `t0` - (optional) fixed constant coefficient
* `t1` - (optional) fixed linear coefficient; requires `t0` to also be set

```python
sage: attach("PATH/Curve_Data.sage")
sage: CurveDataSample(5, 3, 50000000)
(5, 3, 50000000, 0, 27, 3165, 396244, 49600564)
```

The data is output in the following format:

(p, r, sample size, # p-rank 0 curves, # p-rank 1 curves, # p-rank 2 curves, # p-rank 3 curves, # p-rank 4 curves)

Below is an example of how these optional arguments might be used on machines with limited memory capacity.  Instead of having to store 4.9 million tuples of coefficients, the computer only has to store 100,000 tuples at a time.  As $t_0$ ranges over $\mathbb{F}_q$, the results are added together.  

```python
sage: attach("PATH/Curve_Data.sage")
sage: Fq = GF(7**2)
sage: output = [7, 2, 0, 0, 0, 0, 0, 0] 
sage: for i in Fq:
....:    res_i = CurveDataSample(7, 2, 100000, i)
....:    for j in range(2, res_i):
....:        output[j] += res_i[j]
....:
sage: output
[7, 2, 4900000, 0, 41, 1958, 97851, 4800150]
```





# Introduction

This document is just a place for me to write out the math as I work it out, so that I'm not constantly erasing / writing new papers.

Some notes:
- If/when this needs to get properly compiled, will need to load some LaTeX packages, and probably define some new commands.
  - `\usepackage{amsmath}`
  - Sign function: `\DeclareMathOperator{\sign}{sign}`, currently `\text{sgn}`
  - 

## Main formulas



## Chebyshev: Formulas and derivation

### Main equations

Chebyshev:

$$
P( |X-\mu| \geq k\sigma ) \leq \dfrac{1}{k^2}
$$

Finite Chebyshev (double-check this):

$$
P( |x-\bar{x}| \geq ks ) \leq \dfrac{1}{n-1} \left \lfloor{\dfrac{n+1}{n}\left(\dfrac{n-1}{k^2} + 1\right)}\right \rfloor 
$$


Cantelli

$$
P( X - \mu \geq k\sigma ) \leq \dfrac{1}{1+k^2}
$$


Finite Cantelli





### Derivation of finite Cantelli


Let $Y = X - \bar{X}$, then:

$$
\begin{aligned}
\mu_{y} &= \mathbb{E}[Y]  = 0 \\
\sigma_{y}^{2} &= \mathbb{V}[Y] = \dfrac{n-1}{n}\sigma^{2}_{x}
\end{aligned}
$$

Therefore $\sigma_{x} = \sqrt{\dfrac{n}{n-1}}\sigma_{y}$

**Part 1**

Apply regular Cantelli to $Y$:

$$
P( Y - \mu_{y} \geq k_{y} \sigma_{y} ) \leq \dfrac{1}{1+k_{y}^2}
$$

Then swap to the original variance for $X$:

$$
P( X - \bar{X} \geq k_{y} \sqrt{\dfrac{n}{n-1}}\sigma_{x} ) \leq \dfrac{1}{1+k_{y}^2}
$$

And absorb the fraction into $k$, which puts it into denominator. Stats from setting $k_{x} = k_{y}\sqrt{n/(n-1)}$.

$$
P( X - \bar{X} \geq k_{x} \sigma_{x} ) \leq \dfrac{1}{1+ \dfrac{n}{n-1}k_{x}^2}
$$

Comments:
- This bound is **smaller** because the $n/(n-1)$ makes the denominator **larger**.
- This solves one portion, getting $\bar{X}$ instead of $\mu_{x}$. However, it still uses $\sigma$ instead of $s$.


**Part 2**

Next, consider

$$
\begin{aligned}
P( Y \geq ks ) &= \mathbb{E}[ \mathbb{I}( Y \geq ks ) ] \\
&\approx \dfrac{1}{n} \sum_{i} \mathbb{I}( Y_{i} \geq ks )
\end{aligned}
$$

Focusing on the indicator. The third line follows because the RHS of the indicator function is strictly positive. See the table following.

$$
\begin{aligned}
\mathbb{I}( Y_{i} \geq ks ) &= \mathbb{I}( Y_{i} + kcs \geq ks + (ksc) ) \\
&= \mathbb{I}( Y_{i} + kcs \geq ks(1+c) ) \\
&\leq \mathbb{I}( (Y_{i} + kcs)^{2} \geq t^{2}s^{2}(1+c)^{2} ) \\
&\leq \dfrac{(Y_{i} + kcs)^{2}}{k^{2}s^{2}(1+c)^{2}}
\end{aligned}
$$

    a   |a|>|b|   a>=b  a^2 >= b^2 
    ------------------------------
    +      1        1       1
    +      0        0       0
    -      0        0       0
    -      1        0       1

Returning to the previous, note that since $Y = X-\bar{X}$, we have $\sum Y = 0$ and $\sum Y_{i}^{2} = (n-1)s$ (and note $s_{x} = s_{y}$ since it is just scale-shifted).

$$
\begin{aligned}
\dfrac{1}{n} \sum_{i} \mathbb{I}( Y_{i} \geq ks ) &\leq 
\dfrac{1}{n}\sum_{i} \dfrac{(Y_{i} + kcs)^{2}}{k^{2}s^{2}(1+c)^{2}} \\
&= \dfrac{ (\sum_{i} Y_{i}^{2}) + (2kcs \sum_{i}Y_{i}) + (nk^{2}s^{2}c^{2}) }{ n k^{2}s^{2}(1+c)^{2}} \\
&= \dfrac{ (n-1)s^{2} + nk^{2}s^{2}c^{2} }{ nk^{2}s^{2}(1+c)^{2}} \\
&= \dfrac{ (n-1) + nk^{2}c^{2} }{ nk^{2}(1+c)^{2}} \\
\end{aligned}
$$

So this is the *upper bound* on $P( Y \geq ks )$. To get the *smallest upper bound*, we need to minimize this with respect to $c$. This is a calculus problem, using the quotient rule:

If $H(c) = f(c) / g(c)$, then:

$$
H'(c) = \dfrac{f'(c)g(c) - f(c)g'(c)}{g(c)^2}
$$

This is a lot of algebra, but comes out to:

$$
c = \dfrac{n-1}{nk^{2}}
$$

Plugging this into the above upper bound and simplifying obtains:

$$
P( Y \geq ks ) \leq \dfrac{1}{1 + \dfrac{n}{n-1}k^{2}}
$$

So we have:

$$
P( X - \bar{X} \geq ks_{x} ) \leq \dfrac{1}{1 + \dfrac{n}{n-1}k^{2}}
$$


To obtain the value of $k$ needed to achieve a certain coverage, we set the RHS equal to $\alpha$ and obtain:

$$
k = \sqrt{ \dfrac{1-\alpha}{\alpha} \dfrac{n-1}{n} }
$$










## Example 1.7.1, Page 43

The contaminated normal PDF is given by:

$$
f_{\varepsilon}(x) = (1-\varepsilon)\phi(x) + \varepsilon \delta^{-1} \phi(\delta^{-1}x)
$$

The variance for this comes out to $\sigma^{2}_{f} = 1 + \varepsilon(\delta^{2} - 1)$. This result is given on page 27, and makes sense to work out. Then we have:

$$
\int f_{\varepsilon}(x)^{2} = 
\left[(1-\varepsilon)^{2}\phi^{2}(x) \right]
+ 
\left[ \varepsilon^{2} \delta^{-2} \phi^{2}(\delta^{-1}x) \right]
+
\left[2\varepsilon(1-\varepsilon)\delta^{-1} \phi(x) \phi(\delta^{-1}x)\right]
$$

So there are three cases, indicated by the brackets. I'm going to solve this in general for the middle case (which then solves the first case by setting $\delta=1$).

**Case 2**

I puleld the constants, $\varepsilon^{2}$ and $\delta^{2}$ out front first, then squared the normal PDF, hence tehre is no $1/2$ in the exponent. I'm taking one of the normalizing constants out front as well, since the goal is to turn the integral into a normal PDF.

$$
\begin{aligned}
\int \varepsilon^{2} \delta^{-2} \phi^{2}(\delta^{-1}x) &= 
\dfrac{\varepsilon^{2}}{\delta^{2}\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{\sigma^{2}\delta^{2}}} dx
\end{aligned}
$$

Using substitution, I want the exponent to look like $- y^{2} / (2\sigma^{2})$, so I'm setting $y^{2}/2 = x^{2}/\delta^{2}.$ This means the following:

$$
\begin{aligned}
y &= \dfrac{\sqrt{2}}{\delta}x \\
y^{2} &= \dfrac{2}{\delta^{2}}x^{2} \\
dy &= dx  \dfrac{\sqrt{2}}{\delta} \\
dx &= dy  \dfrac{\delta}{\sqrt{2}}
\end{aligned}
$$

Since the bounds of integration are from $-\infty$ to $\infty$, they don't change here. So the integral above turns into:

$$
\begin{aligned}
\int \varepsilon^{2} \delta^{-2} \phi^{2}(\delta^{-1}x) 
&= \dfrac{\varepsilon^{2}}{\delta^{2}\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{\sigma^{2}\delta^{2}}} dx \\
&= \dfrac{\varepsilon^{2}}{\delta^{2}\sqrt{2\pi\sigma^{2}}} 
\int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{y^{2}}{2\sigma^{2}}} dy \dfrac{\delta}{\sqrt{2}} \\
&= \dfrac{\varepsilon^{2}}{2\delta\sqrt{\pi\sigma^{2}}} 
\int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{y^{2}}{2\sigma^{2}}} dy
\end{aligned}
$$

The integral now evaluates to 1, so all we're left with is the constant out front. When $\sigma=1$, that part disappears, so we have:

$$
\dfrac{\varepsilon^{2}}{2\delta\sqrt{\pi}} 
$$

Setting $\delta=3$, which makes the contaminated part of the variance equal to 9, gets the middle term. Replacing $\varepsilon^2$ with $(1-\varepsilon)^2$ and setting $\delta=1$ gets the first component of the original integral. So the first two terms are 

$$
\dfrac{(1-\varepsilon)^2}{2\sqrt{\pi}} 
+
\dfrac{\varepsilon^{2}}{6\sqrt{\pi}} 
$$

which checks out with Example 1.7.1. Where I'm having an issue is in the third term.

**Case 3**


$$
\begin{aligned}
\int 2\varepsilon(1-\varepsilon) \delta^{-1} \phi^{2}(x) \phi^{2}(\delta^{-1}x) 
&= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{2\sigma^{2}} -\dfrac{x^{2}}{2\sigma^{2}\delta^{2}} } dx \\
(*) &= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{2\sigma^{2}}\left( 1 + \delta^{-2} \right) } dx \\
\end{aligned}
$$

This change-of-variable is slightly easier, since the $1/2$ is already there. 

$$
\begin{aligned}
y &= x \sqrt{1 + \delta^{-2}} \\
y^{2} &= x^2 \left( 1 + \delta^{-2} \right) \\
dy &= dx \sqrt{1 + \delta^{-2}} \\
dx &= dy \dfrac{1}{\sqrt{1 + \delta^{-2}}}
\end{aligned}
$$

And again the bounds stay the same. So the integral for Case 3 comes out to:



$$
\begin{aligned}
\int 2\varepsilon(1-\varepsilon) \delta^{-1} \phi^{2}(x) \phi^{2}(\delta^{-1}x) 
&= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{2\sigma^{2}} -\dfrac{x^{2}}{2\sigma^{2}\delta^{2}} } dx \\
(*) &= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{x^{2}}{2\sigma^{2}}\left( 1 + \delta^{-2} \right) } dx \\
&= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{y^{2}}{2\sigma^{2}} } dy \dfrac{1}{\sqrt{1 + \delta^{-2}}} \\
&= \dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \dfrac{1}{\sqrt{1 + \delta^{-2}}} \int \dfrac{1}{\sqrt{2\pi\sigma^{2}}} e^{-\dfrac{y^{2}}{2\sigma^{2}} } dy \\
\end{aligned}
$$

Again the integral comes out to $1$, so all that's left is the constant out front. Simplifying this:

$$
\begin{aligned}
\dfrac{2\varepsilon(1-\varepsilon)}{\delta\sqrt{2\pi\sigma^{2}}} \dfrac{1}{\sqrt{1 + \delta^{-2}}} 
&= \dfrac{2\varepsilon(1-\varepsilon)}{\sqrt{2\pi\sigma^{2}}\sqrt{\delta^2 + 1}} \\
&= \dfrac{ \sqrt{2}\sqrt{2} \varepsilon(1-\varepsilon)}{\sqrt{2\pi\sigma^{2}}\sqrt{\delta^2 + 1}} \\
&= \dfrac{ \sqrt{2} \varepsilon(1-\varepsilon)}{\sqrt{\pi\sigma^{2}}\sqrt{\delta^2 + 1}} \\
&= \dfrac{ \sqrt{2} \varepsilon(1-\varepsilon)}{\sqrt{\pi\sigma^{2}}\sqrt{\delta^2 + 1}} \dfrac{1/\sqrt{2}}{1/\sqrt{2}}\\
&= \dfrac{\varepsilon(1-\varepsilon)}{\sqrt{\pi\sigma^{2}}\sqrt{(\delta^2 + 1)/2}} \\
\end{aligned}
$$







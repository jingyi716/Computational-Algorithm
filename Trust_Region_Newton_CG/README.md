# Trust-Region-Newton-CG
# Trust-Region-Newton-CG
# Trust-Region-Newton-CG
Test problem: Quadratic function:\\
    Consider the convex quadratic function 
    $$
    \text{minimize}\quad f(x)=\frac{1}{2}x^TQx+q^Tx
    $$
    where $q=\mathbb{R}^n$ and $Q\in \mathbb{R}^{n\times n} $ be randomly generated convex quadratic function. We choose $x^{(0)}=20*\text{rand}(n,1)-10$ at the initial value, where n=10 and 1000, depending on the dimension of Q/\\
    We change the dimensions of $Q$ on n=10 and n=1000 and the conditional number of $Q$ ($\kappa=10$ and to $\kappa=1000$) to see how each method performs.

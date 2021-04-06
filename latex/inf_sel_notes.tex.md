


# Some questions

1. Variogram/Covariogram

* Q1: The formulas in the paper are for the general case $U\in\mathbb{R}^{d}$ or for the case $U\in\mathbb{R}^{2}$?

R1: For the definition of $G$ and $C$, this is the general case: $U= \mathbb{R}^{d}$

2. **Equation (5)**

$Y:\Omega\to\left(U\to\mathcal{Y}\right)$ Gaussian random process

$X:\Omega\to U$

$n,m\in\mathcal{N}$, $x\in U^{n}, x^{'}\in U^{m}$

* Q2: $n\ne m$?

R2: not necessarily.

Expected value of signal $\mu:U\to\mathcal{Y}$, $x\to E\left[Y\left[x\right]\right]$ and covariance matrix between points $\textbf{x}_{1},\dots,\textbf{x}_{n}$ and $\textbf{x}^{'}_{1},\dots,\textbf{x}^{'}_{m}$ $\Sigma_{\textbf{x},\textbf{x}^{'},\theta}=Cov\left[Y\left[\textbf{x}\right],\left[\textbf{x}^{'}\right]\right]$.

In the paper we have 
$$f_{Y[\textbf{x}]}\left(\textbf{y}\right)=\frac{1}{2\pi^{n/2}}\exp{\frac{-(\textbf{y}-\mu(\textbf{x}))\Sigma_{\textbf{x},\textbf{x}}^{-1}(\textbf{y}-\mu(\textbf{x}))}{2}}.
$$

* Q3: Most likely I'm mistaken, but I was wondering if (i) in the first fraction the denominator has also the term $|\Sigma|^{n/2}$ where $|\Sigma|$ is the determinant of $\Sigma$, and (ii) the first $(\textbf{y}-\mu(\textbf{x}))$ is transposed, so that we have
$$f_{Y[\textbf{x}]}\left(\textbf{y}\right)=\frac{1}{2\pi^{n/2}|\Sigma|^{n/2}}\exp{\frac{-(\textbf{y}-\mu(\textbf{x}))^{T}\Sigma_{\textbf{x},\textbf{x}}^{-1}(\textbf{y}-\mu(\textbf{x}))}{2}}.
$$
R3: You are right.


3. **Design/sample**

$D$ is the random process for the design.

* Q4: D has domain equal to the power set of $U^{n}$, i.e. $P\left(U^{n}\right)$?

R4: we have made the assumption of a fixed size design, so in this case, $D(\omega)$ has domain equal to  the power set of $U^{n}$, i.e. $P\left(U^{n}\right)$.
$D$ is a random variable: $D:\Omega \to ${the set of all probabilities on $U^{n}$}, $\omega\mapsto D(\omega)$


* Q5: the codomain is the set of probability distributions on $U^{n}$, right? Can we indicate that with $\mathcal{P}_{U^{n}}$?
R5: the codomain of $D$ is the set of all probabilities distributions on $U^{n}$
the codomain of $D(\omega)$ is $[0,1]$

* Q6: If Q4 and Q5 are right, then we should have $D:P\left(U^{n}\right)\to \mathcal{P}_{U^{n}}$

R6: $D(\omega):\mathcal{P}_{U^{n}}\to [0,1]$

4. **Exchangeability condition**

* Q7: Why do we use it? (I guess to have property 2.3)
R7: this is a way to ignore the order in which elements are selected. Adaptive sampling is a counter example

* Q8: It is verified in our case? I mean, it's just for SRS or not?
R8: Yes, we limit ourselves to these cases. The selection proportional to size, this is also the case.

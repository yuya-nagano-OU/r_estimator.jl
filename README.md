# tensor-to-schalar ratio estimator

Here, we are implementing the likelihood described in Section 5.3.3 of [LiteBIRD PTEP](https://arxiv.org/abs/2202.02773).

This code estimates the tensor-scalar ratio r of the observed $C_{\ell}^{BB}$ using the theoretical $C_{\ell}^{BB}$ and the observed $C_{\ell}^{BB}$ of the universe.

The likelihood function is given by

$$
\log{L\left( r \right)} = \sum_{\ell=\ell_\rm{min}}^{\ell_\rm{max}} \log{P_{\ell}\left( r \right)}
$$

Here, $\ell_\rm{min}$ and $\ell_\rm{max}$ are assumed to be used up to scales particularly important for the tensor mode.


And for $\log{P_{\ell}\left( r \right)}$,

$$
\log{P_{\ell}\left( r \right)} = -f_\rm{sky} \frac{2\ell+1}{2}\left[ \right]
$$

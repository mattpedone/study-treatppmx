# info_simulation_sensitivity

These are the results of the simulation studies for evaluating sensitivity with respect to hyper parameter choice.

The simulation study is designed to evaluate different configurations of parameters $(\kappa, \sigma)$, $(\Sigma, S_0)$ and $\sigma^{2}_{0}$.

Data are generated using the file `R/genscen.R` and saved in `data/`. Fore details on scenarios see `data/info-data.md`.

***

Considering that the parameters of the NGG measure are $NGG(\beta\sigma, \sigma)$, where $\kappa = \beta\times\sigma$, for a better understanding of the behavior of NGG for varying parameters, we will consider the parameter $\beta$.
The parameters are set at the following default values: $\beta = 40$, $\sigma = 0.25$, $\Sigma = 10\mathbf{1}$, $S_0 = 0.1\mathbf{1}$, $\sigma^{2}_{0} = 1$, where $\mathbf{1}$ is $3\times 3$ identity matrix.
Keeping all other parameters fixed, the pairs $(\beta, \sigma)$ and $(\Sigma, S_0)$ and the scalar $\sigma^{2}_{0}$ are evaluated over the following grids:

* $\beta=\{1, 10, 40\}$

* $\sigma=\{0.01, 0.05, 0.25\}$

* $diag(\Sigma)=\{1, 10, 50\}$

* $diag(S_0)=\{0.1, 1, 10\}$

* $\sigma^{2}_{0}=\{1, 2, 10\}$

In particular the values for the grids of $\beta$ and $\sigma$ are fixed such that, for $n=75$ observations, the prior expected number of cluster is:

|$\sigma$/$\beta$ | 1 | 10 | 40 |
|:-:|:-:|:-:|:-:|
| **0.01** | 1.10 | 1.53 | 2.84 |
| **0.05** | 1.54 | 3.63 | 8.63 |
| **0.25** | 5.16 | 14.67 | 28.99 |

See `R/sens-scripts/ev_nlu_ngg.R` to compute prior expected number of cluster.

***

In order to obtain the performances in the 21 scenarios six scripts have been employed:

* `kappasigma_sd.R` for eveluating parameters $\kappa, \sigma$ in scenarios 1, 2
* `kappasigma_lgg.R` for eveluating parameters $\kappa, \sigma$ in LGG data
* `SigmaS0_sd.R` for eveluating parameters $\Sigma, S_0$ in scenarios 1, 2
* `SigmaS0_lgg.R` for eveluating parameters $\Sigma, S_0$ in LGG data
* `sigma20_sd.R` for eveluating parameter $\sigma_{0}^{2}$ in scenarios 1, 2
* `sigma20_lgg.R` for eveluating parameter $\sigma_{0}^{2}$ in LGG data

All the scripts are stored in `R/sens-scripts`.

For each scenario the sensitivity of the parameters has been evluated replicating for $30$ times the LOOCV procedure. Each *replica* consists in 50k MCMC iterations, 20k are discarded due to burn-in and thinning is 1 out of 10.

***

sessioninfo()

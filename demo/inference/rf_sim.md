Prediction of random fields
================

## Set-up

  - Population: simulated Gaussian random field with \(\mu=10\) and
    \(C\left(h\right)=\exp{-h^{2}}\), on a grid \(20\times 20\).
  - Variable of interest: *y*.
  - Sampling designs: SRS, PPS % *z*, PWD spread, PWD cluster.
  - Sample size: 150.
  - Estimation method: *reml*.
  - Number of replications: 100.

![](fig/grf.png)

![](fig/rf_sim0.png)

## Set-up

  - Population: simulated random field with \(\mu=10\) and
    \(C\left(h\right)=\left(1-1.5h+0.5h^{3}\right)1_{\left[0,1\right]\),
    on a grid \(20\times 20\).
  - Variable of interest: *y*.
  - Sampling designs: SRS, PPS % *z*, PWD spread, PWD cluster.
  - Sample size: 150.
  - Estimation method: *reml*.
  - Number of replications: 100.

![](fig/sph.png)

![](fig/rf_sim1.png)

## Set-up

  - Population: simulated random field with \(\mu=10\) and Matern
    covariance, on a grid \(20\times 20\).
  - Variable of interest: *y*.
  - Sampling designs: SRS, PPS % *z*, PWD spread, PWD cluster.
  - Sample size: 150.
  - Estimation method: *reml*.
  - Number of replications: 100.

![](fig/mat.png)

![](fig/rf_sim2.png)

## Set-up

  - Population: simulated random field with \(\mu=10\), scale = 0.10 and
    Gaussian covariance on a grid \(20\times 20\).
  - Variable of interest: *y*.
  - Sampling designs: SRS, PPS % *z*, PWD spread, PWD cluster.
  - Sample size: 150.
  - Estimation method: *reml*.
  - Number of replications: 100.

![](fig/sim3.png)

![](fig/rf_sim3.png)

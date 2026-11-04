# Five-Minute Water-Level Forecasts at the Oosterschelde Storm-Surge Barrier

Regression models that forecast the water level at Roompot Buiten five minutes ahead, just before high tide, to support the decision to close the barrier.

**TU/e · Regression Models and Applications (2MBS50) · April 2025 - July 2025 · Grade 9.0/10**

<!-- REPORT_LINK -->

## Problem

Rijkswaterstaat closes the Oosterschelde storm-surge barrier when the water level at Roompot Buiten (RPBU) reaches 3 m. If the level gets there without a decision, the barrier closes automatically and the emergency team loses control, while closing it unnecessarily harms the environment. The decision team therefore needs reliable forecasts a few minutes ahead. The data are water levels every 10 seconds at RPBU and every 10 minutes at offshore stations.

## Approach

- **Harmonic linear regression:** a linear trend plus sine and cosine terms at the three harbour-oscillation periods (545 s, 205 s and 85 s), fitted on windows of 15–120 minutes before high tide, with residual, Q–Q and Durbin–Watson diagnostics.
- **Quadratic trend:** adding a $t^2$ term lifts the adjusted $R^2$ above 0.99 for every window.
- **Non-linear model:** a logistic trend with three sinusoids of free frequency and phase (12 parameters, shown below), fitted by Levenberg–Marquardt (`minpack.lm::nlsLM`) from interpretation-based starting values, together with a proposed self-starter.
- **Upstream stations:** lagged regressions of RPBU on the offshore stations OS11 and VR, choosing the lag that maximises $R^2$, to test the travel times quoted by Rijkswaterstaat.

$$Y_t = \frac{\beta_0}{1+\beta_2 e^{-\beta_1 t}} + \sum_{j=1}^{3}\alpha_j\sin(\omega_j t + \varphi_j) + \varepsilon_t$$

## Results

| Model (fit window) | 5-minute forecast | Observed | Error |
|---|---|---|---|
| Linear with quadratic trend (30 min) | 2.405 m | 2.370 m | 1.49% |
| Non-linear (30 min) | 2.389 m | 2.370 m | 0.80% |
| Non-linear (90 min) | 2.469 m | 2.370 m | 4.20% |

- The linear model fits closely, but the observed level falls outside its 95% intervals, which motivated the non-linear model.
- The non-linear model on a 30-minute window gives the best forecast. Its accuracy drops on longer windows, where Durbin–Watson tests also show significant residual autocorrelation.
- At 10-minute resolution, the estimated lags are 10 min for OS11 → RPBU (Rijkswaterstaat quotes 6), 30 min for VR → RPBU (23) and 10 min for VR → OS11 (17).

## Repository contents

- `src/part_3_1_high_tide_exploration.R`: finds the high-tide days and plots the water level around high tide
- `src/part_3_2_regression_models.R`: the linear, quadratic and non-linear models, their diagnostics and the 5-minute forecasts
- `src/part_3_3_upstream_stations.R`: cross-correlations and lagged regressions against the OS11 and VR stations

The scripts expect `data/RPBU20102011.txt` and `data/data20102011-10min.txt`, which are not included (`data/` is git-ignored). Packages: `zoo`, `minpack.lm`, and `car` for the Durbin–Watson test.

**Tools:** R (`minpack.lm`, `zoo`)

SISIR
================

<!-- badges: start -->
[![CRAN version](https://www.r-pkg.org/badges/version/SISIR)](https://CRAN.R-project.org/package=SISIR)
[![CRAN checks](https://badges.cranchecks.info/summary/SISIR.svg)](https://cran.r-project.org/web/checks/check_results_SISIR.html)
[![CRANLOGS](https://cranlogs.r-pkg.org/badges/SISIR)](https://CRAN.R-project.org/package=SISIR)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:98f3e762fa7fbf55fb72e3f2a6bcc7c5a115e505/)](https://archive.softwareheritage.org/swh:1:dir:98f3e762fa7fbf55fb72e3f2a6bcc7c5a115e505;origin=https://forge.inrae.fr/sfcb/sisir.git;visit=swh:1:snp:c69201df54de7da3973c3104a31386c0cd856e08;anchor=swh:1:rev:8568af1877fec2324de0529df8dea7675961912a)
<!-- badges: end -->

`SISIR` is an **R** package designed to handle function data (e.g., curved
sampled on a discrete time grid). Implemented methods can perform interval
fusion and selection procedures in regression models with functional inputs.
They include a semiparametric approach based on Sliced Inverse Regression
(SIR), as described in
<a href="https://dx.doi.org/10.1007/s11222-018-9806-6">doi:10.1007/s11222-018-9806-6</a>
(standard ridge and sparse SIR are also included in the package) and a random
forest based approach, as described in <a href="https://dx.doi.org/10.1002/sam.11705">doi:10.1002/sam.11705</a>.

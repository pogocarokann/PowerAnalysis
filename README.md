# PowerAnalysis

[leidwesen.github.io/PowerAnalysis](https://leidwesen.github.io/PowerAnalysis/)

JS offers a UI, and uses the [libRMath.js](https://github.com/R-js/libRmath.js) port of R. The library is ~0.5MB.

Python is faster, but must be run locally without UI (edit variables then run through CLI), requires numpy and scipy.

Rebuilt from the skeleton of PiFlavour's [PowerAnalysis](https://github.com/PiFlavour/PowerAnalysis).

## What's a Power Analysis ?

Supposing a discrete probability distribution `p`, we guess it may have changed to a distribution with different weights `ptest`.
The analysis's goal is to find a required sample size `N` that will be able to detect such a deviation with a certain amount of confidence.

A theoretical approach to the analysis would go as such:

0. Declare the confidence we wish to have.

1. Sample `N` times from a `ptest` distribution

2. Determine if the ptest sample is significantly differently from `p`. This is done by summing the chi-squared statistics for each term, then comparing to the threshold our confidence gives.

3. Repeat the test of steps 1&2 a predefined number `Ntests` of times (the more the longer it takes, but the less random the results are)

4. The accuracy of this sample size `N` is thus the # of significant tests / # of total tests.

This can be optimized in practice.

Instead of sampling `N` times then counting the results up, we can instead sample directly from a multinomial distribution with probabilities `ptest` and `x=N`.

The python code leverages numpy to vectorize all the chi-squared computations instead of looping.

For the chi-squared statistic, the traditional formula is `(O-E)^2/E`. `O` is the observed counts (what was sampled from `ptest`), `E` are expected counts under `p` , so `E = p*N`. This can be factored to obtain `N * (O/N - p)^2/p`. I'm not sure if this is faster, but it probably has some numerical stability by ensuring that large numbers aren't being squared/divided, and the `*N` can be factored outside the summation of statistics. 

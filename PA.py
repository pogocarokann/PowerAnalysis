# Power Analysis

# Given the (theoretical) discrete probability distribution p, we postulate a potential deviation ptest from this.
# Our goal is to find the required sample size so that we will "likely" (see below) be able to detect such a deviation

# This is what happens:
# 1] Draw Ntest random samples of size Z from ptest, precounted
# 2] Compute chi^2 of all tests w.r.t p
# 3] Determinine accuracy w.r.t. confidence level
# this can either output the accuracy for a given samplesize, or iteratively figure out the correct size to obtain at least the desired confidence within a given limit.

import numpy as np
from scipy.stats import chi2

################################ Constants (change as you see fit) #######################################################

# Verbose. whether to print alot of intermediate steps to the console.
verbose = True

# Confidence. This is how desired accuracy (ratio of signifcant tests for given p and ptest)
confidence = 0.95

# Numberof tests. The more tests, the more accurate the result
Ntests = 10000

# Precision. How close we want to get to confidence
precision = 0.0001

# Sample Size. Set this and turn singleOutput to true to only print the accuracy for the single N written here. Otherwise, the program will ignore this N and calculate the optimal N for the given confidence.
singleOutput = False
N = 25000

# The theoretical probability distribution which we assume to be the case (null hypothesis)
# For example, use [0.5,0.5] for a fair coin
null_p = 2/5000
p = np.array([null_p,1-null_p])

# The probility distribution that we might potentially expect. For example when some reports of changed rates come up.
# For example, use [0.45,0.55] if you suspect a coin might be slightly biased
test_p = 5/5000
ptest = np.array([test_p,1-test_p])

################################# Input Check (do not modify) ############################################################

def inputCheck(p, ptest):
    dim = len(ptest)
    if (len(p) != dim):
        raise ValueError("p and ptest do not have the same dimension")
    if( np.sum(p) != 1):
        raise ValueError(f"The probabilities in p do not sum to 1, but rather {np.sum(p)}")
    if( np.sum(ptest) != 1):
        raise ValueError(f"The probabilities in ptest do not sum to 1, but rather{np.sum(ptest)}")
    return dim
dim = inputCheck(p, ptest)

################################# Start calculation (do not modify) ######################################################

# Given samplesize N, and global constants ptest, p, Ntests
# Samples Ntests times from a size N ptest multinomial distribution
# Compares results to p multinomial distribution using chi-squared stat.
def accuracyFromSample(N, Ntests=Ntests, ptest=ptest, p=p, confidence=confidence, dim=dim, verbose=verbose):
    if (not dim):
        dim = inputCheck(p, ptest)
    counts = np.random.multinomial(N, ptest, Ntests)
    chisq_arr = N*np.sum((counts/(1.0*N) - p)**2/p, axis=1)
    chisq_sig = chi2.ppf(confidence, dim-1)
    sigCount = np.sum(chisq_arr > chisq_sig)
    acc = sigCount*1.0/Ntests
    if verbose: print(f"significance of {N} = {acc} using {100*confidence}% confidence.")
    return acc

# A "binary" search to find samples needed for confidence
# double until above accuracy
# delta = N/2
# then halve  delta, decrease until below
# then halve delta, increase until above
def sampleForConfidence(confidence=confidence, precision=precision, Ntests=Ntests, ptest=ptest, p=p, dim=dim, verbose=verbose,):
    if (not dim):
        dim = inputCheck(p, ptest)
    N = 1
    accuracy = 0.0
    delta = np.inf
    initIncrease = True
    isIncrease = False
    swapsAtOne = 0
    while (abs(accuracy - confidence) > precision or accuracy < confidence):
        if (initIncrease):
            if (accuracy > confidence):
                delta = N/4.0
                N = N-delta
                initIncrease = False
                if verbose: print(f"inital increase finished, d={delta}")
            else:
                N = N*2
        else:
            if ((isIncrease and accuracy > confidence) or
                (not isIncrease and accuracy < confidence)):
                delta = max(1,delta/2.0)
                if (delta<=1):
                    if (swapsAtOne > 2):
                        print(f"search got stuck around N={N} and accuracy={accuracy} searching for confidence in [{confidence}, {confidence+precision}], try running again with widened precision");
                        break
                    swapsAtOne += 1;
                if verbose: print(f"passed confidence, d={delta}")
                isIncrease = not isIncrease
            N = N + (1.0 if isIncrease else  -1.0) * delta
        accuracy = accuracyFromSample(N, Ntests, ptest, p, confidence, dim, verbose)
    if (swapsAtOne<=2):
        print(f"Number of Tests: {Ntests}")
        print(f"To differentiate between null {p} and {ptest} with {100*confidence}% confidence, use at least")
        print(f"Sample Size: {N} (accuracy = {accuracy})")

if __name__=="__main__":
    if singleOutput:
        accuracyFromSample(N, verbose=True)
    else:
        sampleForConfidence()

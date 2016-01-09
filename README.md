# Comparing the accuracy and precision of measurement modalities without gold standard data

## Introduction

This repository contains R and JAGS files that implement a version of a method described by Hoppin et al. in IEEE TMI 21(5) 2002. The problem addressed is how to assess the accuracy and precision of multiple ways (modalities) of estimating some quantity given a sample of such estimates, but in the absence of gold standard data. This seems intractable, but Hoppin et al. assume linear relationships with iid residuals between the true and estimated values for each modality and write those relationships in terms of the *distribution* of the true values. (Such a distribution could be obtained from the literature or derived using knowledge of the specific thing being estimated.) They then write a maximum likelihood expression for the parameters of the linear models in which they integrate out the unknown true values. The result is estimates of the slopes and intercepts of the line of best fit between the true and estimated values, and the standard deviation of the distributions of the residuals, for each modality. The paper also presents some simulations that demonstrate the method.

The R and JAGS code here implement a very similar method and replicates the first experiment in the paper. The main difference is that a Bayesian approach is taken. Rather than derive a closed form for the maximum likelihood described in the paper, or use numerical integration when that is not possible, here we use use JAGS to write the linear models directly and use MCMC sampling both to sample the posterior distribution of the slops, intercepts, and standard deviations, and to integrate out the unknown true values. The appeal of this implementation is that:

1. The JAGS language allows us to describe the problem in a simple declarative way. This makes it very simple to change the assumed distribution of the gold standard values being estimated (a one-line change in most cases, although it is also sensible to think about choosing sensible priors for the parameters of the linear models given the expected values, accuracies, and precisions of whatever quantity your specific modalities estimate).
2. The JAGS implementation of MCMC makes performing the numerical computations much simpler than if we had to hand-code it. Rather than having to carefully code the integration over the unknown true values, we can place reasonable trust in JAGS's Gibbs sampler to do that integration for us.

Note that while I think the implementation of the method is sound, the replication of the Hoppin et al. experiment is a proof-of-principle implementation. In particular, I have not looked in any detail at any diagnostics. I did briefly experiment with using more aggressive thinning given the degree of autocorrelation in the sampling (see the plots, discussed below), but it made relatively little difference to the results. However, you are encouraged to be more thorough if you use this code for a real analysis!

## Requirements

I have not checked the following requirements extensively, so you may need to use your own common sense. To run the code, I assume you have the following:

1. [RStudio](http://www.rstudio.com) (You can run this code without RStudio, but I think this is currently the best graphical R environment).
2. [JAGS](http://mcmc-jags.sourceforge.net) (see below for details on installing this).
3. The following R packages (follow these [instructions to install R packages using RStudio](http://quietube6.com/v.php/http://www.youtube.com/watch?v=u1r5XTqrCTQ)): runjags and its dependencies (just rjags, I think).

At the time of writing, runjags requires a 4.x version of JAGS.

If you are using a UNIX-like system, you can likely install JAGS from the terminal using whatever packaging system your distribution uses.

On Linux (specifically, Ubuntu, Debian), you can likely run

```
apt-get install jags
```

I have previously installed JAGS on a Mac using [Homebrew](http://brew.sh) via

```
brew install jags
```

but recent attempts to install via Homebrew were unsuccessful. Instead I used the standard Mac OS X installer available from the [JAGS SourceForge page](http://sourceforge.net/projects/mcmc-jags/files/). Windows installers are also available there.


## Running the example code

After installed JAGS and the required R libraries:

1. Using RStudio, open the [comparing-measurement-modalities-without-gold-standard-data.R](comparing-measurement-modalities-without-gold-standard-data.R) file using the File menu.
2. Click the "Source" button in the top-right of the pane containing the above source code.

You should see progress being made in the analysis (it takes about 30 s on my 2014-vintage MacBook Pro). A table will appear in the R Console pane, and graphs appear in the Plots pane.

## Understanding what just happened

To really understand what just happened, you should read the code. However, briefly, we:

* Simulated gold standard data according to the Hoppin et al. paper. In a real analysis this would not be available. It is used in the Hoppin et al. paper, and here, as a way to check that the method is able to recover the values of the parameters of the linear models that relate the true values to the estimates.
* Simulated observed from the gold standard data using three modalities.
* Plotted the simulated gold standard and observed data.
* Ran the version of the Hoppin et al. method on the simulated observed data.
* Summarized the posterior distributions of the linear model parameters and plotted simple diagnostics and the posterior distributions for each of the linear model parameters.

On my computer, the tabular result was as follows (some columns removed):

     |   Lower95 |    Median |   Upper95 |      Mean |        SD
-----|-----------|-----------|-----------|-----------|----------
a[1] |   0.51061 | 0.57114   | 0.63707   | 0.5724    |  0.032384
a[2] |   0.58167 | 0.64167   | 0.70729   | 0.64294   |  0.032488
a[3] |   0.69824 | 0.78931   | 0.88578   | 0.79093   |  0.047924
b[1] |  -0.12524 | -0.093647 | -0.061995 | -0.094383 |  0.016187
b[2] | -0.019156 | 0.015343  | 0.047074  | 0.01444   |  0.017148
b[3] |  0.041432 | 0.091301  | 0.13402   | 0.09009   |  0.023964
s[1] |  0.037298 | 0.04699   | 0.056831  | 0.047102  | 0.0049714
s[2] |    0.0242 | 0.036176  | 0.049022  | 0.03621   | 0.0062333
s[3] |   0.06573 | 0.079162  | 0.093629  | 0.079414  | 0.0071156

The first three rows are for the slope terms in the linear models for the three modalities. These describe how the modalities' estimates deviate from the true values as a function of the value itself. An ideal modality would have a slope term equal to unity. The correct values (used to generate the data) are 0.6, 0.7, and 0.8, respectively. The posterior medians and means are all reasonably similar to those values, and those values lie in the lower- to upper-95% intervals. Inference is perhaps slightly poor for the second modality, which is underestimated.

The second three rows are for the intercept terms. These describe the modalities' accuracies (biases). An ideal modality would have an intercept of zero. The correct values are -0.1, 0.0, and 1.0, respectively. The posterior medians and means are very close.

The second three rows are for the residual standard deviations. These describe the modalities' precisions. An ideal modality would have a standard deviation of zero. The correct values are 0.05, 0.03, and 0.08, respectively. The posterior medians and means are very close.

Several plots are produced (use the left and right arrow buttons in the Plots pane to see them). The first shows three scatter plots, illustrating the relationships between the gold standard and estimated values for the three modalities. Subsequent plots each show: a trace of the MCMC samples drawn from the posterior distribution of a given parameter of the linear model (top left; two chains were used, hence two colors); empirical cumulative distribution functions for the parameter (top right; again, two chains); an autocorrelation plot (bottom right); and a histogram of the posterior distribution for the parameter (bottom left). I'm not entirely sure what the final plot shows; perhaps the covariance matrix for parameters?


## Outline of the JAGS file

The [comparing-measurement-modalities-without-gold-standard-data.jags](comparing-measurement-modalities-without-gold-standard-data.jags) file describes the model. It has two section: a data section, which simply determines the number of observations and modalities, and the model section, which describes how the data are modeled. See the comments for full details, but briefly each observation is modeled using the linear equation proposed by Hoppin et al., we specify a prior for each of the unknown parameters, and also the distribution of the value that the modalities estimate. The latter is the key characteristic that needs to be adapted to new analyses and the Hoppin et al. paper investigates what results can be expected if this distribution does not agree well with the actual distribution. It may also be necessary to adapt the priors on the linear model parameters to reasonably match the values you are trying to estimate.


## Parameterization of distributions used in Bayesian analyses

It is crucial to note that it is typical in Bayesian analyses to parameterize distributions differently to how they are commonly parameterized. For example, in Bayesian analyses the normal distribution is often parameterized by mean and precision (where precision is the reciprocal of variance), rather than mean and standard deviation. This is because when you do Bayesian analyses algebraically (where that is possible), you often end up working with quantities such as inverse variance. Therefore it makes sense to re-parameterize distributions to make the notation simpler. You should read the [JAGS manual](http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/) to check exactly how the distribution you want to use is parameterized, and probably plot it to check that it has the shape you think it should.


## A note on the JAGS language

JAGS derives from the [BUGS project](http://www.mrc-bsu.cam.ac.uk/software/bugs/) and uses essentially the same model specification language. There is a fundamental difference in the interpretation of the language that JAGS uses compared to most common programming languages. Ignoring multithreading, in a language like R (or C, C++, Java, etc.), each statement in a program is (conceptually at least) executed in turn, where the state of the program after a statement depends on the current state of the program. In JAGS, this is not true: instead, the entire model is read and converted to an internal format that is used to determine how sampling should be performed. Even though the JAGS language has a **for** loop, this should not be interpreted in the same way as a **for** loop in conventional languages: it is really just a convenience to avoid having to write similar parts of a model specification many times. So, do not think of a JAGS file as a "program", but rather as a description of a model.


If you want to learn more about Bayesian estimation, I suggest two books:

1. [Bayesian Statistics, by the Open University Course Team](http://www.amazon.co.uk/dp/B00IJ0ORL6). This is a very solid introduction that focuses on the basic mathematics of Baysian statistics and later introduces the use of computation to perform estimation.
2. [Bayesian Data Analysis, by Andrew Gelman, John Carlin, Hal Stern, David Dunson, Aki Vehtari, and Donald Rubin.](http://www.stat.columbia.edu/~gelman/book/). This is arguably *the* book on Bayesian data analysis, but I would suggest starting with the Open University book first.



#!/usr/bin/env python
# coding: utf-8

# # Lab 4:  Fitting a Model to Data
# 
# In this lab, we will begin to build various tools to fit a model to data.  The goal is to understand, implement and compare chi2 and MCMC fitting routines to fake data.    In Lab 5, we will apply these tools to real spectroscopic data from Keck / DEIMOS.
# 
# The goals of the lab are:
# 
# 1. Use python modules to perform $\chi^2$ fits on data
# 2. Write and compare to your own $\chi^2$ algorithm
# 3. Explore sampling algorithms
# 4. Write an MCMC algorithmm and compare to $\chi^2$ results
# 
# 

# ## Part 1:   $\chi^2$ Fitting
# 
# In this lab, we will largely follow two papers — one written by D.W. Hogg (NYU/CCA), Jo Bovy (U Toronto), and Dustin Lang (Perimeter), and the other by Hogg and Dan Foreman-Mackey (CCA).
# 
# The first of these is accessible here:
# https://ui.adsabs.harvard.edu/abs/2010arXiv1008.4686H/abstract

# ### Question 1
# 
# Do Exercise 1 from Hogg (2010) using a standard python routine of your choice.   The data (Table 1) are available in the A330 public page under Data Access as well as in your Git directories (as the file is small).   Please report your best fit values and errors on these values.  These should be very close to those in Figure 1 of the paper.
# 
# ```{tip}
# If you are using np.polyfit, set the `cov=True` keyword in polyfit to have it return not only the fit coefficients, but also the full covariance matrix. The parameter uncertainties, assuming no off-axis covariance, are the square roots of the diagonal terms (of which for a linear fit there will be 2. You can pull the diagonal terms of a square array using `np.diag()`.
# ```

# In[ ]:





# ### Question 2
# 
# Repeat the question above, however, this time write your own script to solve this problem by evaluating chi2 on a grid of m and b values. You should write a `chi2()` function that reads in 4 arguments `m,b,data_x,data_y,unc_y` (though you can call them what you want). You should then write a `fit_line()` or `minimize_chi2()` function that will, across your grid of $m$ and $b$ values, evaluate chi2 using your `chi2()` function. You may use the values above to guide your grid, making sure the grid spans at least 2-sigma in both directions.
# 
# Plot the chi2 values for all grid points. We suggest creating an array, `chi2_image`, which is a shape characterized by the lengthd of your `m_grid` and `b_grid`s. Then, as you double-for-loop over `m` and `b` values and calculate `chi2`, you can set `chi2_image[i,j]` to the output chi2 value. 
# 
# ```{tip}
# Remember, $m$ and $b$ values won't index your array directly. So you'll want to loop via something like `for i,m in enumerate(m_grid):` and `for j,b in enumerate(b_grid):` if you're going to do that.
# ```

# While chi2 fitting is reliable for determining the best-fit values, it is not always easy to estimate errors on these parameters. For example, in the above example, we had to explicitly initialize a grid of parameters to fit on, and as soon as this grid has to get finely spaced, or moves into any number of dimensions > 2, everything gets much more computationally expensive to calculate, and understanding the chi-squared "surface" in multi-D becomes difficult. Additionally, we had to narrow in our range of $m$ and $b$ values to get it to work, but there may actually be a better solution elsewhere in parameter space that we're not accessing. 

# ### Question 3
# 
# Determine the best fit parameters and one-sigma errors from Question 1.2.  The best-fit value can either be the minimum chi2 value or (bonus) by fitting a function to your chi2 values and interpolating the best fit.
# 
# Determine the 1-sigma errors on your best-fit parameters. by noting the surface where chi2 = chi2 +2.3

# ## Part 2:   MCMC Fitting
# 
# While chi2 is a good method for determining best-fitting values, it less reliable in determining errors on those parameters.   If your science question requires good error estimates and/or if your model contains more than a few parameters, Monte Carlo (MCMC) is a popular tool.
# 
# 
# https://ui.adsabs.harvard.edu/abs/2018ApJS..236...11H/abstract
# 

# You will need to install two packages for this work inside of your A330 environment:
# 
# ```
# conda install emcee
# conda install corner
# ```

# ### Question 4
# 
# Read Hogg et al. 2018 and do Problems 1-4.   
# 
# For Problem 1, you are welcome to explore just the mean and varience.
# 
# For Problem 2, you have no choice.  Use python :)
# 
# For Problem 4, I found it easier to do 4b first, then 4a.

# ### Question 5
# 
# While the above problems should give you a sense for how MCMC works, most reseach problems use standard packages to run MCMC.  MCMC packages in astronomy include emcee, MultiNest,  Dynasty.   
# 
# 
# Write an MCMC to evaluate the data in Question 1+2 above using emcee.   We suggest checking out the guide to MCMC here:
# 
# https://prappleizer.github.io/Tutorials/MCMC/MCMC_Tutorial.html
# 
# 
# We suggest 20 walkers and 2000 steps for your sampler.    Plot both the sampler chains and a corner plot of the results.   
# 
# Compare the best fit values for m and b from the chi2 and MCMC, as well as their errors.   

# In[ ]:





# In[ ]:





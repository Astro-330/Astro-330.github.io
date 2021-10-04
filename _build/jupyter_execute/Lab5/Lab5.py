#!/usr/bin/env python
# coding: utf-8

# # Lab 5: Stellar Spectroscopy + Fitting Models to Data
# 
# If imaging data is the 'bread and butter' of astronomy (see Lab 2), then spectrosopy is meat and potatoes.    
# 
# In this lab, we will guide you through reading, plotting and fitting spectra of stars in a Milky Way globular cluster.  The science goal is to determine the velocity and velocity errors for a handful of stars in order to determine if they are members of the globular cluster, or foreground stars in the Milky Way.    The coding goal is to apply both $\chi^2$ fitting and MCMC fitting techniques when the model is more complicated.
# 

# ### Goals of this lab:
# 
# 1. Explore a maintained software package (pypeit).
# 2. Read a complicated fits file and plot a spectrum. 
# 3. Find parameters and errors via chi2 fitting when the model is not an analytic function
# 4. Find parameters and errors via MCMC.
# 5. Fitting polynomials to 2D surfaces, corner plots
# 

# ### Question 1:   Keck DEIMOS data

# We will be working with data from the Keck Telescope's DEIMOS instrument.   All Keck data is publically available on the Keck Data Archive (KOA) website.   While we will not be directly reducing the raw data, let's take a look at these files to get a sense for what the data look like.   We've selected data from the Milky Way globular cluster NGC 7006. 
# 
# Head to the KOA website (https://koa.ipac.caltech.edu/cgi-bin/KOA/nph-KOAlogin) and search for all files take with DEIMOS on the night of June 3, 2011 (20110603).   Search the list for files with `Target Name == n7006` and `Image or Dispersion == mos`.    Find the column named `Quicklook Previews` and click on `[Raw]`.   This is a single exposure of a spectroscopic mask centered on NGC 7006.   You should see a hundred or so spectra in this image.   

# In[13]:


#Include a screen grab here and/or describe in words the image that you  see.


# ### Question 2:  Spectral Reductions with PypeIt

# Using the raw files downloaded from the KOA above , we [the A330 instructors] have run the science and calibration frames through a spectral reduction softare package called `PypeIt`:  https://pypeit.readthedocs.io/en/release/.  The `PypeIt` github repository can be found here:  https://github.com/pypeit/PypeIt

# While we won't actually run `PypeIt` in this lab, we will be using its output files.   This is a software project that is actively being developed, so let's look around at the code and identify some familar pieces:
# 
# On github, take a look in the /pypeit directory and click on a few of the *.py files.  
#    1.  Find one instance of PypeIt using a Class structure  
#    2.  Find one instance of PypeIt not fully/properly populating a doc string :)
#    3.  Find a line of code that you understand and explain what its doing
#    4.  Fine a line of code that you don't understand.
#    5.  How many branches current exist from the main `release` branch?
# 

# In[14]:


# Answers to 5 items above.


# In the data access directory, we have provide a PypeIt output file which contains one-dimensional spectra for all the stars observed in the DEIMOS mask `n7006a` that you viewed above.   Read in the file using the astropy.io.fits commands and view the contents using `hdu.info()`.   State how many spectra are contained in this file.

# In[22]:


from astropy.io import fits
file = 'spec1d_DE.20110603.45055-n7006a_DEIMOS_2011Jun03T123053.021.fits'


# In[24]:


# Code to view file contents
# How many spectra are contained in this file?


# ### Question 3: Plotting 1D PypeIt output spectra and fitting by eye

# We have selected 3 spectra from this file which are high signal-to-noise stars.   From your fits table that you have read in, select extension 121, 135 and 157.   These can also be selected using the names 'SPAT0564-SLIT0560-DET06', 'SPAT1163-SLIT1162-DET06' and 'SPAT0288-SLIT0302-DET07'.   Save the data for each spectrum separately.
# 
# Plot wavelength versus counts/flux for each star.   Please use the optimal extraction results ('OPT_*').  If you need additional guidence for what is in this file, see the PypeIt documentation on spec1d_* files:  https://pypeit.readthedocs.io/en/release/out_spec1D.html

# In[25]:


# For each of these three stars, plot the wavelength versus counts.  Use ylim =  8300-8800 Angstrum


# ### Extra (+0.5)   
# To get a sense for the velocity of each star, you might try measuring a rough velocity 'by eye'.   The three strongest lines in the spectra above are from Calcium II:  8500.36, 8544.44, 8664.52 Angstrum. What velocity do you estimate?   

# ### Question 4: Synthetic model spectra 

# In ASTR 255 and the extra question above, you have measured the velocity of a star by measuring the center of a known absorption line (either by eye or fitting a Gaussian) and comparing to its rest wavelength.   While this process does estimate the star's velocity, it wastes much of the information present in the full spectrum.  To determine more accurate velocities, we turn to "template fitting" where a spectrum with a known velocity is compared to our unknown science spectrum.    A template spectrum can either be empiricial (an observed spectrum of a standard star where the velocity is already known) or synthetic (numerically computed from stellar models).   Here we will use synthetic templates from the PHEONIX library:  https://phoenix.astro.physik.uni-goettingen.de/

# In[ ]:


template_file = 'dmost_lte_5000_3.0_-2.0_.fits'


# In[ ]:


def read_synthetic_spectrum(pfile):
    '''
    Function to load synthetic template file into python using vacuum wavelengths
    
    Parameters
    ----------
    pfile: str
        path to the synthitic fits file to load. 
        
    Returns
    -------
    pwave: float array
        Wavelengths of synthetic spectrum
    pflux: float array
        Flux of sythetic spectrum
    '''

    with fits.open(pfile) as hdu:
        data     = hdu[1].data
        
    pflux = np.array(data['flux']).flatten()
    awave = np.exp((data['wave']).flatten())
    

    # CONVERTING AIR WAVELENGTHS TO VACUUM
    s = 10**4 / awave
    n = 1. + 0.00008336624212083 +             (0.02408926869968 / (130.1065924522 - s**2)) +            (0.0001599740894897 / (38.92568793293 - s**2))

    pwave  = awave*n
    
    return pwave, pflux


# In[ ]:


# Read in synthetic spectra and plot wavelegth versus flux


# ### Question 5: Synthetic model spectra -- Smoothing and Continuum fitting

# We will fit the sythetic spectrum to our science data with the goal of determining the velocity of our science spectrum.   The synthetic spectrum is at zero velocity.   To match the science data, we will need to (1) smooth the synthetic spectrum to the wavelength resolution of the science, (2) shift the sythetic spectrum to the velocity of the science data, and (3) rebin the synthetic spectrum and match continuum levels.

# #### Smoothing the templates 
# 
# We will first address how to smooth the synthetic spectrum to match the data.   We will fit for this value below, but for the moment, let's just choose a number based on a single point estimate. The DEIMOS spectral lines are well fit by a Gaussian with a 1-$\sigma$ line width that is roughly 0.5 Angstrum.   The synthetic spectra have resolution of 0.02 Angstrum.   Thus, we need to smooth the sythetic spectra with a Gaussian kernal that is 0.5/0.02 = 25 pixels.   
# 
# Hint: scipy has functions which do Gaussian filtering in 1D.

# In[26]:


# Write a function to Gaussin smooth the synthtic spectrum, using a smoothing kernal of 25 pixels.


# #### Fitting the Continuum
# 
# We will next address the above step (3), the overall shape and value of the spectrum which we will call the 'continuum'.   Let's fit a function to the synthetic spectrum so that it is approximately the same as the science spectrum. For the section of a spectrum we are working with a **linear function** (i.e., like we fit in lab 4) is sufficient. To do this, we will first rebin the synthetic spectrum in wavelength to the same array as the data. 
# 
# Choose a science spectrum from above and rebin the sythentic template so that it uses the same wavelength array (consider using `np.interp()`). We need this to carry out point by point fits and comparisons between arrays. 
# 
# Next, determine the **linear function** (mx+b) needed to match the continuum of the synthetic spectrum to that of the science. If you wish, you may also try a second order polynomial, if you think it is a better fit.
# 
# ```{tip}
# If you just try to fit the spectrum, the absorption lines and other features will "drag" your fit away from the level of the continuum. This is easy to see by eye. There's a few ways around this. First, we could mask regions of deep absorption. Second, we could run something like `np.percentile()` on the spectrum, and remove all points farther than, say, 1 or 2-sigma from the median when doing the fit. Third, we could do an iterative fit (see the extra below). 
# 
# For this problem, we'll allow you to estimate the linear fit to the continuum by eye, or by any of the methods above. 
# ```

# In[29]:


# Write a function to rebin the synthetic template to the data wavelength array and fit the continuuum.


# OK, now run both functions (your smoothing function and your rebin/continuum function) on the sythetic spectrum and plot the results.  

# In[1]:


# Run both functions (smooth + rebin/continumm) and plot your smoothed, continuum normalized synthetic spectrum
# Compare this to one of your science spectra.


# ### Extra (1.0) 
# When fitting continua, we usually want to avoid the "features" in the spectrum. We could mask them out, or drop percentiles of the data far from the median... but we could also iteratively remove them. To do this, you would fit your chosen function to the data as a start, then iterate, throwing out 3 (or 5 or whatever works) sigma distant points and re-fitting. This works because emission and absorption lines have data points far from the continuum value. Try fitting your continuum this way to get a better estimate. 

# In[ ]:


# iterative fit method


# ### Question 6: $\chi^2$ fitting to find velocity

# The science and synthetic spectra above should roughly match-- except for an unknown velocity shift.   You can shift the synthetic template in velocity by changing its wavelength array *before smoothing*.  Recall that $\delta \lambda = \lambda * v/c$.   

# Write a $\chi^2$ code to find the best-fit velocity for each of the three stars above.   Look up the velocity of the globular cluster NGC 7006 to justify the range of velocities to search over.   Consider the wavelength resolution of your science data to determine the spacing of your grid.

# In[31]:


# Write a chi2 algorithm to determine the best fitting velocity and error.


# ### Question 7:  $\chi^2$ fitting with more parameters

# In Question 6, we fixed the smoothing value to 25 pixels and used a single function to match the sythentic to science continuum.   Next, let's redo $\chi^2$, but now including these values in the fit.   This will be a 2 parameter chi2 fit. 

# In[ ]:


# Repeat $chi^2$ fitting searching over 2 (and bonus 4) parameters: 
#             velocity, smoothing, and (bonus) continuum value (m,b)
# If you use 4 parameters, this will get ugly.  
#
# Calculate errors from your chi2 contours on the velocity only.


# ### Question 8:  MCMC with to find velocity

# Repeat Question 7 but this time fitting with MCMC.  We suggest writing a single function `make_model`  which creates a single synthetic model spectrum given an input velocity and smoothing.
# Report your best fit velocity and errors.
# 
# You can chose to fit 2 parameters (velocity and smoothing), or as a bonus all 4 parameters (velocity, smoothing and continuum fit values).

# In[33]:


# MCMC to find velocity only.  Report your best fit velocity and errors.
# Plot full corner plots for all fitted parameters.


# ```{note}
# In the context of MCMC, you'll often hear people talk about "marginalization". This is a classic example. Marginalization is the process of fitting for parameters we care about, plus "nuisance parameters" that we don't (like the smoothing and continuum values), and then "marginalizing out" the nuisance parameters by taking the 1D posterior spread only of the parameter of interest.
# ```

# ### Question 9:  MCMC convergence

# Confirm that your MCMC above converged and that you are discarding the appropriate number of samples when determining your parameters (that is the burnin number).

# In[34]:


# Confirm convergence 


# ### Question 10:  Science
# 
# And finally, some science questions:
# 1.  Do velocities agree between chi2 and mcmc within error?
# 2.  Are the velocity errors the same?
# 3.  Are these three stars part of NGC 7006?

# In[ ]:


# Answers to the 3 questions above


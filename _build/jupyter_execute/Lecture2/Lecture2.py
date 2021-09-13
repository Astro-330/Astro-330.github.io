#!/usr/bin/env python
# coding: utf-8

# # 2. Object Oriented Programming 

# In this lecture, we'll go over the basics of writing object oriented code. Because this is a large shift in programming syntax (and thinking), we'll be spending several lectures building up intuition for it. 

# ## Part I: Methods and Attributes

# *Writing* classes can be intimidating at first, and doing so involves some general boilerplate that can be confusing at first. 

# *Using* classes, on the other hand, is generally fairly intuitive. One reason for this: 
# 
# ***Everything in Python is an Object***. This means we are fairly used to using objects, even if we don't know it.

# Let's begin with a simple example. I'm going to initialize a basic `numpy` array. 

# In[2]:


import numpy as np

arr = np.ones((5,6))*15.5 
arr


# This simple array, `arr`, is an *object*. Objects in python have both *attributes* and *methods*. Both are accessed via what we call *dot notation*. Using dot notation, we can access an attribute of an object as follows:

# In[3]:


arr.shape


# `shape` is an attribute of all `numpy` arrays. By executing `arr.shape`, we retrieved this attribute. For documented classes, available attributes are generally listed in the documentation. The are also generally available for tab completion. 
# 
# Attributes are called as such because they are not *callable*. We can query them, but they are more like "properties" of the object, which are often used to make our code more robust (e.g., instead of hard-coding numbers into a code, we can flexibly adjust based on the sizes and shapes of the arrays within).

# Alternatively, `methods` of an object are callable functions that the object knows how to run, generally upon itself (with or without additional information). Returning to our `arr` array, we can call a method that reshapes the array:

# In[9]:


arr.reshape(6,5)


# We use the same dot notation as for attribute access, but we then give parenthesis and arguments as necessary, as we would when calling a function. We can see here that the output of this method is a reshaping of the array. Most Python objects (especially the basic datatypes, such as `lists`, `arrays`, `dictionaries`, and even `floats`, `ints`, and `strings`) have many useful attributes and methods attached to them, and knowing them can save you time in your coding. 

# In examples previous, we've seen plenty of quick method and attribute use --- for example, taking a string and running `.startswith()` on it, or `.strip()`. 
# 
# But how do we give ourselves the ability to write `classes` which give us this sort of attribute/method access flexibility?

# ## Part II: Defining Classes

# Classes are defined via the `class` call, and all contain, at minimum, a method known as `__init__`. 

# In[10]:


class MyObject():
    def __init__(self):
        pass


# Methods of a class are defined just like normal functions, but are all indented within the block of the class. They also all contain `self` as the first argument of the method. Why is this?

# When we created an array above, we did something known as *instantiating an object*. Somewhere in the `numpy` library is a master template for what all `numpy arrays` look like, and we put in some arguments to determine what our specific instance of the class would look like. (In our context, this is the value of the array, 15.5, and its shape). 

# After this, we were able to to call methods of the class via its name (that we set it to in our code) plus the dot notation. The mapping of the *internal* methods of the class to the *external* method call is handled by the use of the `self` argument. 

# As long as we remember to set `self` as the first argument of every method (followed by any other arguments you would typically add), our classes should work properly. There's one other critical element of using self:

# In[11]:


class MyObject():
    def __init__(self,parameter):
        self.parameter = parameter       


# We've now given our class it's first initialized attribute, `parameter`. This parameter is fed in when we *instantiate* the object, not as any additional call: 

# In[12]:


obj1 = MyObject(parameter=5)


# In[13]:


obj1.parameter


# By setting the *class attribute* `self.parameter = parameter`, we make it accessible both anywhere within the scope of our class (even in other methods) as well as outside, as a dot-accessible attribute. Any variable otherwise created in any class method that doesn't have `self.` in front of it won't be accessible outside the class as an attribute. Notice we never make reference to `__init__()` when interacting with the object from the outside.

# Let's move to a slightly more physical example:

# In[16]:


class Body():
    def __init__(self,name,mass,radius):
        self.name = name
        self.mass = mass
        self.radius = radius


# In[17]:


import astropy.units as u
sun = Body('Sun',1*u.Msun,1*u.Rsun)


# In[18]:


sun.mass


# In[19]:


sun.radius


# Our new `Body` class takes the inputs to the initialization function and decides that they should immediately be class attributes (accessible anywhere in the class or outside, as demonstrated above.) 
# 
# We need not set every input to `__init__` to *be* a class attribute. For example, if a class reads in a string pointer to some `fits` file, we may load the file and then not care about having access to that input string anymore. 

# ## Creating Methods

# The init method just handles stuff we need to put into to the class to get it up and running. Let's define our first "real" method. 

# In[26]:


class Body():
    def __init__(self,name,mass,radius):
        self.name = name
        self.mass = mass
        self.radius = radius
    def calc_density(self):
        vol = (4/3.) * np.pi * self.radius**3
        return (self.mass / vol).to((u.g/u.cm**3))


# In[27]:


sun = Body('Sun',1*u.Msun,1*u.Rsun)
sun.calc_density()


# You may be thinking that we could, technically, define a class *attribute* called density since it depends only on the previously defined parameters. It's true that we could've done this, but, in this case, changing the mass attribute wouldn't propogate to the density. By forcing the calculation each time, we ensure that changes to the mass or radius then are reflected in a *query* of the density.

# ### Exercise
# 
# Copy the class below, add an input attribute $x$, representing position. Assume it is a vector containing the x, y, z position of the object. 
# Then, add a method to the class `calc_grav`, which takes as input another instance of a `Body` class, and uses its own position and mass, and the other body's position and mass, to calculate the force $F_{\rm grav}$ between them via Newton's law. Create two bodies and test your methods!

# In[ ]:


class Body():
    def __init__(self,mass,radius):
        self.mass = mass
        self.radius = radius
    def calc_density(self):
        vol = (4/3.) * np.pi * self.radius**3
        return self.mass / vol


# # PSF Photometry

# We're all familiar (from 255) with aperture photometry. Aperture photometry sets a radius within which to sum the flux from some astronomical source. 
# 
# One question to ask when performing aperture photometry is what radius to choose. Too large, and you introduce noise from pixels that are sky. Too small, and you miss some flux from the actual star. 
# 
# We also know that the true "point source" of a star is spread out by the point spread function of the telescope/optics/detector. If we can *estimate* the shape (width) of the PSF in an image of sources, then we can construct a model which weights each pixel based on it's estimated fractional contribution to the stellar flux. Thus, pixels near the center of the stellar distribution are summed in full, while those near the edge are not included (or only just barely).

# To construct a pipeline that uses this technique, we'll need to 
# 
# - Read in an astronomical image with stars (and other point sources)
# - Perform the releant cleanings and calibrations to the image
# - Algorithmically locate the peaks (point sources) in the image, and carefully calculate their centers 
# - Estimate the PSF and construct the normalized model
# - Multiply the weight model by each source and calculate the sum using a weighting formula (using cutouts and our known positions) 
# 

# *Jamboard Discussion*: Based on the list above, discuss some reasons why moving to `classes` might be beneficial in trying to create such a pipeline. 
# 
# bit.ly/330jamboard

# ### Formulas Needed 

# **Centroid** 
# 
# $$
# x_{\rm com} = \frac{\sum x_i \hat{f}_i}{\sum \hat{f}_i}
# $$
# 
# $$
# y_{\rm com} = \frac{\sum y_i \hat{f}_i}{\sum \hat{f}_i}
# $$
# 
# where the $x_i$ are the pixel numbers and $\hat{f}_i$ are the fluxes in those pixels. This formula amounts to a light weighted mean --- pixels with larger fluxes will push the mean toward their pixel value. (This is also just a center of mass formula). 

# **2nd Moment** 
# 
# The centroid is the first moment of the light distribution. The second moment gives us an estimate of the width of a distribution. (To use stats terms, a mean is a first moment measurement, the standard deviation is a second moment measurement). 
# 
# We will use the following formulas in the lab to determine the x and y widths of distributions: 
# 
# $$
# \sigma_x = \left[\frac{\sum x_i^2 \hat{f}_i}{\sum \hat{f}_i} - x_{\rm com}^2\right]^{1/2}
# $$
# 
# $$
# \sigma_y = \left[\frac{\sum y_i^2 \hat{f}_i}{\sum \hat{f}_i} - y_{\rm com}^2\right]^{1/2}
# $$

# **PSF Photometry** 
# 
# Finally, the actual step of PSF photometry is performed via 
# 
# $$
# f_{\rm PSF} = \frac{\sum \hat{f_i} p_i}{\sum p_i^2},
# $$
# 
# where $\hat{f}_i$ are our pixel fluxes from before, and $p_i$ is our model for the PSF. One can model the 2D PSF with several functional forms. The Moffat profile is one of the closest matches to real stars, but it requires more than one shape parameter. Since our second moment only gives us one, we'll use a 2D Gaussian in the lab.

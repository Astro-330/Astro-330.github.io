#!/usr/bin/env python
# coding: utf-8

# # Lab 7: Package Structure

# In this lab, we're going to be doing less coding, and more exploration of how projects/packages are organized. We'll look at several packages and evaluate how the developer designed the UI (user interface) for the code. 

# ## Problem 1: Artpop
# 
# Artpop (co-written by Johnny Greco and Shany Danieli) is a package that lets you simulate images of a galaxy or stellar population under realistic observing conditions. There are several pieces that go into this, from the creation of a stellar population from isochrones, spatially distributing those stars, and the "simulating" the obsevation under different conditions and bands. 
# 
# ### Problem 1.1 
# 
# Navigate to the [documentation](https://artpop.readthedocs.io/en/latest/index.html) for Artpop. Follow the instructions to install it in your a330 environment using the normal pip install. 
# 
# Next, go to the **Quick Start** Tutorial, and follow along with it in your own notebook. Write a function that reads in an artpop `src` object and does the imaging and observing steps shown (we'll leave those alone). Then, create 3-5 sources with different stellar population inputs than the example provided and plot them up. We'll give a bonus point for the most interesting looking galaxy created. 
# 
# Below, I've shown an output I made.
# 

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.visualization import make_lupton_rgb
get_ipython().run_line_magic('matplotlib', 'inline')

# Project import
import artpop

# Code here


# In[ ]:





# ### Problem 1.2 
# 
# Returning to the home page, look at the list of Modules at the bottom. These are the outermost modules of the package. This is easier visualized by looking at the [Github Repo](https://github.com/ArtificialStellarPopulations/ArtPop) directly. 
# 
# Go to the repo and click into the `setup.py` file. You'll see there's a bit more code there than in our in-class example, but the actual setup block looks similar, with some additional dependancies listed. 
# 
# In the cell below, write a setup block based on the one from class / the one here for artpop, but for YOUR final project idea. You won't know your dependancies and such yet, but fill in whatever you can. 
# 
# 

# In[ ]:





# ### Problem 1.3 
# 
# Now click into the folder labeled `src/artpop`. This is where the actual source code for artpop is stored. You'll see several folders and several files, along with the all important `__init__.py`. Johnny and Shany have elected to have some of the classes/functions directly as files here, while some others are inside submodules with a few files and their own init files. 
# 
# Click into the `image` folder and copy/paste below the import statements listed in that submodule's init py. Based on the form of imports there, how would the user import the `ArtImager` class from `imager.py`? Write out a valid import that a user could use. 
# 
# Next, go back to the outer directory (`src/artpop`) and look at the init file there. This one is a bit more complicated, but only the imports at the bottom matter to us. Based on THOSE imports, is `ArtImager` (stored at `artpop.image.imager.ArtImager`, formally) directly importable as `from artpop import ArtImager`? Why or why not?
# 

# In[1]:





# In[2]:





# ### Problem 1.4
# 
# Based on your answer above and looking back at the lecture slides, list 2 pros and 2 cons of using the particular import schema Artpop uses in its `__init__.py` files. 

# In[ ]:





# ### Problem 1.5 
# 
# ArtPop has some tests set up in the `/tests` folder. Click into it, and click on `test_space.py`. This suite of tests is set up to ensure the spatial star distributions are being created correctly by the code. 
# 
# The test consists of one class with three methods. Explain what each is doing in a sentence. Then pick one of the latter two to copy below and explain in more detail what it is checking and how. You may find the documentation for [unittest](https://docs.python.org/3/library/unittest.html) helpful. 
# 
# Unit testing is not required for your final project, but having a few tests is highly recommended (and will be worth some bonus credit). 

# In[ ]:





# ### Problem 1.6 
# 
# In class, we learned that in an `__init__.py`, you can say `from .foo import *` to import all classes/functions/variables from that script. It turns out, we can actually get a bit more control over what gets imported. In the Artpop Github, navigate into `/image/imager.py`. Near the top, you'll see the following:

# In[ ]:


__all__ = ['IdealObservation', 'ArtObservation', 'IdealImager', 'ArtImager']


# The `__all__` command is special, and actually defines "all" the things that are going to be seen by the `*` import statement above when pulling from, e.g., `imager.py`. If an `__all__` is defined, only the things listed will by imported by default.
# 
# To practice this, go ahead and find your `exampy` package from lecture (or whatever you called it), and in the `imaging.py` file where you pasted `implot()`, add an `__all__` as shown in Artpop that lists it. 
# 
# At the end of this lab, you'll copy your package (the outermost directory) into your github repo for this assignment so we can go look at it.

# In[ ]:





# ### Problem 1.7
# 
# Using the provided documentation and tutorials, create a galaxy with artpop that has stellar populations similar to question 1.1 above..... but in the shape of a ring! 
# 
# All the information and pieces needed to do this are there in the various docs on Artpop's page. As you search for what you need, think about the structure of the tutorials and guides in the context of your own upcoming projects. 
# 
# After creating your "ring galaxy" and displaying it, answer the following: What was one helpful and one confusing thing about the way artpop's docs are laid out? Is there anything you would change or add? 

# In[32]:


# Your ring should look something like this 


# ### Bonus Problem (+1) 
# 
# Following more of the tutorials present, inject your ring galaxy into a real astronomical image background.

# In[39]:





# In[40]:


# Here's the blank sky we'll be injecting into 


# In[64]:


# OOOoooooOOOhhhhh, Ring Galaxyyyyyy!


# Tada!

# 
# 
# 

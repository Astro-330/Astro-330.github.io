#!/usr/bin/env python
# coding: utf-8

# # 7. Building Installable Packages

# In this lecture, we'll be talking about how to take code and turn it into something that can be installed on our system via `pip`, and `imported` anywhere on our system within Python. 
# 
# Note: Much of what I know about packaging up code was learned by reading Jo Bovy's helpful [resources](https://pythonpackaging.info/index.html) on this topic, and the example presented in this lecture draws heavily from the linked guide.

# ## Why Write a Package 

# A `package` is, in the abstract, just code. But most of the code we write (in particular, research code), doesn't constitute a package. *Generally* this is because our research code is written to handle a hyper-specific set of commands for a unique project. 
# 
# This doesn't mean we shouldn't make our research code publically accessible. It's scary, but as Hogg intimated during our conversation, it's good for everyone involved, and for the science. 
# 
# A package, meanwhile, is what happens when you discover some piece of code you wrote is more *generally* applicable, to either yourself later or to others.

# For your **final projects**, you will be making a package. By definition, you are thinking of topics that prompt the creation of a general set of tools such as a package. But I want to emphasize that even when you're just, say, writing research code, you may write a general function that you want to keep easily accessible for use later. Both are valid reasons to create a package.

# *Example* The `implot()` function we wrote a while back for quickly displaying astronomical images? That's a super handy function! I actually use it myself in my research code now. But rather than copying and pasting the function each time (which is both time consuming and prone to "version control" nightmares), I simply "installed" it on my computer.

# Finally: Well maintained, well documented, easy to use code that is hosted on GitHub is a *hallmark green flag* for anyone looking to hire a data scientist or computer scientist. Taking the bit of extra time to showcase your work in this manner is a great idea, and may even spawn new collaborations!

# *"In my opinion, the single-most important reason behind a code’s success in gaining wide-spread use is being **easy and intuitive to use**. This should be the main driver behind any decision regarding details of implementation, documentation, and installation of your code. Always put the user first" - Jo Bovy* 
# 
# *Hint: This is exactly why `emcee` has over 5000 citations (a pretty ridiculous number in astronomy*). 

# ## Package Design (abstract)

# Here are a few things to consider when making code a package for others to install and use:
# 
# - Limit dependencies as much as possible, ideally to just widespread ones like `numpy`. If you only use one function from a specialized library, can you just implement that one thing yourself?
# - Try to make any weird dependancies optional if at all possible, having `try: except` blocks to handle whether the user has them or not
# - Try to make installation as easy as `pip install packagname` or `git clone package_address; cd packagename; pip install .`. 
# - Have good documentation; have a readthedocs site if possible 
# - Keep a changelog somewhere public, especially for any major changes to the code's structure or use
# - In the docs, report which versions of python and numpy etc. your code has been most recently tested on (and try to keep it updated to the latest few versions) 
# - If possible (this is fancy), have a test suite that runs continuously when changes are pushed to your codebase, to check nothing breaks.

# ## Structure of a Package

# It's time for the nitty gritty. How do we actually go about making a package? Beyond basic python code, a package has three key components 
# - Structure (handled via folders and files within folders) 
# - `__init__.py` files which define the import behavior 
# - a `setup.py` file which tells standard install libraries how to install the package on a computer. 

# During this lecture, we're going to create a package together. This package is going to be a place for you all to keep handy functions/classes you write over time so you can easily get to them. 
# 
# I have one of these, I call mine `myutils`. For example, in any python script/terminal, I can say `from myutils.plotting import implot` to import that nice image plotter we wrote earilier this semester. 
# 
# Everyone decide now what you'd like to call your personal package. Make a folder somewhere on your system to store this package. I recommend somewhere easy to get to (I have a folder `/home/Software/` where I keep my own software as well as any I need to install from github rather than PyPI). There's going to be a top-level directory for your package that holds everything. The name of this one doesn't matter, it can also be your package name if you wish. 

# In[15]:


ls


# As we saw, I've made a folder `example/` that is going to serve as my example in this lecture. This `example/` is the top level directory; my package name for this is going to be `exampy`, a name I've shamelessly stolen from Jo Bovy. 
# 
# Inside your toplevel folder (for me, `example/`, we're going to have any README and documention files, a LICENSE (tells people how to cite/attribute our code and allowed uses), a CHANGELOG, and other ancilliary info. We will also, crucially, have our `setup.py` here. 
# 
# For now, let's focus on the code. Inside your toplevel folder, make a new directory that is the name of your desired package:

# In[2]:


import os 
os.mkdir('example/exampy/')


# In[4]:


ls example


# Great. This folder is going to store all the modules/submodules of our package. 
# 
# A simple package may have only one python file inside it, a more complicated one may have several files, and an even more complicated one may have sub *folders* with files. 
# 
# Each folder and file's name is going to contribute to how people use our code. For example, my `myutils` package *could* just be one script with a ton of functions in it. But I could also make a separate script called `plotting.py` and put `implot()` in there. Now instead of `from myutils import implot`, it's `from myutils.plotting import implot`. 
# 
# Is one better? Not really. For a personal utility library, it's probably fine to shove them all in one. But it can be nice to organize things, especially if you're trying to find where things are later to work on them. 

# Inside your package folder (for me, `example/exampy/`), we're going to need to make a file called `__init__.py`. It's allowed to be empty, so we could just do the following:

# In[5]:


get_ipython().system('touch example/exampy/__init__.py')


# In[6]:


ls example/exampy


# When someone types `import exampy` later, it's *that* init file that gets run. So we either need to put our functions/classes in there, or we can make separate files for them, and just import those files *inside the init*, leading to the same outcome. This is the standard practice, because it makes it easier for a user to navigate your github and find relevant functions.

# For the sake of this example, let's say we know we want to have a completely different submodule, `plotting`, where we'll put all the handy plotting functions we come up with over the years. 
# 
# Let's further say that we predict we might want multiple python files with different sub-ideas within `plotting`, for example, `imaging` and `histograms`. 
# 
# This tells us `plotting` should be a folder-level name, with `imaging` and `histograms` being `.py` files inside. (Again, this is just a choice). 
# 
# In your `exampy/` folder, make a new folder called `plotting/`, and within **it**, make two empty python files, one called `imaging` and one called `histograms`

# In[7]:


os.mkdir('example/exampy/plotting')


# In[8]:


get_ipython().system('touch example/exampy/plotting/imaging.py')


# In[9]:


get_ipython().system('touch example/exampy/plotting/histograms.py')


# For every folder-level thing we make in a package, we need an `__init__.py`. Add one of these to your `plotting/` folder. 

# In[10]:


get_ipython().system('touch example/exampy/plotting/__init__.py')


# In[11]:


ls example/exampy


# In[12]:


ls example/exampy/plotting


# Awesome! Believe it or not, this is actually pretty much it when it comes to the overall structure of things within a simple package.
# 
# The next question is, what do we put in each `__init__.py`? This turns out to depend on how we want the import statements to look when people use our code. 
# 
# Let's start with the inner-nested init. What if we want the organizational convenience of having histogram and imaging stuff in different files, but still want to import them as `from exampy.plotting import implot` or `from exampy.plotting import my_cool_hist_func`?
# 
# In this case, we would add the lines 
# `from .histograms import *` and `from .imaging import *` 
# 
# to that inner init (note the dot). By importing `*` we pull all (usually bad practice, you could also specify them by name). Either way, the `from ___ import ____` structure dumps them into this init's namespace without the prepended python module names (histograms and imaging in this case). The dot says to search for those modules in the current directory.

# Alternatively, if in that init, I wrote `from . import histograms` and `from . import imaging`, those modules then get imported as *those names* when the `__init__` runs, meaning the user would have to run, e.g., `exampy.plotting.imaging.implot()` or `exampy.plotting.histograms.my_cool_hist_func.`. Of course, when the USER imports into THEIR code, they can renane those long calls to something more convenient. 
# 
# The behavior you choose between those two is entirely up to you. Going back to Bovy's advice, it comes down to what you think the best balance is between organization/clutter and ease of use. If you have *hundreds* of plotting functions, maybe keeping them more separate is good. If it's only a few in each category, maybe pull them all together in the init so the user calls them directly under `plotting.`. Note the following convenience: if each big function is in its own file and you're using git, changes are *very* easy to track.
# 
# Following along, add one or the other version of the above imports to your `__init__.py`. Then for the *outer* `__init__.py`, do the same, but for `plotting`. In this case, you probably want to use `from . import plotting` rather than `from .plotting import *`. (It's going to know where to look because plotting has its own init file).

# In[16]:


get_ipython().system('more example/exampy/__init__.py')


# In[17]:


get_ipython().system('more example/exampy/plotting/__init__.py')


# We may not have anything interesting to import yet inside our actual `.py` files, but we are actually done setting up the package as far as organization is concerned! 
# 
# Two steps to go: add some code, and make our `setup.py`. Go grab your beautiful `implot()` function from earlier, and paste it into `imaging.py`. Be sure to grab the necessary imports (`numpy, astropy wcs, matplotlib`, etc) for it too, and put those at the top of the `imaging.py` file. 
# 
# 

# We are now ready to make `setup.py`. Create this as a new file in your toplevel directory (for me, `example/`). A very barebones working setup script looks like this:
# 
# ```
# import setuptools
# 
# setuptools.setup(
#     name="exampy",
#     version="0.1",
#     author="Imad Pasha",
#     author_email="imad.pasha@yale.edu",
#     description="A small example Python package",
#     packages=["exampy","exampy/plotting"]
# )
# ```

# (I'll send this snippet to slack). `setuptools` is something that ships with python, and handles most setup things for us. Your `setup.py` for a self-installed, "only for me" package may not look any more complicated than this, but when publishing a package for the public, there's some more bits to add which we won't cover here.

# As a note, where you have to list all packages and subpackages, we could actually have setuptools find these for us via:
# ```
# packages=setuptools.find_packages(include=['exampy','exampy.*'])
# ```
# in the `setuptools.setup()` call, which will search `exampy/` for any submodules. 

# Some other useful things to set include:
# ```
#  python_requires='>=3'
# ```
# to specify the python version be >3. I recommend this, since these days very few maintainers of *new* code bother supporting python 2. There is also 
# 
# ```
# install_requires=["numpy","scipy"]
# ```
# The above argument lets you specify pip-installable packages that are dependancies of your code. If `pip` goes to install your code and these aren't present, it will attempt to install them via `pip` first. As with the python version above, I could specify `numpy>=1.7` to specify versions.
# 

# ## Installing the Package

# At this point, we are ready to install our package. The preferred way to do this these days is with `pip`, which can be used both to retrieve things from PyPI online, or to install things locally. 
# 
# When developing code, it is best to install it in *developer mode*, which makes it easier for your to implement and test changes. This is done via 
# 
# ```
# pip install -e .
# ``` 
# With this, the source is not copied to the installation directory, but rather an entry is made in the installation directory to find the code back in the original directory. This means that any changes you make are immediately available system-wide without requiring a re-installation. Of course, if you have the package already loaded in a Python session, you still have to exit and re-start the session (or use `importlib.reload`). 

# When you are ready to do a full install, all it takes is 
# 
# ```
# pip install .
# ```
# 

# So. I'm going to leave it here, and see if we can get everyone in class to a point where your `implot` function is installed and ready to import!

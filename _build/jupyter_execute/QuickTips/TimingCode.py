#!/usr/bin/env python
# coding: utf-8

# # Timing your code

# Throughout the semester, we will highlight ways to make code faster.  It is almost always faster to stay away from using `for loops`; solutions generally involve using clever options associated with arrays to do quick math. Below are some examples to illustrate this from the first problem in Lab 1.

# In[1]:


import numpy as np
import time


# In[2]:


mean  = 10
sigma = 2


# In Lab 1, Q1 asks you to generate a 1000x1000 array of random numbers.   There are many ways to do this correctly, but some are faster than others.   
# 
# Here are 3 methods for generating the array in Q1.   I will use the command '%%timeit' which runs a cell several times to determine the average run time. Note than when you run a cell using this comand, you'll have to wait ~factor of 10 longer than its normal runtime to give the `timeit` module enough loops to get a steady average.

# ## Method 1

# In[18]:


get_ipython().run_cell_magic('timeit', '', 'arr = np.empty((1000,1000))    \nfor x in range (1000):\n    val = np.random.normal(mean,sigma,1000)\n    arr[x] = val\n    ')


# ## Method 2

# In[14]:


get_ipython().run_cell_magic('timeit', '', 'array = np.array([[np.random.normal(mean, sigma) for i in range(1000)] for j in range(1000)])')


# ## Method 3

# In[15]:


get_ipython().run_cell_magic('timeit', '', '\narray = np.random.normal(mean, sigma, size = (1000,1000))')


# Method 1 looks like a `for loop`, yet it is only marginally slower than Method 3 which doesn't use a `for loop`.   Method 2 looks like it could be fast, but the double for loop makes it 100 times slower!!
# 
# One reason for this is that Method one loops over only one axis (rows, in this case), and creating a 1000 length long random array is very fast, even if you do it 1000 times. But Method 2, which loops over every pixel, has to run the generator 1000x1000 = 1,000,000 times! 
# 
# For pieces of code that are hard to isolate, you can also determine run time using `import time`.   Below I use both to explicitly show what %%timeit is doing.   The cell below is run many times, and you can see the different runtimes for the same command.

# In[20]:


get_ipython().run_cell_magic('timeit', '', "t0 = time.time()\narray = np.random.normal(mean, sigma, size = (1000,1000))\nt1 = time.time()\nprint('Time = {:0.5f}'.format(t1-t0))")


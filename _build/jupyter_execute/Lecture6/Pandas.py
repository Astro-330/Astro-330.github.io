#!/usr/bin/env python
# coding: utf-8

# # 6. Pandas

# What is `pandas`?
# 

# - Lot of answers, depends on your use case. 
# - Primarily, pandas is a package built around the DataFrame
# - DataFrames are a more flexible, SQL-like data container

# ## Numpy vs Pandas
# 
# When would we want to use `pandas` instead of a trusty `numpy` array?

# - When we need to index columns by name
# - When we need to index *rows* by name 
# - When we want to store heterogeneous data (e.g., one column of strings, next column of floats)
# - When we have multiple tables we want to join together
# - When we want to perform complex filtering operations on a single table. 
# - When the I/O for some data on disk is frustrating using `numpy`

# Note that the first bullet point (indexing columns of an array by name) can be accomplished via some other data containers, including `Astropy` `Table`s and `numpy` `record arrays`. 
# 
# *However*, many of the other bullet points are not easily accomplished with these methods. *Furthermore*, `pandas` is the industry standard in data science and other fields. Learning it is a transferrable skill.

# ## Introduction to Pandas

# For the purposes of this introduction, we need to learn about 2 `pandas` objects: 
# - `pd.Series`, and 
# - `pd.DataFrame`. 
# 
# DataFrames are constructed out of Series objects, but as it turns out, we don't need to define them explicitly on most occasions. 

# Let's set up a data frame:

# In[2]:


import numpy as np 
import pandas as pd 

array1 = np.array([1,2,3])
array2 = ['Good','Bad','Ugly']

df = pd.DataFrame({'measurement':array1,'quality':array2})
df


# As you can see, one valid initialization is feeding in a dictionary with the desired column names as keys and the desired columns as values.

# Here's another way:

# In[3]:


array1 = np.array([1,2,3])
array2 = ['Good','Bad','Ugly']

df = pd.DataFrame()
df['measurement'] = array1
df['quality'] = array2
df


# i.e., we can create empty `DataFrames`, and dynamically set new columns the way we are used to setting the values of dictionaries for a given key. 

# It is worth noting that this "dictionary-style" indexing of DataFrame columns is not the only way. If our column names have no spaces, we can also use dot notation:

# In[4]:


df['measurement']


# In[5]:


df.measurement


# As you can see, both forms return the same array. Notice, though, that this "column" isn't just an array... here we see our first pandas `Series`. 

# You can think of a `Series` as a sort of wrapper/container around data arrays that provides extra functionality. 

# In[6]:


df.quality


# For example, we can see that this Series object has a `Name` (which it tells us when we print it), and it has an `index` --- that's the 0,1,2 over on the left. 
# As a note, if we ever want to extract *just* the data values of a series, we run

# In[7]:


df.quality.values


# I've found that the dot-notation indexing of DataFrames can be very handy -- `.thing` is fewer characters to type than `['thing']`, and doesn't require your fingers to find the `[]` keys. 
# 
# Let's come back to our DataFrame, `df`. What's going on with this index on the left? 

# In[8]:


df


# In[10]:


df.loc[0]


# The **index** part of a DataFrame is a way to index specific *rows*, rather than specific *columns*. As you can see above, I used the `.loc[]` syntax, which tells pandas I'm inserting a row index, rather than a "column index". Pandas then spits out each column and associated value for that row.

# We can set the index of a DataFrame to be whatever we want. For example, let's say our measurements were taken on 3 dates:

# In[11]:


df.index = ['2021-10-07','2021-10-08','2021-10-09']


# In[12]:


df


# Now the `[0,1,2]` is replaced by our date strings. 

# Now, if I want to pull the measurement and measurement quality from `2021-10-08`, I would use

# In[13]:


df.loc['2021-10-08']


# What if I know the numerical index of the row I want, but don't know the name of that row? Pandas doesn't *love* it when we are in this state, but there is a solution:

# In[14]:


df.iloc[1]


# The `df.iloc[]` syntax allows you to index a row by true index number, rather than whatever the name of that row is. 

# ## Filtering DataFrames

# So now that we know how to create a simple DataFrame, how do we go about filtering it? Instead of creating external masks or using `np.where()`, `pandas` has a built in set of methods for accessing parts of a dataframe.

# Sticking with our measurement df, what if I want a view of my DataFrame, but only for the rows where the `quality` is either `'Bad'` or `'Ugly'`. We can do that as follows:

# In[17]:


df.loc[(df.quality=='Bad')|(df.quality=='Ugly')]


# We can see that by putting our conditions inside `df.loc[]` the way we normally would in `np.where()`, we return a version of the DataFrame where only the bad and ugly measurement remains.

# What if I want to perform the same filtering, but I don't want the whole DataFrame out, I just want a specific column? We can specify this as well:

# In[18]:


df.loc[(df.quality=='Bad')|(df.quality=='Ugly'),'measurement']


# We add a comma after our conditionals, and then specify the name of the column we want.

# This is extensible to multiple columns; we only have two here but it would look like this:

# In[19]:


df.loc[(df.quality=='Bad')|(df.quality=='Ugly'),['measurement','quality']]


# i.e., we provide a list of column names we want to index, only where the conditions are met. 

# Let's say I want to set any 'Bad' or 'Ugly' measurements to `np.nan`:

# In[22]:


df.loc[(df.quality=='Bad')|(df.quality=='Ugly'),'measurement'] = np.nan
df


# **Warning**. A common "gotcha" with DataFrames is to not use the `col_indexer` as above when setting values, but rather something like the following:

# In[23]:


df.loc[(df.quality=='Bad')|(df.quality=='Ugly')].measurement = np.nan


# We get this error because running `.measurement` or indexing with `['measurement']` creates a copy of that slice of the dataframe for view -- meaning that if we try to set values in it, we are **not** affecting values in the underlying dataframe.

# ## Joining DataFrames

# A powerful strength of DataFrames is the ability to join/merge multiple dataframes on a given column. This behavior is inhereted from the way Databases (such as SQL) work. 
# Let's create another DataFrame that has some anciliary data also taken on our three dates:

# In[25]:


df2 = pd.DataFrame()
df2['weather'] = ['sunny','cloudy','cloudy']
df2['seeing'] = np.array([0.55,2.3,3.2])
df2.index = ['2021-10-07','2021-10-08','2021-10-09']
df2


# We can easily combine these two DataFrames:

# In[31]:


new_df = pd.merge(df,df2,left_index=True,right_index=True)
new_df


# When merging dataframes, we have to select which columns to merge on. Here, I've chosen a full join on the index. Let's talk about the options:

# From database parlance, we can join/merge two tables as long as there is a *unique identifer* column that is shared across both tables. In our case, the *index* is that "column". 
# 
# - A left (outer) join means that we take *all* rows of `df`, and if df2 has no entry for some rows, we'll fill it with NaNs. 
# - A right (outer) join means we will keep *all* rows of `df2`, and if df has no entry for some rows, we'll fill it with NaNs.
# - A full (outer) join means we keep *all* rows of **both** dataframes, filling in NaNs where either is missing data. 
# 
# Alternatively, we could ask for just the rows that have matches in both dataframes, dropping all other columns. This would be an inner join.

# Let's see this in action. I'm going to add a date, measurement, and quality to df:

# In[34]:


df.loc['2021-10-10'] = [5,'Good']
df


# Now when we go to do our merge, `df2` won't have a row for that last index. But if I use `df.join()`, I can specify an outer join:

# In[40]:


new_df = df.join(df2,how='outer')
new_df


# Notice that I get NaN for the weather and seeing, as these weren't present in `df2`. But I do keep that row overall. 

# In[41]:


new_df = df.join(df2,how='inner')
new_df


# Alternatively, if I use an "inner" join, I only keep indices that were in *both* dataframes. So I lose that entire row. 

# In[42]:


new_df = df.join(df2,how='right')
new_df


# This time, I've specified I want a "right outer" join. Because the "righthand" dataframe (`df2`) doesn't have a row for '2021-10-10', it's not part of the final df. 

# In[43]:


new_df = df.join(df2,how='left')
new_df


# Wheras if I use 'left', I keep all rows of `df` including the new one, filling in NaNs for df2 where it has no values. 

# You may have noticed that in the above examples, 'left' and 'outer' produced the same results. If so, good eye! They're the same output in this case because `df2` has no rows *not* present in `df`. But if it did, then the 'outer' join would keep both that row and '2021-10-10' from df, while 'left' would ignore the extra row in df2.

# As a note, the join's we've been doing have *implicitly been using the index* as the join column, which is the default behavior. But if there were an actual column in the data that had the necessary properties, we could use the `pd.merge(df,df2,on='name')` syntax to join on that column instead.

# ## Operations on Columns

# We can carry out operations on columns of dataframes like we would with numpy arrays. For example, let's say I have a dataframe as follows:

# In[67]:


df3 = pd.DataFrame()
df3['galaxy_name'] = ['NGC 4300','NGC 3055','NGC 3235','NGC 6532']
df3['flux'] = np.array([1.1e-12,1.5e-13,4.6e-13,1.8e-12])
df3['distance'] = np.array([10.5,33.7,105,22])
df3


# I can calculate the luminosity of these galaxies and make a new column as follows

# In[68]:


import astropy.units as u
lums = (df3.flux.values*u.erg/u.s/u.cm**2* 4 * np.pi * (df3.distance.values*u.Mpc)**2).to(u.Lsun).value
df3['luminosity'] = np.log10(lums)
df3


# In this case, I use the `.values` to be sure things play nicely with `astropy`. 

# ### Dropping rows or columns 
# 
# So far, we've seen how to add things into dataframes and modify their values. But what if we want to drop a row or column entirely?
# 
# For this, we need `df.drop()`. This command returns a version of the dataframe with the row or column dropped. As a note, the default behavior is to drop *rows*, so if you're specifying a column, you also need to add `axis=1` as follows:

# In[76]:


df2.drop('2021-10-09')


# We specified an index there, let's try a column:

# In[77]:


df2.drop('seeing',axis=1)


# Note that `df2.drop()` doesn't overwrite df2. To do this, we'd need do do either of the following:

# In[ ]:


df2 = df2.drop('seeing',axis=1)

#OR 

df2.drop('seeing',inplace=True,axis=1)


# I suggest the first, as it is easier to read and remember.

# ## Other Common Operations

# Here I'm going to provide a grab bag of common pandas methods that are useful to know:

# We can return the unique elements of a column (and the number of them) via 
# `df['col2'].unique()`
# and 
# `df['col2].nunique()`
# 
# We can count how many times each value appears in a given column via 
# `df['col2'].value_counts()` 
# 
# We can sort a dataframe by the values in some column via
# `df.sort_values('col2')`
# 

# ## Data I/O

# A final use for `pandas` is that it was designed to work with *messy* data. When loading data from, e.g., ascii files, there are often a lot of subtleties surrounding the spacings and columns, mixes of strings and floats, etc. 
# 
# Traditional `numpy` functions like `loadtxt()` or `genfromtxt()` struggle with these kinds of files, but often, `pandas` has no trouble with it. At the very least, then, we can use pandas to get the data into python, where we know how to mess with it until it meets our requirements.

# In[69]:


df = pd.read_csv('table2.dat')


# In[71]:


df.head()


# Above, I've run a one line, no modification call of pandas' `read_csv()` function. With no issue, it loaded up this file, which contains column names (first row), and a variety of data types. 

# Let's try the same thing with numpy:

# In[72]:


x = np.loadtxt('table2.dat')


# Yep. Numpy doesn't like it. In part, having columns with strings and also columns with floats leaves us dead in the water. We want pandas for this.
# 
# There are cases where pandas is annoying about reading things in and numpy is actually easier. My advice? Know both, and pick whichever one to read in the data causes the least headache. You can always turn an array into a dataframe (and vice versa, if the data are all the same type).

# We're going to use the following dataset in the lab, so I'd like to briefly familiarize you with it:

# In[74]:


df.head()


# In[75]:


len(df)


# In[94]:


df.groupby(['name','ion','wl']).flux.describe()


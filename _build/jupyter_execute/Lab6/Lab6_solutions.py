#!/usr/bin/env python
# coding: utf-8

# # Lab 6: Solutions

# In this lab, we'll work through some of the basics of using Pandas, using a few different tabular data sets. Ultimately, one need not do anything particularly fancy with DataFrames for them to be useful as data containers. But we would like to highlight a few extra abilities these objects have, that illustrate situations where we may actually have a strong reason to use pandas over another library. 

# ## Problem 1: HII regions + Planetary Nebulae measurements in M81
# 
# For our first data set, we're going to look at a file (`table2.dat`), which contains measurements of the flux and intensity of various ions' line emission from a set of known emitting objects (PNs and HII regions) in the M81 galaxy. 
# 
# The columns of this file are `name`, `ion`, `wl`, `flux`, and `I` (intensity). Two of the columns are string-valued (name and ion), three are numerical-values (wl, flux, I). This mix of strings and floats tells us before we even decide how to read in this file that `numpy` data structures won't be usable, as they demand all values in an array to have the same `dtype`. 
# 
# ### Problem 1.1 
# 
# Using the `pd.read_csv()` function shown in the lecture, read this data file into a dataframe called `df`, and print it. 
# ```{hint}
# You can get a "pretty" visualization of a dataframe by simply typing its name into a jupyter cell -- as long as it's the last line of the cell, the dataframe will print more nicely than typing `print(df)`. This does not work outside of notebooks.
# ```

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
df = pd.read_csv('table2.dat')


# In[5]:


df


# ### Problem 1.2 
# 
# Though it doesn't show up in the clean representation above, the strings associated with the name and ion columns above have trailing and leading spaces that we don't want. 
# 
# Use a *list comprehension* to modify the data frame such that each value in the name and ion columns are replaced with a `.strip()` version of themselves.

# In[12]:


df['name'] = [i.strip() for i in df.name]
df['ion'] = [i.strip() for i in df.ion]


# ### Problem 1.3 
# 
# Write a function `select_object` which takes in as an argument the name of an HII region or planetrary nebula, and filters the dataframe for only the entries for that object using `df.loc[]`. Consider having the dataframe be an optional argument you set to `df`, the dataframe we are 
# working with.
# 
# Have your function take in an optional argument `drop_empty=True` which additionally selects only those rows where the flux/intensity is **not** zero.

# In[14]:


def select_object(name,df=df,drop_empty=True):
    if drop_empty:
        out = df.loc[(df.name==name)&(df.flux!=0)&(df.I!=0)]
    else:
        out = df.loc[(df.name==name)]
    return out


# In[15]:


select_object('PN3m')


# ### Problem 1.4 
# 
# Write a function `select_ion_by_wavelength()` which takes in the name of an ion and its wavelength (and a dataframe), and returns the filtered data frame for all objects, but only ions for the selected wavelengths. 
# 
# As before, have a `drop_empty` optional argument to not include entries where the flux and intensity are zero.
# 
# Additionally, as the index is now uniquely identified by the name of the PN/HII region, set the index to be the name column.

# In[26]:


def select_ion_by_wavelength(ion,wavelength,df=df,drop_empty=True): 
    if drop_empty:
        out = df.loc[(df.ion==ion)&(df.wl==wavelength)&(df.flux!=0)&(df.I!=0)]
    else:
        out = df.loc[(df.ion==ion)&(df.wl==wavelength)]
    out = out.set_index('name')
    return out


# In[27]:


select_ion_by_wavelength('[OII]',3727).head()


# ### Problem 1.5
# It will be helpful to know for a given ion which wavelengths are avalable and have data in the dataframe. Write a function `get_wavelenghs_by_ion()` that determines for a given input ion, which wavelengths are available. 
# 
# **Bonus + 0.5: Some ions of forbidden transitions like `[OII]` have brackets in the name. Add a bit to your get_wavelengths_by_ion code that allows the user to enter either `"[OII]"` or `"OII"` and get the same answer.**
# 
# Additionally, make a convenience function `get_ions()` that just returns the full list of ions represented in the dataframe.
# 
# ```{hint}
# The `.unique()` method in pandas will be useful here.
# ```

# In[21]:


def get_ions(df=df):
    return df.ion.unique()

def get_wavelengths_by_ion(ion,df=df,drop_empty=True):
    return df.loc[df.ion==ion,'wl'].unique()


# In[22]:


get_ions()


# In[24]:


get_wavelengths_by_ion('[NII]')


# In[25]:


get_wavelengths_by_ion('[OII]')


# ### Problem 1.6 
# 
# Rather than having all these convenience functions littered around our code, let's go ahead and make a class, `FluxTable`, which initializes with our dataframe, and then has all of the functions created above as methods. The input DataFrame, `df`, should be accessible as an attribute as well. 
# 
# When you're done, you should be able to do something like the following 
# ```
# ions = FluxTable(df)
# print(ions.df)
# ions.get_ions()
# ions.get_wavelengths_by_ions('[OII]')
# PN3m = ions.select_object('PN3m')
# ```

# In[29]:


class FluxTable():
    def __init__(self,df):
        self.df=df
    def select_object(self,name,drop_empty=True):
        df=self.df
        if drop_empty:
            out = df.loc[(df.name==name)&(df.flux!=0)&(df.I!=0)]
        else:
            out = df.loc[(df.name==name)]
        return out
    def select_ion_by_wavelength(self,ion,wavelength,drop_empty=True): 
        df=self.df
        if drop_empty:
            out = df.loc[(df.ion==ion)&(df.wl==wavelength)&(df.flux!=0)&(df.I!=0)]
        else:
            out = df.loc[(df.ion==ion)&(df.wl==wavelength)]
        out = out.set_index('name')
        return out
    def get_ions(self):
        df=self.df
        return df.ion.unique()

    def get_wavelengths_by_ion(self,ion,drop_empty=True):
        df=self.df
        return df.loc[df.ion==ion,'wl'].unique()


# In[30]:


ft = FluxTable(df)


# In[31]:


ft.get_ions()


# In[32]:


ft.get_wavelengths_by_ion('HeI')


# In[33]:


ft.df.head()


# In[34]:


ft.select_object('PN3m')


# ### Bonus: (+1) Problem 1.7 
# 
# Finally, let's add one final method to our class. This method will be called `query`, and it will act as a bit of a "catch all", allowing the user to query for certain conditions as desired. 
# 
# Your `query` method should take as it's primary argument a string containing a comma separated list of desired columns. It should then have optional arguments for `name`, `ion`, and `wl`, which are by default set to `None`. For name and ion, the goal is to allow the user to specify specific ones. For `wl`, we'll go one step further and allow either a specific wavelength, or a range of wavelengths input as a string of the form `>4343` or `<3050` or `2010-5000`. 
# 
# The usage of this method will look something like 
# 
# ```
# ft.query('name,flux',ion='[OII]',wl='3000-5000')
# ```
# 
# You will of course need to do some string checking (particularly with wl) to figure out what the user wants, and then you can use your filtering methods you already wrote to successfully construct a result dataframe to return.

# In[ ]:





# ## Problem 2
# 
# In this problem, we're going to use the [3DHST catalog](https://archive.stsci.edu/prepds/3d-hst/), which contains numerous measurements of galaxies at high redshift taken with HST. This program was led at Yale(!) and the "3D" refers to the fact that beyond just imaging, spectroscopic information was also obtained.  
# 
# We'll be using a subset of the catalog set up for the GOODS-South field, a well studied patch of the sky. The data are split across four `fits` files -- but we'll be using `pandas` to join them together! 
# 
# - the `.cat` file contains the primary catalog 
# - the `.fout` file contains the output of `FAST`, a quick template fitting program for galaxies. 
# - the `.RF` file contains the Rest Frame (de-redshifted) colors of the galaxies in common bands
# - the `.zout` file contains redshift estimates for all galaxies (either spec-z or photo-z) made using the EAZY redshift fitting code (also Yale!) 
# 
# ### Problem 2.1
# 
# **Load the four datasets into Python, and create dataframes for each.** For ease of following the solutions, we suggest you name these
# 
# - `cat_df` for the catalog
# - `fast_df` for the fast output
# - `rf_df` for the RF file 
# - `z_df` for the redshifts
# 
# **Examine each of these dataframes to see what types of columns they have.**
# 
# ```{hint}
# Remember that the default extension for tabular data (as this is) will be 1, not 0 as we are used to for images. You can run pd.DataFrame() directly on the data attribute of this extension
# ```

# In[2]:


from astropy.io import fits


# In[3]:


with fits.open('goodss_3dhst.v4.1.cat.FITS') as hdu:
    cat_df=pd.DataFrame(np.array(hdu[1].data).byteswap().newbyteorder()) 


# In[6]:


cat_df


# In[4]:


with fits.open('goodss_3dhst.v4.1.fout.FITS') as hdu:
    fast_df = pd.DataFrame(np.array(hdu[1].data).byteswap().newbyteorder()) 


# In[8]:


fast_df


# In[5]:


with fits.open('goodss_3dhst.v4.1.master.RF.FITS') as hdu:
    rf_df = pd.DataFrame(np.array(hdu[1].data).byteswap().newbyteorder()) 


# In[10]:


rf_df


# In[6]:


with fits.open('goodss_3dhst.v4.1.zout.FITS') as hdu:
    z_df = pd.DataFrame(np.array(hdu[1].data).byteswap().newbyteorder()) 


# In[12]:


z_df


# ### Problem 2.2 
# 
# You should notice that every one of these tables has 50507 rows. This is a relief! It means the creators were consistent, and that each object has a row in each table. 
# 
# You should also notice one column, `id`, is a unique identifier for each row in each table (and is consistent across tables). Pandas has assigned its own index (it's off by 1 from the id column), but we might as well use `id`. 
# 
# **For each of your four dataframes, set 'id' to be the index column.**

# In[7]:


cat_df = cat_df.set_index('id')
fast_df = fast_df.set_index('id')
rf_df = rf_df.set_index('id')
z_df = z_df.set_index('id')


# In[14]:


rf_df.head()


# ### Problem 2.3 
# 
# Instead of working with these four dataframes separately, let's join them all into one MEGA DATAFRAME. 
# 
# By setting 'id' as the index for each, we'll be able to merge the dataframes using `pd.merge` with the `left_index` and `right_index` parameters set to true. 
# 
# **Create your mega dataframe, which we'll simply call `df`, by consecutively merging the first two, then third and fourth dataframes.**
# 
# ```{hint}
# You should end up with 215 columns in your final dataframe.
# ```

# In[8]:


df = cat_df.merge(fast_df,how='outer',left_index=True,right_index=True)
df = df.merge(rf_df,how='outer',left_index=True,right_index=True)
df = df.merge(z_df,how='outer',left_index=True,right_index=True)


# In[9]:


tmp_df  = pd.merge(cat_df, fast_df, left_index=True,right_index=True)
tmp2_df = pd.merge(tmp_df, rf_df,   left_index=True,right_index=True)
mega_df = pd.merge(tmp2_df, z_df,   left_index=True,right_index=True)


# In[18]:


mega_df


# ### Problem 2.4 
# 
# Let's take a look at the redshift distribution of sources in the catalog.
# 
# There are several estimates of the photometric redshift, the one we want to use is `z_peak` for photometry, or, if available, `z_spec_x`, from spectroscopy. 
# 
# **What percentage of the catalog have measured spectroscopic redshifts? The `z_spec_x` column is set to `-1` if no spectroscopic redshift exits.**

# In[92]:


(len(df.loc[df.z_spec_x!=-1,'z_spec_x']) / len(df))*100


# ### Problem 2.5 
# 
# Write a function `get_redshift()` which takes in an object ID, and returns `z_spec_x` for that source if it's not -1, otherwise returns `z_peak`. Because `id` is a special word in python, we suggest using `objid` as the input, and setting df=df as an optional argument. You can make this a memory lite function by using `df.loc[]` to pull the row and only the two columns you need for this.
# 
# There are two additional "flagged" values: -99.0 and -99.9 -- Have your function output np.nan if this is the value in the table. 
# 
# Your function should return the redshift as well as a flag (string) 's' or 'p', or 'f' for spectroscopic or photometric (or fail, if you're returning nan). 

# In[10]:


def get_redshift(objid,df=df):
    options = df.loc[objid,['z_spec_x','z_peak']]
    if options.z_spec_x != -1:
        if options.z_spec_x != -99.0 and options.z_spec_x != -99.9:
            return options.z_spec_x, 's'
        else:
            return np.nan,'f'
    else:
        if options.z_peak != -99.0 and  options.z_peak != -99.9:
            return options.z_peak, 'p'
        else:
            return np.nan, 'f'


# In[11]:


get_redshift(150)


# ### Problem 2.6 
# 
# Now that we can get the best redshift for each row, use a list comprehension to grab these values for every object. You can index the output tuple of your function at 0 to drop the flag for now. 
# 
# Once you have this, plot a histogram of the redshifts, using `fig, ax = plt.subplots`. Make your plot nice!
# 
# ```{note}
# My list comprehension takes ~15 seconds to run. It's a lot of data! If you wish, you may try to find a more optimized solution built on a full column operation rather than a loop. One possibility is to take the spec-z column, mask any bad values, and then replace those entries in the z-phot column and plot that... 
# ```

# In[109]:


all_z = [get_redshift(i)[0] for i in df.index]


# In[110]:


all_z = np.array(all_z)


# In[113]:


fig, ax = plt.subplots(figsize=(9,9))
ax.hist(all_z[~np.isnan(all_z)],bins=50)
ax.set_xlabel('redshift',fontsize=15)
ax.set_ylabel('Number of galaxies',fontsize=15)
ax.tick_params(direction='in',top=True,right=True,length=8);


# ### Problem 2.7
# 
# Now do the same, but separately plot the distributions of redshift for those with spectroscopic redshifts and those that only have photometric. For this, you'll want to set `density=True` in your `hist`, and play with the linestyles and opacities so that both are visible. 
# 
# **Bonus (+0.5): Use KDE from seaborn to smoothly represent the distribution**. 

# In[114]:


spec_z = [get_redshift(i)[0] for i in df.index if get_redshift(i)[1]=='s']; spec_z = np.array(spec_z)


# In[115]:


phot_z = [get_redshift(i)[0] for i in df.index if get_redshift(i)[1]=='p']; phot_z = np.array(phot_z)


# In[143]:


fig, ax = plt.subplots(figsize=(9,9))
ax.hist(spec_z,bins=40,density=True,alpha=0.9,label='Spectroscopic')
ax.hist(spec_z,bins=40,density=True,alpha=0.5,histtype='step',color='k',lw=3)

ax.hist(phot_z,bins=50,density=True,alpha=0.65,histtype='stepfilled',label='Photometric')
ax.hist(phot_z,bins=50,density=True,alpha=0.5,histtype='step',color='k',lw=3)

ax.set_xlabel('redshift',fontsize=15)
ax.set_ylabel('Number of galaxies',fontsize=15)
ax.legend(prop={'size': 16})
ax.tick_params(direction='in',top=True,right=True,length=8);


# **Do the differences between the two distributions make sense? Why or why not?**

# ```{note}
# Yes! We can see that the spectroscopic redshifts are clustered toward low redshift, while the tail of the photometric redshifts extends to higher redshift. This is because obtaining a spectrum is many times more observationally expensive than obtaining photometry, especially for many objects, and objects far enough away become too faint to measure in a survey of a set depth. 
# ```

# ## Problem 3
# 
# The "UVJ diagram" is a useful diagnostic tool in extragalactic astronomy which allows one to relatively accurately diagnose a galaxy as star forming or quiescent, even in the presence of dust. It is composed by plotting the "U-V" color of the galaxy against the "V-J" colors. You'll likely know U and V (these are from the Johnsons-Cousin's filter system). J is a filter in the near infrared. 
# 
# In this problem, we're going to write a function that can create a UVJ diagram for subsets of our data, cutting on mass and redshift. 
# 
# You'll need to access the following columns in the data (besides redshift, which you've already handled):
# 
# - stellar mass: the mass of all the stars in the galaxy. (column: `lmass`, flagged value of -1)
# - star formation rate: rate at which galaxy is forming stars. (column: `lsfr`, flagged value of -99.0)
# - U band flux density (column: `l153`, flagged value of -99.0)
# - V band flux density (column: `l155`, flagged value of -99.0)
# - J band flux density (column: `l161`, flagged value of -99.0)
# 
# ### Problem 3.1 
# 
# For step one, we need to be able to filter our dataframe for particular mass and redshift ranges. 
# 
# Write a function, `select_galaxies()`, which takes as arguments `M, delta_M, z, delta_z` (and our dataframe). 
# It should then return a df of only systems between M-deltaM to M+deltaM, and z-deltaz to z+deltaz. The columns it should return are the ones specified above.
# 
# There is actually a column in `rf_df` called `z`, that contains the spec_z if available or the peak z if not. At the time of writing, I cannot determine why this column was not included in the merge. In any case, set `df['z']` equal to `rf_df.z` before continuing, as you'll use it below.
# 
# 
# ```{note}
# All masses and sfrs are in log units. 
# ```
# 
# Try your function out using a mass of 10, delta M of 0.5 (i.e., a bin from 9.5 - 10.5), a redshift of 1, and a delta z of 0.25.

# In[15]:


df['z'] = rf_df['z']


# In[12]:


def select_galaxies(M,delta_M,z,delta_z,df=df):
    Mmin = M-delta_M
    Mmax = M+delta_M
    zmin = z-delta_z
    zmax = z+delta_z
    return df.loc[((df.lmass>Mmin)&(df.lmass<Mmax)&(df.z>zmin)&(df.z<zmax)),['lmass','lsfr','l153','l155','l161']]


# In[215]:


select_galaxies(10,0.5,1,0.25)


# ### Problem 3.2
# 
# Great, we can now get subsamples in mass/redshift bins. This is important because the UVJ diagram actually changes as a function of mass and redshift. 
# 
# Next we need to get the colors out. Write a function `get_colors()` which takes the same arguments as your `select_galaxies()` function above. Inside, it should run `select_galaxies` passing through the arguments, and then from the resulting data frame, calculate U-V and U-J (see below). Add these as columns to said dataframe with names 'U-V' and 'V-J' and return it. 
# 
# Run this function with the same mass/redshift bin from above and look at it.
# 
# ```{warning}
# As noted above, the U,V, and J band data are in Fnu (flux densities). Thus, a color is computed via -2.5*log10(Lfilter1/Lfilter2)
# ```
# 

# In[13]:


def get_colors(M,delta_M,z,delta_z,df=df):
    minidf = select_galaxies(M,delta_M,z,delta_z,df=df)
    minidf['U-V'] = -2.5*np.log10(minidf.l153/minidf.l155)
    minidf['V-J'] = -2.5*np.log10(minidf.l155/minidf.l161)
    return minidf


# In[219]:


get_colors(10,0.5,0.75,0.25)


# ### Problem 3.3 
# 
# Now that we can easily grab U-V and V-J colors for subsamples, we're ready to plot! 
# 
# Next, set your xlim and ylim to (0,2) in V-J (x axis) and (0,2.8) in U-V (y axis).
# 
# Once you have the distribution plotted nicely, use the definitions of the bounding box provided in [Whitaker et al. 2011](https://iopscience.iop.org/article/10.1088/0004-637X/735/2/86/pdf) (Eqns 15, Fig 17 for example) to draw the box where quiescent galaxies sit. 

# In[28]:


fig, ax = plt.subplots(figsize=(9,9))
ax.set_box_aspect(1)
ob = get_colors(9.5,2,1.0,0.5)
ax.plot(ob['V-J'],ob['U-V'],'.',color='k',alpha=0.1);
ax.plot([0,0.7],[1.3,1.3],'r',lw=3)
ax.plot([1.6,1.6],[2.0,2.8],'r',lw=3)
ax.plot([0.7,1.6],[1.3,2.0],'r',lw=3)
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=20)
ax.set_xlim(0,2)
ax.set_ylim(0,2.8)
ax.set_xlabel('V - J color',fontsize=20)
ax.set_ylabel('U - V color',fontsize=20)
ax.tick_params(direction )


# ### Bonus (+2) 
# 
# (+1) Now that you can easily plot up a UVJ diagram for a given mass/redshift bin, make a plot with 6 panels. Across the 3, plot the UVJ diagram in the redshift bins 0-0.5, 0.5-1.0, 1.0-2.0. In the top panel, include galaxies in the mass range 8-9.5, and in the bottom, 9.5-11. 
# 
# (+1) Feeling extra fancy? Use the conditions on the UVJ quiescent box to color the quiescent galaxies red and the star forming galaxies blue in you plots. 
# 
# 

# In[ ]:





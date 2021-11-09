import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from spectral_cube import SpectralCube
import astropy.units as u 
import copy 

st.title('FITS Cube Viewer')


if 'channel' not in st.session_state:
    st.session_state['channel'] = 0


@st.cache(allow_output_mutation=True)
def load_data():
    cube = SpectralCube.read('ngc1333_12co.fits')
    cube = cube.with_spectral_unit(u.km / u.s)  
    return cube   

data = load_data() 
st.session_state.channel = st.slider('Choose Channel',min_value=0,max_value=len(data),value=None,step=1)
fig, ax = plt.subplots(figsize=(7,7))
ax.imshow(data[st.session_state.channel].data) 
st.pyplot(fig)

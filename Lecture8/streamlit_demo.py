import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from spectral_cube import SpectralCube
import astropy.units as u 
import copy 


st.title('FITS Cube Viewer')

st.header('Channel Viewer')
if 'channel' not in st.session_state:
    st.session_state['channel'] = 0

if 'oldchannel' not in st.session_state:
    st.session_state['oldchannel'] = 0 

if 'velocity' not in st.session_state:
    st.session_state['velocity'] = 0


@st.cache(allow_output_mutation=True)
def load_data():
    cube = SpectralCube.read('ngc1333_12co.fits')
    cube = cube.with_spectral_unit(u.km / u.s)  
    return cube   

data = load_data() 

if 'momrange' not in st.session_state:
    st.session_state['momrange'] = (0,len(data)) 

col1,col2,col3 = st.columns(3)
if col1.button('Previous Channel'):
    if st.session_state.channel-1<=0:
        st.write('Reached channel 0')
    else:    
        
        st.session_state.channel-=1
        st.session_state.oldchannel = copy.copy(st.session_state.channel)

        
if col3.button('Next Channel'):
    if st.session_state.channel + 1 >= len(data):
        st.write('Reached last channel!')
    else:
        
        st.session_state.channel+=1
        st.session_state.oldchannel = copy.copy(st.session_state.channel)


st.session_state.channel = st.slider('Choose Channel',min_value=0,max_value=len(data),value=st.session_state.oldchannel,step=1)

fig, ax = plt.subplots(figsize=(7,7))
ax.imshow(data[st.session_state.channel].data) 
st.pyplot(fig)
st.session_state.velocity = float(data.spectral_axis[st.session_state.channel].value)
col2.metric('Velocity [km/s]', f'{st.session_state.velocity:.2f}', delta=None, delta_color="normal")

st.header('Moment Map Creator')

st.session_state.momrange = st.slider('Select a range of channels',130, 145, (0, len(data)))

fig2, ax2 = plt.subplots(figsize=(7,7))
img = np.mean(data.unmasked_data[st.session_state.momrange[0]:st.session_state.momrange[1]],axis=0)
ax2.imshow(img) 
st.pyplot(fig2)

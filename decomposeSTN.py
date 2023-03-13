#!/usr/bin/env python

# Leonardo Fierro 
# Acoustics Lab, Aalto University, 2021
# Last modified: 14.09.2022 - L. Fierro 
# -------------------------------------------------------------------------
import numpy as np
from scipy import signal as sig
from scipy import ndimage
import numpy as np

def decSTN(x,Fs,nWin1,nWin2):

    nHop1 = nWin1*7//8
    nHop2 = nWin2*7//8
    win1 = sig.windows.hann(nWin1,sym=False)
    win2 = sig.windows.hann(nWin2,sym=False)
    filter_length_t = 200e-3 # ms
    filter_length_f = 500 # Hz

    # Round 1
    f,t,X1 = sig.stft(x, fs=Fs, window=win1, nperseg=nWin1, noverlap=nHop1, return_onesided=True)
    
    nH_1 = int(np.round(filter_length_t * Fs // (nWin1-nHop1)))
    nV_1 = int(np.round(filter_length_f * nWin1 // Fs))
    S1,T1,N1 = fuzzySTN(X1,0.7,0.8,nH_1,nV_1)

    t,xs = sig.istft(S1*X1, fs = Fs, window=win1, nperseg=nWin1, noverlap=nHop1, input_onesided=True)
    t,xres = sig.istft((T1+N1) * X1, fs = Fs, window=win1, nperseg=nWin1, noverlap=nHop1, input_onesided=True)
    
    # Round 2
    
    f,t,X2 = sig.stft(xres, fs=Fs, window=win2, nperseg=nWin2, noverlap=nHop2, return_onesided=True)
    nH_2 = int(np.round(filter_length_t * Fs // (nWin2-nHop2)))
    nV_2 = int(np.round(filter_length_f * nWin2 // Fs))
    
    S2,T2,N2 = fuzzySTN(X2,0.75,0.85,nH_2,nV_2)

    t,xt = sig.istft(T2 * X2, fs = Fs,window=win2, nperseg=nWin2, noverlap=nHop2, input_onesided=True)
    t,xn = sig.istft((S2+N2) * X2, fs = Fs, window=win2, nperseg=nWin2, noverlap=nHop2, input_onesided=True)   

    return xs,xt,xn

# Extract STN masks from transientness of input STFT
def fuzzySTN(X,G2,G1,nMedianH,nMedianV):

    X_h_median = ndimage.median_filter(np.absolute(X),size=(1,nMedianH+1),mode='constant')
    X_v_median = ndimage.median_filter(np.absolute(X),size=(nMedianV+1,1),mode='constant')    
    
    Rt = np.divide(X_v_median,(X_v_median + X_h_median)+np.finfo(float).eps)
    #Rt[np.isnan(Rt)] = 0

    Rs = 1-Rt  

    S = np.sin(np.pi*(Rs-G2)/(2*(G1-G2)))**2    
    S[Rs>=G1] = 1
    S[Rs<G2] = 0

    T = np.sin(np.pi*(Rt-G2)/(2*(G1-G2)))**2
    T[Rt>=G1] = 1
    T[Rt<G2] = 0

    return S,T,1-S-T

# ------------------------------------------------
# Sum and Difference Related Functions (Not used)
# ------------------------------------------------

def normalize(x):
    return x/np.max(np.abs(x))

def sdTransform(x,N):
    xb_l = buffer(x[:,0],N)
    xb_r = buffer(x[:,1],N)
    frame_b = np.zeros(xb_l.shape[1])
    frame_c = np.zeros(xb_l.shape[1])
    S = np.zeros(xb_l.shape)
    D = np.zeros(xb_l.shape)
    for i in range(xb_l.shape[1]):        
        # Stereo balance
        frame_b[i] = stereoBalance(xb_l[:,i],xb_r[:,i])
        # Phase coherence
        frame_c[i] = phaseCoherence(xb_l[:,i],xb_r[:,i])

        S[:,i] = xb_l[:,i] + xb_r[:,i]
        if frame_b[i] > 0:
            D[:,i] = xb_l[:,i] - xb_r[:,i]
        else:
            D[:,i] = xb_r[:,i] - xb_l[:,i]

    y = np.transpose(np.array([np.matrix.flatten(np.transpose(S)),np.matrix.flatten(np.transpose(D))]))

    return y, frame_b, frame_c

def sdReconstruct(S,D,frame_b,frame_c,N):
    Sb = buffer(S,N)
    Db = buffer(D,N)
    L = np.zeros(Sb.shape)
    R = np.zeros(Sb.shape)
    frame_b_new = np.zeros(Sb.shape[1])
    frame_c_new = np.zeros(Sb.shape[1])
    for i in range(Sb.shape[1]):
        if frame_b[i] > 0:
            L[:,i] = 0.5*(Sb[:,i] + Db[:,i]);
            R[:,i] = 0.5*(Sb[:,i] - Db[:,i]);
        else:
            R[:,i] = 0.5*(Sb[:,i] + Db[:,i]);
            L[:,i] = 0.5*(Sb[:,i] - Db[:,i]);
        frame_b_new[i] = stereoBalance(L[:,i],R[:,i])
        frame_c_new[i] = phaseCoherence(L[:,i],R[:,i])
    y = np.transpose(np.array([np.matrix.flatten(np.transpose(L)),np.matrix.flatten(np.transpose(R))]))
    return y, frame_b_new, frame_c_new

def stereoBalance(L,R):
    Ln = np.divide(np.abs(L),np.max([np.abs(L),np.abs(R)]))
    Rn = np.divide(np.abs(R),np.max([np.abs(L),np.abs(R)]))
    res = Ln-Rn
    return np.mean(res)

def phaseCoherence(L,R):
    coherence = np.divide(np.multiply(L,R),np.max([np.abs(L),np.abs(R)]))
    coherence[coherence > 0] = 1
    coherence[coherence < 0] = -1
    return np.mean(coherence)

def buffer(X, n, p=0):
    d = n - p
    m = len(X)//d

    if m * d != len(X):
        m = m + 1

    Xn = np.zeros(d*m)
    Xn[:len(X)] = X

    Xn = np.reshape(Xn,(m,d))
    Xne = np.concatenate((Xn,np.zeros((1,d))))
    Xn = np.concatenate((Xn,Xne[1:,0:p]), axis = 1)

    return np.transpose(Xn[:-1])

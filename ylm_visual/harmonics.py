"""
github.com/itsmoosh/fun/ylm_visual/harmonics.py

Functions for real spherical harmonics.

Author: Marshall 'Moosh' Styczinski

Last updated: July 31, 2018
"""

import numpy as np

l0_norm = 1./np.sqrt(4.*np.pi)
l1_norm = np.sqrt(3./4./np.pi)
l2_norm = np.sqrt(15./16./np.pi)

def Y00(theta,phi):
	return l0_norm

def Y1_1(theta,phi):
	return l1_norm * np.sin(theta) * np.sin(phi)

def Y10(theta,phi):
	return l1_norm * np.cos(theta)

def Y11(theta,phi):
	return l1_norm * np.sin(theta) * np.cos(phi)

def Y2_2(theta,phi):
	return l2_norm * np.sin(theta)**2 * np.sin(2.*phi)

def Y2_1(theta,phi):
	return l2_norm * np.sin(2.*theta) * np.sin(phi)

def Y20(theta,phi):
	return l2_norm * np.sqrt(3.) * ( 2 - 3*np.sin(theta)**2 )

def Y21(theta,phi):
	return l2_norm * np.sin(2.*theta) * np.cos(phi)

def Y22(theta,phi):
	return l2_norm * np.sin(theta)**2 * np.cos(2.*phi)

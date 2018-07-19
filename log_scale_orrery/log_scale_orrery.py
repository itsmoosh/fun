"""
github.com/itsmoosh/fun/log_scale_orrery/log_scale_orrery.py

Required arguments:
	None

Please read the README file for more information.

Author: Marshall 'Moosh' Styczinski

Last updated (v1.0): July 18, 2018
"""

try:
	import ephem
except:
	print("Install the pyephem module to set date.")
	print("Showing the planets for the day Marshall married Nichole.")
	anyephem = False
else:
	anyephem = True

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

###########################
# Set the desired date here
# 	'YYYY/MM/DD'
###########################

date = '2000/01/01'
#date = '2000/01/01'
#date = '2012/07/21'	# Day Marshall married Nichole
#date = '2009/05/22'	# Day Marshall met Nichole
day_we_met = False	# Only used if pyephem is not installed

# Note: Elliptical orbits are drawn in the J2000 epoch,
#	so dates far from 2000 will have ellipses that don't
#	line up with the planets. Draw circular orbits by
#	setting circles = True below.

####################
# toggle log scaling
####################

scaling = True

#######################
# force circular orbits
#######################

circles = False

#############################
# number of asteroids to plot
#############################

num_ast = 300

#############################
# plot the Sun at the center?
#############################

show_sun = False


# Reference direction for ephemeris longitudes is the autumnal equinox
#	("out along ascending node of instantaneous plane
#	of the Earth's orbit and the Earth's mean equator
#	at the reference epoch")
# But I like having the x-axis be the vernal equinox.
# Adjust this value to rotate the entire diagram.

asc_node_long = np.pi


# plotting constants
deg = 180./np.pi	# Easily convert rads to degrees
log_scale = np.exp(1)	# Set the base for the logarithmic scaling here
draw_orbits = True	# Whether or not to draw circles/ellipses showing orbits


# The 'scale' parameter affects mostly just Mercury's orbit relative
#	to the others. If set too small, Mercury's orbit will not appear
#	as an ellipse. If scale = 1.0, Mercury, Venus, and Earth's orbits
#	will all be screwy because taking logs of numbers less than 1.0 (AU)
#	gives negative radii.
# The exact values chosen are based on aesthetics of Mercury's orbit.
if not circles:
	scale = log_scale*2.0
else:
	scale = log_scale*1.1


planet_color = 'black'
orbit_color = 'black'
sun_color = 'black'
ast_size = 0.5	# Asteroid size, in pt^2
planet_size = 40.0	# Planet size, in pt^2
sun_size = 100.0	# Sun size, in pt^2
orbit_thk = 1.5	# Line width for orbit paths, in pt
n_phis = 1000	# Number of angular points for orbit ellipses
fig_size = (6,6)	# In inches


# Source for orbital parameters: JPL Horizons ephemerides
#	https://ssd.jpl.nasa.gov/horizons.cgi

# orbital distances (semi-major axis, AU)
a_mer = 3.870987034266712E-01 * scale
a_ven = 7.233310387890227E-01 * scale
a_ear = 1.0 * scale
a_mar = 1.523671730968075E+00 * scale
a_jup = 5.201476244539330E+00 * scale
a_sat = 9.522205580101991E+00 * scale
a_ura = 1.925653891147859E+01 * scale
a_nep = 3.024327031338852E+01 * scale

# Eccentricities
e_mer = 2.056319033249083E-01
e_ven = 6.796560462447656E-03
e_ear = 1.628067641339369E-02
e_mar = 9.341829248052123E-02
e_jup = 4.877139014378658E-02
e_sat = 5.353004780964420E-02
e_ura = 4.374185069508277E-02
e_nep = 1.077861641590599E-02

# Longitude of periapsis calculated by argument of periapsis + longitude of ascending node
lpe_mer = ( 4.831944327526732E+01 + 2.915249832977428E+01 ) / deg + asc_node_long
lpe_ven = ( 7.665251797189694E+01 + 5.477679598608930E+01 ) / deg + asc_node_long
lpe_ear = ( 1.736244166246548E+02 + 2.916111755121159E+02 ) / deg + asc_node_long
lpe_mar = ( 4.952707806943124E+01 + 2.865227487815606E+02 ) / deg + asc_node_long
lpe_jup = ( 1.005195587183916E+02 + 2.742436131424251E+02 ) / deg + asc_node_long
lpe_sat = ( 1.135853142443672E+02 + 3.371436250079411E+02 ) / deg + asc_node_long
lpe_ura = ( 7.399148799178532E+01 + 9.818779753984455E+01 ) / deg + asc_node_long
lpe_nep = ( 1.317271854532969E+02 + 2.423914339165788E+02 ) / deg + asc_node_long

# Calculate latis rectum for each orbit, for convenient plotting later
p_mer = a_mer * (1 - e_mer**2)
p_ven = a_ven * (1 - e_ven**2)
p_ear = a_ear * (1 - e_ear**2)
p_mar = a_mar * (1 - e_mar**2)
p_jup = a_jup * (1 - e_jup**2)
p_sat = a_sat * (1 - e_sat**2)
p_ura = a_ura * (1 - e_ura**2)
p_nep = a_nep * (1 - e_nep**2)

# xy locations in AU, with x along vernal equinox direction
#	and y toward Earth northern summer
# longitudes from Earth vernal equinox location, in radians
if anyephem:
	# Fetch planet objects from PyEphem for the date we selected
	mer = ephem.Mercury(date)
	ven = ephem.Venus(date)
	ear = ephem.Sun(date)
	mar = ephem.Mars(date)
	jup = ephem.Jupiter(date)
	sat = ephem.Saturn(date)
	ura = ephem.Uranus(date)
	nep = ephem.Neptune(date)

	# Set solar distances and longitudes to variables we share with
	#	backup method
	r_mer = mer.sun_distance * scale
	r_ven = ven.sun_distance * scale
	r_ear = ear.earth_distance * scale
	r_mar = mar.sun_distance * scale
	r_jup = jup.sun_distance * scale
	r_sat = sat.sun_distance * scale
	r_ura = ura.sun_distance * scale
	r_nep = nep.sun_distance * scale

	l_mer = mer.hlon + asc_node_long
	l_ven = ven.hlon + asc_node_long
	l_ear = ear.hlon + asc_node_long
	l_mar = mar.hlon + asc_node_long
	l_jup = jup.hlon + asc_node_long
	l_sat = sat.hlon + asc_node_long
	l_ura = ura.hlon + asc_node_long
	l_nep = nep.hlon + asc_node_long

	# Find xy pos for plotting
	x_mer, y_mer = r_mer*np.cos(l_mer), r_mer*np.sin(l_mer)
	x_ven, y_ven = r_ven*np.cos(l_ven), r_ven*np.sin(l_ven)
	x_ear, y_ear = r_ear*np.cos(l_ear), r_ear*np.sin(l_ear)
	x_mar, y_mar = r_mar*np.cos(l_mar), r_mar*np.sin(l_mar)
	x_jup, y_jup = r_jup*np.cos(l_jup), r_jup*np.sin(l_jup)
	x_sat, y_sat = r_sat*np.cos(l_sat), r_sat*np.sin(l_sat)
	x_ura, y_ura = r_ura*np.cos(l_ura), r_ura*np.sin(l_ura)
	x_nep, y_nep = r_nep*np.cos(l_nep), r_nep*np.sin(l_nep)
	
else:
	if day_we_met:
		x_mer, y_mer = scale*-1.780678943214287E-01, scale*-4.290535224972859E-01
		x_ven, y_ven = scale* 6.734667596381133E-02, scale*-7.240360093685332E-01
		x_ear, y_ear = scale*-4.920363863124096E-01, scale*-8.846875849621418E-01
		x_mar, y_mar = scale* 1.383154748944253E+00, scale*-1.105786647570067E-01
		x_jup, y_jup = scale* 3.549324192629269E+00, scale*-3.615057601550965E+00
		x_sat, y_sat = scale*-9.281015323705194E+00, scale* 1.503535460396754E+00
		x_ura, y_ura = scale* 1.995674267159221E+01, scale*-2.369553099228205E+00
		x_nep, y_nep = scale* 2.441839661213937E+01, scale*-1.747935899530811E+01
	else:
		x_mer, y_mer = scale* 9.402137907119251E-02, scale*-4.438464690952168E-01
		x_ven, y_ven = scale* 6.076610278378199E-01, scale*-3.989783422976876E-01
		x_ear, y_ear = scale* 4.843988229200177E-01, scale*-8.932127376140822E-01
		x_mar, y_mar = scale*-1.051232539438123E+00, scale*-1.148032208668816E+00
		x_jup, y_jup = scale* 2.572466427514478E+00, scale* 4.305405530432761E+00
		x_sat, y_sat = scale*-8.511126719144398E+00, scale*-4.734832602431626E+00
		x_ura, y_ura = scale* 1.996887962522076E+01, scale* 1.962190311014376E+00
		x_nep, y_nep = scale* 2.632225010302883E+01, scale*-1.438045072575671E+01

	# Calculate current solar distance from xy ephemeris data
	r_mer = np.sqrt(x_mer**2 + y_mer**2)
	r_ven = np.sqrt(x_ven**2 + y_ven**2)
	r_ear = np.sqrt(x_ear**2 + y_ear**2)
	r_mar = np.sqrt(x_mar**2 + y_mar**2)
	r_jup = np.sqrt(x_jup**2 + y_jup**2)
	r_sat = np.sqrt(x_sat**2 + y_sat**2)
	r_ura = np.sqrt(x_ura**2 + y_ura**2)
	r_nep = np.sqrt(x_nep**2 + y_nep**2)

	# Calculate current solar longitude from Earth vernal equinox
	l_mer = np.arctan2(y_mer,x_mer) + asc_node_long
	l_ven = np.arctan2(y_ven,x_ven) + asc_node_long
	l_ear = np.arctan2(y_ear,x_ear) + asc_node_long
	l_mar = np.arctan2(y_mar,x_mar) + asc_node_long
	l_jup = np.arctan2(y_jup,x_jup) + asc_node_long
	l_sat = np.arctan2(y_sat,x_sat) + asc_node_long
	l_ura = np.arctan2(y_ura,x_ura) + asc_node_long
	l_nep = np.arctan2(y_nep,x_nep) + asc_node_long

# Array up the values for plotting
r_vals = [ r_mer, r_ven, r_ear, r_mar, r_jup, r_sat, r_ura, r_nep ]
l_vals = [ l_mer, l_ven, l_ear, l_mar, l_jup, l_sat, l_ura, l_nep ]
x_vals = [ x_mer, x_ven, x_ear, x_mar, x_jup, x_sat, x_ura, x_nep ]
y_vals = [ y_mer, y_ven, y_ear, y_mar, y_jup, y_sat, y_ura, y_nep ]

# Initialize arrays for asteroids
r_ast, l_ast, x_ast, y_ast = np.ones(num_ast), np.ones(num_ast), np.ones(num_ast), np.ones(num_ast)

# Make asteroids, which are centered around (in AU) 2.35 +- 0.1, 2.6 +- 0.05, 2.75 +- 0.03, 2.9 +- 0.03, 3.1 +- 0.1
# Source: eyeball estimates based on https://en.wikipedia.org/wiki/Asteroid_belt#/media/File:Kirkwood_Gaps.svg
for i in range(0, int(num_ast*0.4)):
	r_ast[i] = np.random.normal(loc=2.35, scale=0.1) * scale
for i in range(int(num_ast*0.4), int(num_ast*0.65)):
	r_ast[i] = np.random.normal(loc=2.6, scale=0.05) * scale
for i in range(int(num_ast*0.65), int(num_ast*0.75)):
	r_ast[i] = np.random.normal(loc=2.75, scale=0.03) * scale
for i in range(int(num_ast*0.75), int(num_ast*0.8)):
	r_ast[i] = np.random.normal(loc=2.9, scale=0.03) * scale
for i in range(int(num_ast*0.8), int(num_ast)):
	r_ast[i] = np.random.normal(loc=3.1, scale=0.1) * scale

# For each asteroid...
for i in range(num_ast):
	if scaling:
		# Adjust radial distances to log scale if we are doing that
		r_ast[i] = np.log(r_ast[i]) / np.log(log_scale)
	# Place at uniformly random longitudes
	l_ast[i] = np.random.rand()*2*np.pi
	# Find xy pos for plotting
	x_ast[i] = r_ast[i]*np.cos(l_ast[i])
	y_ast[i] = r_ast[i]*np.sin(l_ast[i])

# For each planet...
for i in range(len(r_vals)):
	if scaling:
		# Adjust radial distances to log scale if we are doing that
		r_vals[i] = np.log(r_vals[i])/np.log(log_scale)
	# Find xy pos for plotting
	x_vals[i] = r_vals[i]*np.cos(l_vals[i])
	y_vals[i] = r_vals[i]*np.sin(l_vals[i])

if circles:
	# Force planetary orbits to be circles, sized so that the planet's
	#	current location falls along the circle.
	mer_orbit = Circle((0.0,0.0), r_vals[0], edgecolor=orbit_color, fill=False)
	ven_orbit = Circle((0.0,0.0), r_vals[1], edgecolor=orbit_color, fill=False)
	ear_orbit = Circle((0.0,0.0), r_vals[2], edgecolor=orbit_color, fill=False)
	mar_orbit = Circle((0.0,0.0), r_vals[3], edgecolor=orbit_color, fill=False)
	jup_orbit = Circle((0.0,0.0), r_vals[4], edgecolor=orbit_color, fill=False)
	sat_orbit = Circle((0.0,0.0), r_vals[5], edgecolor=orbit_color, fill=False)
	ura_orbit = Circle((0.0,0.0), r_vals[6], edgecolor=orbit_color, fill=False)
	nep_orbit = Circle((0.0,0.0), r_vals[7], edgecolor=orbit_color, fill=False)


# Make the figure
plt.figure(figsize=fig_size)
fig = plt.gcf()
ax = fig.gca()
# Make it square
ax.set_aspect('equal')
# Plot planets
orrery = plt.scatter( x_vals , y_vals , color=planet_color, s=planet_size )
# Plot asteroids
asteroids = plt.scatter( x_ast , y_ast , color=planet_color, s=ast_size )

if show_sun:
	sun = plt.scatter( [0.0] , [0.0] , color=sun_color, s=sun_size )

if draw_orbits:
	if circles:
		# Add orbit circles to plot
		ax.add_artist(mer_orbit)
		ax.add_artist(ven_orbit)
		ax.add_artist(ear_orbit)
		ax.add_artist(mar_orbit)
		ax.add_artist(jup_orbit)
		ax.add_artist(sat_orbit)
		ax.add_artist(ura_orbit)
		ax.add_artist(nep_orbit)
	else:
		# Create and draw ellipses from the planets' orbital elements

		# Create a linear array of longitudes to draw planet ellipses
		phis = np.arange(0.0, 2.05*np.pi, 2*np.pi/n_phis)
		# Calculate planet distance for each longitude
		mer_rs = [ p_mer/(1+e_mer*np.cos(phi-lpe_mer)) for phi in phis ]
		ven_rs = [ p_ven/(1+e_ven*np.cos(phi-lpe_ven)) for phi in phis ]
		ear_rs = [ p_ear/(1+e_ear*np.cos(phi-lpe_ear)) for phi in phis ]
		mar_rs = [ p_mar/(1+e_mar*np.cos(phi-lpe_mar)) for phi in phis ]
		jup_rs = [ p_jup/(1+e_jup*np.cos(phi-lpe_jup)) for phi in phis ]
		sat_rs = [ p_sat/(1+e_sat*np.cos(phi-lpe_sat)) for phi in phis ]
		ura_rs = [ p_ura/(1+e_ura*np.cos(phi-lpe_ura)) for phi in phis ]
		nep_rs = [ p_nep/(1+e_nep*np.cos(phi-lpe_nep)) for phi in phis ]

		# Initialize arrays so we can reference/set specific elements
		mer_xs, mer_ys, ven_xs, ven_ys, ear_xs, ear_ys, mar_xs, mar_ys, jup_xs, jup_ys, sat_xs, sat_ys, ura_xs, ura_ys, nep_xs, nep_ys = np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis), np.ones(n_phis)

		for j in range(n_phis):
			if scaling:
				# Apply logarithmic scaling to the r(θ) values we derived
				mer_rs[j] = np.log( mer_rs[j] ) / np.log( log_scale )
				ven_rs[j] = np.log( ven_rs[j] ) / np.log( log_scale )
				ear_rs[j] = np.log( ear_rs[j] ) / np.log( log_scale )
				mar_rs[j] = np.log( mar_rs[j] ) / np.log( log_scale )
				jup_rs[j] = np.log( jup_rs[j] ) / np.log( log_scale )
				sat_rs[j] = np.log( sat_rs[j] ) / np.log( log_scale )
				ura_rs[j] = np.log( ura_rs[j] ) / np.log( log_scale )
				nep_rs[j] = np.log( nep_rs[j] ) / np.log( log_scale )

			# Find x pos for plotting
			mer_xs[j] = mer_rs[j] * np.cos(phis[j])
			ven_xs[j] = ven_rs[j] * np.cos(phis[j])
			ear_xs[j] = ear_rs[j] * np.cos(phis[j])
			mar_xs[j] = mar_rs[j] * np.cos(phis[j])
			jup_xs[j] = jup_rs[j] * np.cos(phis[j])
			sat_xs[j] = sat_rs[j] * np.cos(phis[j])
			ura_xs[j] = ura_rs[j] * np.cos(phis[j])
			nep_xs[j] = nep_rs[j] * np.cos(phis[j])
			# and y
			mer_ys[j] = mer_rs[j] * np.sin(phis[j])
			ven_ys[j] = ven_rs[j] * np.sin(phis[j])
			ear_ys[j] = ear_rs[j] * np.sin(phis[j])
			mar_ys[j] = mar_rs[j] * np.sin(phis[j])
			jup_ys[j] = jup_rs[j] * np.sin(phis[j])
			sat_ys[j] = sat_rs[j] * np.sin(phis[j])
			ura_ys[j] = ura_rs[j] * np.sin(phis[j])
			nep_ys[j] = nep_rs[j] * np.sin(phis[j])

		# Add calculated orbital ellipses to the plot
		mer_orbit = plt.plot( mer_xs, mer_ys, linewidth=orbit_thk, color=planet_color )
		ven_orbit = plt.plot( ven_xs, ven_ys, linewidth=orbit_thk, color=planet_color )
		ear_orbit = plt.plot( ear_xs, ear_ys, linewidth=orbit_thk, color=planet_color )
		mar_orbit = plt.plot( mar_xs, mar_ys, linewidth=orbit_thk, color=planet_color )
		jup_orbit = plt.plot( jup_xs, jup_ys, linewidth=orbit_thk, color=planet_color )
		sat_orbit = plt.plot( sat_xs, sat_ys, linewidth=orbit_thk, color=planet_color )
		ura_orbit = plt.plot( ura_xs, ura_ys, linewidth=orbit_thk, color=planet_color )
		nep_orbit = plt.plot( nep_xs, nep_ys, linewidth=orbit_thk, color=planet_color )

# Set axes to show everything
ax.set_xlim(-r_vals[7]*1.1, r_vals[7]*1.1)
ax.set_ylim(-r_vals[7]*1.1, r_vals[7]*1.1)
# Don't show scale markers or plot boundary
plt.axis('off')

# Let's see those planets!
plt.show()

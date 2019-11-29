# Load the necessary packages
import numpy as np
import os
from sympy import *
from sympy.geometry import *
import scipy.io as sio

######### OPTIONS ############
# Directory for the results
os.chdir('D:\\Dropbox')

# Dimension subvolume
k_sub = 40

# Length of rock surface line
k_line = 200

###############################################################################
############################ Functions needed #################################

### Function for Running mean/Moving average
def runningMean(x, N):
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
         y[ctr] = np.sum(x[ctr:(ctr+N)])
    return y/N

### Function to find the contact points from the segmented slice
def contactPoint(Segm):
	Segm_r = Segm.shape[0]
	Segm_c = Segm.shape[1]
	Oil = np.zeros([Segm_r,Segm_c])
	Brine = np.zeros([Segm_r,Segm_c])
	# Find points of rock confining with oil
	for i in range(1,Segm_c-1):
		for j in range(1,Segm_r-1):
			if Segm[j][i] == 3:
				if Segm[j][i-1] == 1 or Segm[j-1][i-1] == 1 or Segm[j+1][i-1] == 1 or Segm[j][i+1] == 1 or Segm[j+1][i+1] == 1 or Segm[j-1][i+1] == 1 or Segm[j-1][i] == 1 or Segm[j+1][i] == 1:
					Oil[j][i] = 1
	for i in range(1,Segm_c-1):
		for j in range(1,Segm_r-1):
			if Segm[j][i] == 3:
				if Segm[j][i-1] == 2 or Segm[j-1][i-1] == 2 or Segm[j+1][i-1] == 2 or Segm[j][i+1] == 2 or Segm[j+1][i+1] == 2 or Segm[j-1][i+1] == 2 or Segm[j-1][i] == 2 or Segm[j+1][i] == 2:
					Brine[j][i] = 1
	# Save in IPx and IPy the coordinates of rock/rest interface (points where
	# IP=2)
	IPx = []
	IPy = []
	IP = Oil+Brine; # Sum in order to find points confining with both o and w
	for i in range(np.shape(IP)[1]):
		for j in range(np.shape(IP)[0]):
			if IP[j][i] == 2:
				IPx.append(i)
				IPy.append(j)
	# Keep only one coordinate couple per each three phase point
	cancel = []
	it = 0
	for i in range(1,len(IPx)):
		if np.sqrt((IPx[i]-IPx[i-1])**2 + (IPy[i]-IPy[i-1])**2) == 1:
			cancel.append(i)
	if cancel != []:
		for i in range(len(cancel)):
			IPx = np.delete(IPx, cancel[i-it])
			IPy = np.delete(IPy, cancel[i-it])
			it = it + 1
	return IPx, IPy;

def contactAngle(Segm, IP, k_line):
	Segm_r = Segm.shape[0]
	Segm_c = Segm.shape[1]
	# Find the interface pixels between rock and the rest
	L = np.zeros([Segm_r, Segm_c])
	for i in range(1, Segm_r - 1):
		for j in range(1, Segm_c - 1):
			if Segm[j][i] == 3:
				if Segm[j][i-1] != 3 or Segm[j][i+1] != 3 or Segm[j+1][i] != 3 or Segm[j-1][i] != 3:
					L[j][i] = 1
	# Save the coordinates of the rock/rest interface
	Lx = []
	Ly = []
	for i in range(np.shape(L)[0]):
		for j in range(np.shape(L)[1]):
			if L[i][j] == 1:
				Lx.append(j)
				Ly.append(i)
	Lx = np.array(Lx)
	Ly = np.array(Ly)
	IPx = np.array(IP[0][0])
	IPy = np.array(IP[1][0])
	k_line = np.array(k_line)

	# Use the defined length for the rock/rest contact line
	Lx_new = Lx[(Lx > IPx - k_line) & (Lx < IPx + k_line) & (Ly > IPy - k_line) & (Ly < IPy + k_line)]
	Ly_new = Ly[(Lx > IPx - k_line) & (Lx < IPx + k_line) & (Ly > IPy - k_line) & (Ly < IPy + k_line)]

	# Calculate the interpolating contact line
	L_line = np.polyfit(Lx_new,Ly_new,1);
	m_line = L_line[0]
	q_line = L_line[1]
	Ly_fitted = m_line*Lx_new + q_line

	# Find the interface pixels between brine and rest
	C = np.zeros([Segm_r, Segm_c])
	for i in range(1, Segm_r - 1):
		for j in range(1, Segm_c - 1):
			if Segm[j][i] == 2:
				if Segm[j][i-1] == 1 or Segm[j][i+1] == 1 or Segm[j+1][i] == 1 or Segm[j-1][i] == 1:
					C[j][i] = 1
	# Save the coordinates of the interface
	Cx = []
	Cy = []
	for i in range(np.shape(C)[0]):
		for j in range(np.shape(C)[1]):
			if C[i][j] == 1:
				Cx.append(j)
				Cy.append(i)
	Cx = np.array(Cx)
	Cy = np.array(Cy)

	# Fit a circle to brine-oil contact points

	from scipy      import optimize

	x = Cx
	y = Cy

	# coordinates of the barycenter
	x_m = np.mean(x)
	y_m = np.mean(y)

	def calc_R(c):
		""" calculate the distance of each 2D points from the center c=(xc, yc) """
		return np.sqrt((x-c[0])**2 + (y-c[1])**2)

	def f_leastsq(c):
		""" calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
		Ri = calc_R(c)
		return Ri - np.mean(Ri)

	# Estimate the circle radius and center
	# Center
	center_estimate = x_m, y_m
	center, ier = optimize.leastsq(f_leastsq, center_estimate)
	xc, yc = center
	# Radius
	Ri       = calc_R(center)
	R        = np.mean(Ri)
	# Compute the residuals
	residu   = sum((Ri - R)**2)
	residu2  = sum((Ri**2-R**2)**2)

	line_start_p = Lx_new[0], Ly_fitted[0]
	line_end_p   = Lx_new[-1], Ly_fitted[-1]
	circle_centre_p = Point(xc, yc)
	line = Line(line_start_p, line_end_p)
	circle = Circle(circle_centre_p, R)
	line_circle_intersection = intersection(circle, line)
	for i in range(2):
		if line_circle_intersection[i].evalf(10)[0] > IPx-k_line and line_circle_intersection[i].evalf(10)[0] < IPx+k_line:
					   x_interception = float(line_circle_intersection[i].evalf(10)[0])
	y_interception = L_line[0] * x_interception + L_line[1]
	theta_secant = np.arctan(m_line)

	# Compute the contact angle

	cos_theta_circle = (x_interception - xc) / R
	sin_theta_circle = (y_interception - yc) / R
	theta_circle = np.arctan(sin_theta_circle / cos_theta_circle)

	theta_tangent = theta_circle + pi/2

	if theta_secant < 0:
		theta_secant = theta_secant + pi

	theta_contact_radiant = abs(theta_tangent - theta_secant)

	if theta_contact_radiant > pi/2:
		theta_contact_radiant = pi - theta_contact_radiant
		
	# Check the points inside the circle to see if the angle is acute or obtuse
	oil = 0
	water = 0
	area = 0
	for x in range(Segm_c):
		for y in range(Segm_r):
			if (abs(x - xc)**2 + abs(y - yc)**2)**2 < R**2:
				area = area + 1
				if Segm[y][x] == 1:
					oil = oil + 1
				if Segm[y][x] == 2:
					water = water + 1

	if water > oil:
		theta_contact_radiant = pi - theta_contact_radiant

	theta_contact_degree=float(deg(theta_contact_radiant))

	return (theta_contact_degree, R)

###############################################################################

seg_data = hx_project.get('Seged_data')

# Find three-phase points with label interface (getmultiphase) module
lab_interfaces = hx_project.create('getmultiphase')
lab_interfaces.ports.inputLabelImage.connect(seg_data)
lab_interfaces.ports.neighborhood._set_state('value 0')
lab_interfaces.ports.onlyBlackVoxels._set_state('values 1')
lab_interfaces.execute()
interfaces = lab_interfaces.results[0]


# Save the interfaces in a 3-dimensional matrix
IP_matrix = interfaces.get_array()

# Remove Avizo objects
hx_project.remove(lab_interfaces)
hx_project.remove(interfaces)

# Save interface points in variable IP
size_IP = np.shape(IP_matrix)
IP = [(x,y,z) for x in range(size_IP[0]) for y in range(size_IP[1]) for z in range(size_IP[2]) if IP_matrix[x][y][z] == 1]
IP_x = [IP[n][0] for n in range(np.shape(IP)[0])]
IP_y = [IP[n][1] for n in range(np.shape(IP)[0])]
IP_z = [IP[n][2] for n in range(np.shape(IP)[0])]

# Compute running mean on IP
# IP_x_movavg3 = runningMean(IP_x, 3)
# IP_y_movavg3 = runningMean(IP_y, 3)
# IP_z_movavg3 = runningMean(IP_z, 3)

# Compute normal directions
plane_nx = [IP_x[n + 1] - IP_x[n] for n in range(len(IP_x) - 1)]
plane_ny = [IP_y[n + 1] - IP_y[n] for n in range(len(IP_y) - 1)]
plane_nz = [IP_z[n + 1] - IP_z[n] for n in range(len(IP_z) - 1)]

k_ROI = 20;

# Compute the running mean inside small regions, to be more precise
IPx_R = np.zeros(len(IP_x))
IPy_R = np.zeros(len(IP_y))
IPz_R = np.zeros(len(IP_z))
nx_R = np.zeros(len(IP_x))
ny_R = np.zeros(len(IP_y))
nz_R = np.zeros(len(IP_z))

it = 0;
for i in range(len(IP_x)):
	# Consider only a small region
	# Define the borders (North South East West Upper Lower)
	x_E=IP_x[i]-k_ROI
	x_W=IP_x[i]+k_ROI
	y_S=IP_y[i]-k_ROI
	y_N=IP_y[i]+k_ROI
	z_L=IP_z[i]-k_ROI
	z_U=IP_z[i]+k_ROI
	if x_E < 0:
		x_E = 0
	if y_S < 0:
		y_S = 0
	if z_L < 0:
		z_L = 0
	if x_W > size_IP[0]:
		x_W = size_IP[0]
	if y_N > size_IP[1]:
		y_N = size_IP[1]
	if z_U > size_IP[2]:
		z_U = size_IP[2]
	
	IP_ROI = [(x,y,z) for x in range(x_E,x_W) for y in range(y_S,y_N) for z in range(z_L,z_U) if IP_matrix[x][y][z] == 1]
	IP_ROIx = [IP_ROI[n][0] for n in range(np.shape(IP_ROI)[0])]
	IP_ROIy = [IP_ROI[n][1] for n in range(np.shape(IP_ROI)[0])]
	IP_ROIz = [IP_ROI[n][2] for n in range(np.shape(IP_ROI)[0])]
	
	# if the contact point is not close to other contact points, it is not possible to compute normal directions
	try:
		# Compute running mean on IP_ROI
		IP_ROIx_movavg3 = runningMean(IP_ROIx, 5)
		IP_ROIy_movavg3 = runningMean(IP_ROIy, 5)
		IP_ROIz_movavg3 = runningMean(IP_ROIz, 5)
		
		# Compute normal directions
		plane_ROInx = [IP_ROIx_movavg3[n + 1] - IP_ROIx_movavg3[n] for n in range(len(IP_ROIx_movavg3) - 1)]
		plane_ROIny = [IP_ROIy_movavg3[n + 1] - IP_ROIy_movavg3[n] for n in range(len(IP_ROIy_movavg3) - 1)]
		plane_ROInz = [IP_ROIz_movavg3[n + 1] - IP_ROIz_movavg3[n] for n in range(len(IP_ROIz_movavg3) - 1)]
		
		IPx_R[it] = IP_ROIx[int(round(len(IP_ROIx_movavg3)/2))]
		IPy_R[it] = IP_ROIy[int(round(len(IP_ROIy_movavg3)/2))]
		IPz_R[it] = IP_ROIz[int(round(len(IP_ROIz_movavg3)/2))]
		nx_R[it] = plane_ROInx[int(round(len(IP_ROIx_movavg3)/2))]
		ny_R[it] = plane_ROIny[int(round(len(IP_ROIy_movavg3)/2))]
		nz_R[it] = plane_ROInz[int(round(len(IP_ROIz_movavg3)/2))]
		
		it = it + 1;
		
	except:
		it = it;
			
# Compute the contact angle for each position
total_points = len(IPx_R)
R = np.zeros(total_points)
theta = np.zeros(total_points)
conversion = np.zeros(total_points)
for i in range(total_points):
	# Extract subvolume
	subvolume = hx_project.create('HxLatticeAccess')
	subvolume.ports.data.connect(seg_data)
	subvolume.ports.boxMin._set_state('values 3 %.4f %.4f %.4f' % (IPx_R[i] - k_sub/2, IPy_R[i] - k_sub/2, IPz_R[i] - k_sub/2))
	subvolume.ports.boxSize._set_state('values 3 %.4f %.4f %.4f' % (k_sub, k_sub, k_sub))
	subvolume.execute()
	subvolume.viewer_mask = 0
	seg_data_view = subvolume.results[0]

	# Get the slice
	slice = hx_project.create('HxFilteredObliqueSlice')
	slice.viewer_mask = 0
	slice.ports.data.connect(seg_data_view)
	slice.fire()
	slice.ports.planeNormal._set_state('values 3 %.4f %.4f %.4f' % (nx_R[i], ny_R[i], nz_R[i]))
	slice.ports.planePoint1._set_state('points 1  %.4f %.4f %.4f' % (IPx_R[i], IPy_R[i], IPz_R[i]))
	slice.fire()

	# Extract image from the slice
	image = hx_project.create('HxCreateImage')
	image.ports.data.connect(slice)
	image.fire()
	imaged = image.results[0]
	Segm = imaged.get_array()

	# Find the three phase contact points from segmented images
	IP = np.zeros(2)
	IP = contactPoint(Segm)
	# Compute the contact angle
	try:
		(theta[i], R[i]) = contactAngle(Segm, IP, k_line)
	except:
		(theta[i], R[i]) = (np.NAN, np.NAN)

		
	# Remove Avizo objects
	hx_project.remove(subvolume)
	hx_project.remove(seg_data_view)
	hx_project.remove(slice)
	hx_project.remove(image)
	hx_project.remove(imaged)

sio.savemat('Contact_angles_ROI.mat',{'theta_p':theta})
sio.savemat('IPx.mat',{'IPx': IPx_R})
sio.savemat('IPy.mat',{'IPy': IPy_R})
sio.savemat('IPz.mat',{'IPz': IPz_R})
sio.savemat('nx.mat',{'nx': nx_R})
sio.savemat('ny.mat',{'ny': ny_R})
sio.savemat('nz.mat',{'nz': nz_R})


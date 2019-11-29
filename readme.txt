0) Install the sympy package in Avizo, using Enstaller. The instructions to do so can be found in the Avizo User's guide at the page "Python Package Manager: Enstaller".

1) Obtain segmented data, with the following specifications:

	- Rock phase labelled as phase 3
	- Denser phase labbelled as phase 2
	- Lighter phase labelled as phase 1
	
	With the crop editor of Avizo please set:
	- Minimum coordinates: 	0	0	0
	- Voxel size:			1	1	1
	
	- In Avizo, rename segmented data as 'Seged_data'
	
2) Copy and paste the code in the Python console of Avizo, with these optional changes:

	- Line 10: Provide the path of your working folder (where data will be saved)
	- Line 13: Provide the dimension of the subvolumes to be extracted around each three-phase contact point (called 'b' in Scanziani et al., 2017)
	- Line 16: Provide the length of the straight line representing the rock surface (called 'l' in Scanziani et al., 2017)
	
3) Wait until Avizo finishes its calculations. While the code is running, Avizo is freezed and does not respond. Computational time for a dataset with 2 billion voxels is about 1 day (24h).

4) Once the calculations are completed, you will find in your working folder 7 matlab databases:

	- Contact_angles_ROI.mat, the vector of contact angles.
	- IPx.mat, the x-coordinates of each three-phase contact point where the contact angles are computed.
	- IPy.mat, the y-coordinates of each three-phase contact point where the contact angles are computed.
	- IPz.mat, the z-coordinates of each three-phase contact point where the contact angles are computed.
	- nx.mat, the x-components of the unit vectors identifying the normal directions of the planes where the contact angles are computed.
	- ny.mat, the y-components of the unit vectors identifying the normal directions of the planes where the contact angles are computed.
	- nz.mat, the z-components of the unit vectors identifying the normal directions of the planes where the contact angles are computed.
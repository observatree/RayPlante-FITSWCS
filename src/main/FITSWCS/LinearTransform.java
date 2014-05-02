/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS;

import FITSWCS.exceptions.SingularMatrixException;

/**
 *   This class provides support for linear coordinate
 *   transformations used by the FITS "World Coordinate System" (WCS) 
 *   convention. <p>
 *
 *   The FITSWCS package was translated from the WCSLIB C library
 *   This original library was written in support for coordinate 
 *   systems used by astronomical data stored in FITS format.  For more 
 *   information on these coordinate systems, refer to the paper by Greisen 
 *   and Calabretta at:
 *   <blockquote>
 *       ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z 
 *   </blockquote>
 *
 *   Note that the functionality the linprm data structrue in the 
 *   WCSLIB C library has been incorporated into the data portion of this
 *   class; similarly the functionality of linset() has been put into
 *   the constructor.
 *
 *   <b> Nomenclature </b><p>
 *
 *   In WCSLIB the "forward" direction is from (lng,lat) celestial
 *   coordinates to (x,y) coordinates in the plane of projection.  This
 *   accords with the notion that spherical projections are a projection of
 *   the sphere onto a plane, the "reverse" direction is therefore that of
 *   deprojection from plane to sphere. <p>
 *   
 *   Unfortunately, this is opposite to what is generally understood to be
 *   the forward direction for FITS, namely that of transforming pixel
 *   coordinates to world coordinates.  However, the ordering of function
 *   argument lists should make it clear what is intended.
 *   
 *   <hr>
 *
 *   <b> COPYRIGHT NOTICE </b><p>
 *
 *   This library is free software; you can redistribute it and/or modify it
 *   under the terms of the GNU Library General Public License as published
 *   by the Free Software Foundation; either version 2 of the License, or (at
 *   your option) any later version. <p>
 *
 *   This library is distributed in the hope that it will be useful, but
 *   WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library
 *   General Public License for more details. <p>
 *
 *   You should have received a copy of the GNU Library General Public License
 *   along with this library; if not, write to the Free Software Foundation,
 *   Inc., 675 Mass Ave, Cambridge, MA 02139, USA. <p>
 *
 *   Correspondence concerning WCSLIB may be directed to:<pre>
 *      Internet email: mcalabre@atnf.csiro.au
 *      Postal address: Dr. Mark Calabretta,
 *                      Australia Telescope National Facility,
 *                      P.O. Box 76,
 *                      Epping, NSW, 2121,
 *                      AUSTRALIA</pre>
 *   Correspondence concerning the Java implementation may be directed
 *   to Raymond L. Plante (rplante@ncsa.uiuc.edu).
 *
 * @author Mark Calabretta
 * @version 2.4
 *---------------------------------------------------------------------------*/
public class LinearTransform {

    protected int naxis=0;
    protected double[] crpix, pc, cdelt;

    // Intermediates
    protected double[] piximg;
    protected double[] imgpix;

    /**
     * create a LinearTransform object given the matrices describing the 
     * coordinate system.
     * @param naxis Number of image axes
     * @param crpix naxis-length array containing the coordinate 
     *              reference pixel (CRPIXn).
     * @param pc    (naxis*naxis)-length array containing the elements 
     *              of the PC (pixel coordinate) transformation matrix;
     *              the elements should be arranged such that the first 
     *              axis is the most rapidly varying.  
     * @param cdelt naxis-length array containing the coordinate 
     *              increments (CDELTn).<p>
     * @exception ArrayIndexOutOfBoundsException if naxis < 1, 
     *              crpix.length < naxis, cdelt.length < naxis, or
     *              pc is not an naxis by naxis 2D array
     * @exception SingularMatrixException if pc represents a singular
     *              matrix.
     */
    public LinearTransform(int naxis, double[] crpix, double[][] pc, 
			   double[] cdelt) 
	throws ArrayIndexOutOfBoundsException, SingularMatrixException
    {
	int i, j;
	if (naxis < 1) throw new ArrayIndexOutOfBoundsException("naxis: " +
	                                                        naxis);
	double[] use = new double[naxis*naxis];
	if (pc == null) {
	    pc = new double[1][1];
	    pc[0][0] = 1.0;
	}
	
	for(i=0; i < pc.length; i++) {
	    for(j=0; j < pc[i].length; j++) {
		use[i*naxis+j] = pc[i][j];
	    }
	    for(j=pc[i].length; j < naxis; j++) {
		use[i*naxis+j] = (i == j) ? 1.0 : 0.0;
	    }
	}
	for(i=pc.length; i < naxis; i++) {
	    for(j=0; j < naxis; j++) {
		use[i*naxis+j] = (i == j) ? 1.0 : 0.0;
	    }
	}

	init(naxis, crpix, use, cdelt);
    }
	
    /**
     * create a LinearTransform object given the matrices describing the 
     * coordinate system.
     * @param naxis Number of image axes
     * @param crpix naxis-length array containing the coordinate 
     *              reference pixel (CRPIXn).
     * @param pc    (naxis*naxis)-length array containing the elements 
     *              of the PC (pixel coordinate) transformation matrix;
     *              the elements should be arranged such that the first 
     *              axis is the most rapidly varying.  
     * @param cdelt naxis-length array containing the coordinate 
     *              increments (CDELTn).<p>
     * @exception ArrayIndexOutOfBoundsException if naxis < 1, 
     *              crpix.length < naxis, cdelt.length < naxis, or
     *              pc.lenght < naxis*naxis.
     * @exception SingularMatrixException if pc represents a singular
     *              matrix.
     */
    public LinearTransform(int naxis, double[] crpix, double[] pc, 
			   double[] cdelt) 
	throws ArrayIndexOutOfBoundsException, SingularMatrixException
    {
	init(naxis, crpix, pc, cdelt);
    }

    /**
     * create a LinearTransform object given the matrices describing the 
     * coordinate system.
     * @param naxis Number of image axes
     * @param crpix naxis-length array containing the coordinate 
     *              reference pixel (CRPIXn).
     * @param cdelt naxis-length array containing the coordinate 
     *              increments (CDELTn).<p>
     * @exception ArrayIndexOutOfBoundsException if naxis < 1, 
     *              crpix.length < naxis, cdelt.length < naxis, or
     *              pc.lenght < naxis*naxis.
     * @exception SingularMatrixException if pc represents a singular
     *              matrix.
     */
    public LinearTransform(int naxis, double[] crpix, double[] cdelt) 
	throws ArrayIndexOutOfBoundsException, SingularMatrixException
    {
	this(naxis, crpix, (double[][])null, cdelt);
    }

    private void init(int naxis, double[] crpix, double[] pc, 
			   double[] cdelt) 
	throws ArrayIndexOutOfBoundsException, SingularMatrixException
    {
	int i, ij, j;

	// set naxis
	if (naxis < 1) throw new ArrayIndexOutOfBoundsException("naxis: " +
	                                                        naxis);
	this.naxis = naxis;

	// copy crpix: set missing elements to a "sensible" default (1)
	this.crpix = new double[naxis];
	System.arraycopy(crpix, 0, this.crpix, 0, 
			 Math.min(crpix.length, this.crpix.length));
	for(i=crpix.length; i < this.crpix.length; i++) {
	    this.crpix[i] = 1;
	}
	
	// copy cdelt: set missing elements to a "sensible" default (1)
	this.cdelt = new double[naxis];
	System.arraycopy(cdelt, 0, this.cdelt, 0, 
			 Math.min(cdelt.length, naxis));
	for(i=cdelt.length; i < naxis; i++) {
	    this.cdelt[i] = 1;
	}

	// copy pc matrix: set missing elements to a "sensible" default (1)
	j = naxis*naxis;
	this.pc = new double[j];
	System.arraycopy(pc, 0, this.pc, 0, Math.min(pc.length, j));
	for(i=pc.length; i < j; i++) {
	    this.pc[i] = 1;
	}

	// set piximg, imgpix matrices

	// Allocate memory for internal arrays. 
	piximg = new double[naxis*naxis];

	// Compute the pixel-to-image transformation matrix. 
	for (i = 0, ij = 0; i < naxis; i++) {
	    for (j = 0; j < naxis; j++, ij++) {
		piximg[ij] = cdelt[i] * pc[ij];
	    }
	}

	// Compute the image-to-pixel transformation matrix. 
	try {
	    imgpix = matinv(naxis, piximg);
	} catch (SingularMatrixException ex) {
	    throw new SingularMatrixException("PC matrix is singular");
	}
    }

    /**
     * An auxiliary matrix inversion routine, matinv().  It uses
     *   LU-triangular factorization with scaled partial pivoting.
     * @param n    number of dimensions of the input matrix
     * @param mat  (naxis*naxis)-length array containing the input matrix
     * @return double[] (naxis*naxis)-length array containing the 
     *             inversion of mat
     * @exception ArrayIndexOutOfBoundsException if n < 1 or
     *              mat.length < naxis*naxis.
     * @exception SingularMatrixException if mat represents a singular
     *              matrix.
     */
    public final static double[] matinv(int n, double[] mat) 
	throws ArrayIndexOutOfBoundsException, SingularMatrixException
    {
	int i, ij, ik, j, k, kj, pj;
	int  itemp, pivot;
	int[] mxl, lxm;
	double colmax, dtemp;
	double[] lu, rowmax;
	double[] inv;

	if (n < 1) throw new ArrayIndexOutOfBoundsException("n: " + n);
	if (mat.length < n*n) 
	    throw new ArrayIndexOutOfBoundsException("mat: " + 
		                                     mat.length + " < n*n");

	// Allocate memory for internal arrays. 
	mxl = new int[n];
	lxm = new int[n];
	rowmax = new double[n];
	lu = new double[n*n];
	inv = new double[n*n];

	// Initialize arrays. 
	for (i = 0, ij = 0; i < n; i++) {

	    // Vector which records row interchanges. 
	    mxl[i] = i;

	    rowmax[i] = 0.0;

	    for (j = 0; j < n; j++, ij++) {
		dtemp = Math.abs(mat[ij]);
		if (dtemp > rowmax[i]) rowmax[i] = dtemp;

		lu[ij] = mat[ij];
	    }

	    // A row of zeroes indicates a singular matrix. 
	    if (rowmax[i] == 0.0) throw new SingularMatrixException();

	}

	// Form the LU triangular factorization using scaled partial pivoting. 
	for (k = 0; k < n; k++) {
	    
	    // Decide whether to pivot. 
	    colmax = Math.abs(lu[k*n+k]) / rowmax[k];
	    pivot = k;

	    for (i = k+1; i < n; i++) {
		ik = i*n + k;
		dtemp = Math.abs(lu[ik]) / rowmax[i];
		if (dtemp > colmax) {
		    colmax = dtemp;
		    pivot = i;
		}
	    }

	    if (pivot > k) {
		
		// We must pivot, interchange the rows of the design matrix. 
		for (j = 0, pj = pivot*n, kj = k*n; j < n; j++, pj++, kj++) {
		    dtemp = lu[pj];
		    lu[pj] = lu[kj];
		    lu[kj] = dtemp;
		}

		// Amend the vector of row maxima. 
		dtemp = rowmax[pivot];
		rowmax[pivot] = rowmax[k];
		rowmax[k] = dtemp;

		// Record the interchange for later use. 
		itemp = mxl[pivot];
		mxl[pivot] = mxl[k];
		mxl[k] = itemp;
	    }

	    // Gaussian elimination. 
	    for (i = k+1; i < n; i++) {
		ik = i*n + k;

		// Nothing to do if lu[ik] is zero. 
		if (lu[ik] != 0.0) {

		    // Save the scaling factor. 
		    lu[ik] /= lu[k*n+k];

		    // Subtract rows. 
		    for (j = k+1; j < n; j++) {
			lu[i*n+j] -= lu[ik]*lu[k*n+j];
		    }
		}
	    }
	}

	// mxl[i] records which row of mat corresponds to row i of lu.  
	// lxm[i] records which row of lu  corresponds to row i of mat. 
	for (i = 0; i < n; i++) {
	    lxm[mxl[i]] = i;
	}

	// Determine the inverse matrix. 
	for (i = 0, ij = 0; i < n; i++) {
	    for (j = 0; j < n; j++, ij++) {
		inv[ij] = 0.0;
	    }
	}

	for (k = 0; k < n; k++) {
	    inv[lxm[k]*n+k] = 1.0;

	    // Forward substitution. 
	    for (i = lxm[k]+1; i < n; i++) {
		for (j = lxm[k]; j < i; j++) {
		    inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
		}
	    }

	    // Backward substitution. 
	    for (i = n-1; i >= 0; i--) {
		for (j = i+1; j < n; j++) {
		    inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
		}
		inv[i*n+k] /= lu[i*n+i];
	    }
	}

	return inv;
    }

    /**
     * Compute pixel coordinates from image coordinates.  Note that where
     * celestial coordinate systems are concerned the image coordinates
     * correspond to (x,y) in the plane of projection, not celestial (lng,lat).
     * @param imgcoord array representing the input position; missing elements 
     *                 will default to 1.0;
     * @return double[] array representing the transformed position
     */
    public double[] fwd(double[] imgcoord) {

	double[] pixcrd = new double[naxis];
	double[] imgcrd;
	if (imgcoord.length < naxis) {
	    imgcrd = new double[naxis];
	    System.arraycopy(imgcoord, 0, imgcrd, 0, 
			     Math.min(imgcoord.length,naxis));
	    for(int k=imgcoord.length; k < naxis; k++) 
		imgcrd[k] = 1.0;
	}
	else {
	    imgcrd = imgcoord;
	}

	int i, ij, j;

	for (i = 0, ij = 0; i < naxis; i++) {
	    pixcrd[i] = 0.0;
	    for (j = 0; j < naxis; j++, ij++) {
		pixcrd[i] += imgpix[ij] * imgcrd[j];
	    }
	}

	for (j = 0; j < naxis; j++) {
	    pixcrd[j] += crpix[j];
	}

	return pixcrd;
    }

    /**
     * Compute image coordinates from pixel coordinates.  Note that where
     * celestial coordinate systems are concerned the image coordinates
     * correspond to (x,y) in the plane of projection, not celestial (lng,lat).
     * @param pixcoord array representing the input position; missing elements 
     *                 will default to 1.0;
     * @return double[] array representing the transformed position
     */
    public double[] rev(double[] pixcoord) {

	double[] imgcrd = new double[naxis];
	double[] pixcrd;
	if (pixcoord.length < naxis) {
	    pixcrd = new double[naxis];
	    System.arraycopy(pixcoord, 0, imgcrd, 0, 
			     Math.min(pixcoord.length,naxis));
	    for(int k=pixcoord.length; k < naxis; k++) 
		imgcrd[k] = 1.0;
	}
	else {
	    pixcrd = pixcoord;
	}

	int i, ij, j;
	double temp;

	for (i = 0; i < naxis; i++) {
	    imgcrd[i] = 0.0;
	}

	for (j = 0; j < naxis; j++) {
	    temp = pixcrd[j] - crpix[j];
	    for (i = 0, ij = j; i < naxis; i++, ij+=naxis) {
		imgcrd[i] += piximg[ij] * temp;
	    }
	}

	return imgcrd;
    }
	
}

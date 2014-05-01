/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS.projections;

import FITSWCS.*;
import FITSWCS.exceptions.*;

/**
 *   This class provides support for the zenithal/azimuthal polynomial
 *   projection (ZPN) used by the FITS "World Coordinate System" (WCS) 
 *   convention. <p>
 *
 *   The FITSWCS package was translated from the WCSLIB C library.
 *   This original library was written in support for coordinate 
 *   systems used by astronomical data stored in FITS format.  For more 
 *   information on these coordinate systems, refer to the paper by Greisen 
 *   and Calabretta at:
 *   <blockquote>
 *       ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z 
 *   </blockquote>
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
 *   argument lists should make it clear what is intended.<p>
 *   
 *   <b> Accuracy </b><p>
 *   
 *   Closure to a precision of at least 1.0-10 degree of longitude and latitude
 *   has been verified for typical projection parameters on the 1 degree grid
 *   of native longitude and latitude (to within 5 degrees of any latitude
 *   where the projection may diverge).
 *
 *   Notwithstanding this, absolutely no claim is made for the accuracy or
 *   reliability of these routines.  They are supplied as is, with no warranty
 *   of fitness for any purpose.
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
public class ZPNProjection extends Projection {

    protected final static double tol = 1.0e-13;

    /**
     * Create an ZPNProjection object
     * @param r0 sphere radius
     * @param p  up to 10 (and at least 1) projection parameters
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public ZPNProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(r0, p);
    }

    /**
     * Create an ZPNProjection object
     * @param p  up to 10 (and at least 1) projection parameters
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public ZPNProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public ZPNProjection() { r0 = R2D; }

    /**
     * sets up intermediate data
     * <pre>
     *   r0    reset to 180/pi if 0.
     *   n     degree of polynomial
     *   w[0]  Co-latitude of the first point of inflection (N > 2).
     *   w[1]  Radius of the first point of inflection (N > 2).
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	int   i, j, k;
	double d, d1, d2, r, zd, zd1, zd2;

	if (p == null || p.length < 10) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "10 projection parameters");

	if (r0 == 0.0) r0 = R2D;
	this.r0 = r0;
	this.w = new double[2];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	// Find the highest non-zero coefficient. 
	for (k = 9; k >= 0 && p[k] == 0.0; k--);
	if (k < 0) throw new 
	    BadProjectionParameterException("All coefficients = 0");
	n = k;

	if (k >= 3) {

	    // Find the point of inflection closest to the pole. 
	    zd1 = zd2 = d2 = zd = 0.0;
	    d1  = p[1];
	    if (d1 <= 0.0) throw new 
		BadProjectionParameterException("p[1] <= 0");

	    // Find the point where the derivative first goes negative. 
	    for (i = 0; i < 180; i++) {
		zd2 = i*D2R;
		d2  = 0.0;
		for (j = k; j > 0; j--) {
		    d2 = d2*zd2 + j*p[j];
		}

		if (d2 <= 0.0) break;
		zd1 = zd2;
		d1  = d2;
	    }

	    if (i == 180) {

		// No negative derivative -> no point of inflection. 
		zd = PI;
	    } else {
		
		// Find where the derivative is zero. 
		for (i = 1; i <= 10; i++) {
		    zd = zd1 - d1*(zd2-zd1)/(d2-d1);

		    d = 0.0;
		    for (j = k; j > 0; j--) {
			d = d*zd + j*p[j];
		    }

		    if (Math.abs(d) < tol) break;

		    if (d < 0.0) {
			zd2 = zd;
			d2  = d;
		    } else {
			zd1 = zd;
			d1  = d;
		    }
		}
	    }

	    r = 0.0;
	    for (j = k; j >= 0; j--) {
		r = r*zd + p[j];
	    }
	    w[0] = zd;
	    w[1] = r;
	}
    }	

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
    {
	int   j;
	double r, s;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	s = (90.0 - theta)*D2R;
	r = 0.0;
	for (j = 9; j >= 0; j--) {
	    r = r*s + p[j];
	}
	r = r0*r;

	out[0] =  r*TrigD.sin(phi);
	out[1] = -r*TrigD.cos(phi);
	return out;
    }

    /**
     * same as fwd(phitheta[0], phitheta[1])
     */
    public double[] fwd(double[] phitheta) { 
	return fwd(phitheta[0], phitheta[1]); 
    }

    /**
     * Compute native spherical coordinates (phi,theta) from the 
     * (x,y) coordinates in the plane of projection. 
     * @return double[] a two-element array containing phi,theta
     */
    public double[] rev(double x, double y) 
	throws PixelBeyondProjectionException
    {
	int i, j, k;
	double a, b, c, d, lambda, r, r1, r2, rt, zd=0, zd1, zd2;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	r = Math.sqrt(x*x + y*y)/r0;

	if (n == 1) {

	    // Linear. 
	    zd = (r - p[0])/p[1];

	} else if (n == 2) {

	    // Quadratic. 
	    a = p[2];
	    b = p[1];
	    c = p[0] - r;

	    d = b*b - 4.0*a*c;
	    if (d < 0.0) throw new
		PixelBeyondProjectionException("ZPN: (x,y) = (" + x +
					       ", " + y + ")");
	    d = Math.sqrt(d);

	    // Choose solution closest to pole. 
	    zd1 = (-b + d)/(2.0*a);
	    zd2 = (-b - d)/(2.0*a);
	    zd  = (zd1<zd2) ? zd1 : zd2;
	    if (zd < -tol) zd = (zd1>zd2) ? zd1 : zd2;
	    if (zd < 0.0) {
		if (zd < -tol) throw new
		    PixelBeyondProjectionException("ZPN: (x,y) = (" + x +
						   ", " + y + ")");
		zd = 0.0;
	    } else if (zd > PI) {
		if (zd > PI+tol) throw new
		    PixelBeyondProjectionException("ZPN: (x,y) = (" + x +
						   ", " + y + ")");
		zd = PI;
	    }
	} else {

	    // Higher order - solve iteratively. 
	    zd1 = 0.0;
	    r1  = p[0];
	    zd2 = w[0];
	    r2  = w[1];

	    if (r < r1) {
		if (r < r1-tol) throw new
		    PixelBeyondProjectionException("ZPN: (x,y) = (" + x +
						   ", " + y + ")");
		zd = zd1;
	    } else if (r > r2) {
		if (r > r2+tol) throw new
		    PixelBeyondProjectionException("ZPN: (x,y) = (" + x +
						   ", " + y + ")");
		zd = zd2;
	    } else {

		// Disect the interval. 
		for (j = 0; j < 100; j++) {
		    lambda = (r2 - r)/(r2 - r1);
		    if (lambda < 0.1) {
			lambda = 0.1;
		    } else if (lambda > 0.9) {
			lambda = 0.9;
		    }

		    zd = zd2 - lambda*(zd2 - zd1);

		    rt = 0.0;
		    for (i = n; i >= 0; i--) {
			rt = (rt * zd) + p[i];
		    }

		    if (rt < r) {
			if (r-rt < tol) break;
			r1 = rt;
			zd1 = zd;
		    } else {
			if (rt-r < tol) break;
			r2 = rt;
			zd2 = zd;
		    }

		    if (Math.abs(zd2-zd1) < tol) break;
		}
	    }
	}

	if (r == 0.0) {
	    out[0] = 0.0;
	} else {
	    out[0] = TrigD.atan2(x, -y);
	}
	out[1] = 90.0 - zd*R2D;

	return out;
    }

    /**
     * set the sphere radius 
     */
    public void setR0(double r0) {
	this.r0 = (r0 == 0) ? R2D : r0;
	try {
	    if (p != null) init(r0, p);
	} catch (ArrayIndexOutOfBoundsException ex) {
	    p = null;
	} catch (BadProjectionParameterException ex) {
	    p = null;
	}
    }

    /**
     * set the projection parameters 
     * @param p array of projection parameters; a null value cause no updates
     *          to be made.
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public void setProjParm(double[] p)
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	if (p == null) return;
	if (p.length < 10) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "10 projection parameters");
	this.p = null;
	init(r0, p);
    }
}



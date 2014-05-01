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
 *   This class provides support for the Airy's zenithal 
 *   projection (AIR) used by the FITS "World Coordinate System" (WCS) 
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
public class AIRProjection extends Projection {

    protected final static double tol = 1.0e-4;

    /**
     * Create an AIRProjection object
     * @param r0 sphere radius; if 0, it defaults to 180/PI 
     * @param p  Latitude theta_b within which the error is minimized,
     *           in degrees.
     * @exception BadProjectionParameterException if p is outside of
     *            the range (-90,90].
     */
    public AIRProjection(double r0, double p) 
	throws BadProjectionParameterException
    {
	double[] use = { 0, p };
	init(r0, use);
    }

    /**
     * Create an AIRProjection object
     * @param p  Latitude theta_b within which the error is minimized,
     *           in degrees.
     * @exception BadProjectionParameterException if p is outside of
     *            the range (-90,90].
     */
    public AIRProjection(double p) 
	throws BadProjectionParameterException
    {
	this(0, p);
    }

    /**
     * Create an AIRProjection object
     * @param r0 sphere radius; if 0, it defaults to 180/PI 
     * @param p  projection parameters, where p[1] is the 
     *           AIR distance parameter, mu in units of r0
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] is outside of
     *            the range (-90,90].
     */
    public AIRProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(r0, p);
    }

    /**
     * Create an AIRProjection object
     * @param p  projection parameters, where p[1] is the 
     *           Latitude theta_b within which the error is minimized,
     *           in degrees
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] is outside of
     *            the range (-90,90].
     */
    public AIRProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	this(0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public AIRProjection() { r0 = R2D; }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   ln(cos(xi_b))/tan(xi_b)**2, where xi_b = (90-theta_b)/2 
     *   w[1]   1/2 - w[0]                                              
     *   w[2]   r0*w[1]                                                 
     *   w[3]   tol, cutoff for using small angle approximation, in     
     *	             radians.                                           
     *	 w[4]   w[1]*tol                                                
     *	 w[5]   (180/pi)/w[1]                                           
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	double cxi, r, txi, xi;
	if (p == null || p.length < 2) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "2 projection parameters");

	if (r0 == 0.0) r0 = R2D;
	this.r0 = r0;
	this.w = new double[6];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	if (p[1] == 90.0) {
	    w[0] = -0.5;
	    w[1] =  1.0;
	} else if (p[1] > -90.0) {
	    cxi = TrigD.cos((90.0 - p[1])/2.0);
	    w[0] = Math.log(cxi)*(cxi*cxi)/(1.0-cxi*cxi);
	    w[1] = 0.5 - w[0];
	} else {
	    throw new BadProjectionParameterException("p[1] out of range " +
						      "(-90,90]: " + p[1]);
	}

	w[2] = r0 * w[1];
	w[3] = tol;
	w[4] = w[1]*tol;
	w[5] = R2D/w[1];
    }	

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
	throws PixelBeyondProjectionException
    {
	double cxi, r, txi, xi;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	if (theta == 90.0) {
	    r = 0.0;
	} else if (theta > -90.0) {
	    xi = D2R*(90.0 - theta)/2.0;
	    if (xi < w[3]) {
		r = xi*w[2];
	    } else {
		cxi = TrigD.cos((90.0 - theta)/2.0);
		txi = Math.sqrt(1.0-cxi*cxi)/cxi;
		r = -r0*(Math.log(cxi)/txi + w[0]*txi);
	    }
	} else {
	    throw new PixelBeyondProjectionException("AIR: angle out of " +
						     "bounds: theta = " +
						     theta);
	}

	out[0] =  r*TrigD.sin(phi);  // x
	out[1] = -r*TrigD.cos(phi);  // y
	return out;
    }

    /**
     * Compute native spherical coordinates (phi,theta) from the 
     * (x,y) coordinates in the plane of projection. 
     * @return double[] a two-element array containing phi,theta
     */
    public double[] rev(double x, double y) 
	throws PixelBeyondProjectionException
    {
	int   j;
	double cxi, lambda, r, r1, r2, rt, txi, x1, x2, xi;
	double tol = 1.0e-12;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	r = Math.sqrt(x*x + y*y)/r0;
	if (r == 0.0) {
	    xi = 0.0;
	} else if (r < w[4]) {
	    xi = r*w[5];
	} else {

	    // Find a solution interval. 
	    x1 = x2 = 1.0;
	    r1 = r2 = 0.0;
	    for (j = 0; j < 30; j++) {
		x2 = x1/2.0;
		txi = Math.sqrt(1.0-x2*x2)/x2;
		r2 = -(Math.log(x2)/txi + w[0]*txi);

		if (r2 >= r) break;
		x1 = x2;
		r1 = r2;
	    }
	    if (j == 30) 
		throw new PixelBeyondProjectionException(
		            "AIR: No solution interval for (x,y)");

	    cxi=0;
	    for (j = 0; j < 100; j++) {

		// Weighted division of the interval. 
		lambda = (r2-r)/(r2-r1);
		if (lambda < 0.1) {
		    lambda = 0.1;
		} else if (lambda > 0.9) {
		    lambda = 0.9;
		}
		cxi = x2 - lambda*(x2-x1);

		txi = Math.sqrt(1.0-cxi*cxi)/cxi;
		rt = -(Math.log(cxi)/txi + w[0]*txi);

		if (rt < r) {
		    if (r-rt < tol) break;
		    r1 = rt;
		    x1 = cxi;
		} else {
		    if (rt-r < tol) break;
		    r2 = rt;
		    x2 = cxi;
		}
	    }
	    if (j == 100) 
		throw new PixelBeyondProjectionException(
		    "AIR: Weighted division for solution interval not found");

	    xi = TrigD.acos(cxi);
	}

	if (r == 0.0) {
	    out[0] = 0.0;
	} else {
	    out[0] = TrigD.atan2(x, -y);
	}
	out[1] = 90.0 - 2.0*xi;
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
	if (p.length < 2) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "2 projection parameters");
	this.p = null;
	init(r0, p);
    }
}



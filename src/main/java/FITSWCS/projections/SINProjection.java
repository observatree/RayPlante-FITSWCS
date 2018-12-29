/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================
* History:
*  Aug 1996     rlp   Original translation from WCSLIB V2.3
*                     fixed bad placement of } in C code for sinrev 
*  Sept 5 1996  rlp   folded in changes from WCSLIB V2.4:
*                       check for pixel out of range for Orthographic proj.
*===========================================================================*/

package FITSWCS.projections;

import FITSWCS.*;
import FITSWCS.exceptions.*;

/**
 *   This class provides support for the orthographic/synthesis 
 *   projection (SIN) used by the FITS "World Coordinate System" (WCS) 
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
public class SINProjection extends Projection {

    protected final static double tol = 1.0e-13;

    /**
     * Create an SINProjection object (with no obliqueness)
     */
    public SINProjection() {
	double[] use = { 0, 0, 0 };
	init(0, use);
    }

    /**
     * Create an SINProjection object with a given obliqueness
     * @param opa obliqueness parameters alpha
     * @param opb obliqueness parameters beta
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public SINProjection(double opa, double opb) 
    {
	double[] use = { 0, opa, opb };
	init(0, use);
    }

    /**
     * Create an SINProjection object
     * @param r0 sphere radius
     * @param opa obliqueness parameters alpha
     * @param opb obliqueness parameters beta
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public SINProjection(double r0, double opa, double opb) 
    {
	double[] use = { 0, opa, opb };
	init(r0, use);
    }

    /**
     * Create an SINProjection object
     * @param r0 sphere radius
     * @param p  projection parameters specifying the obliqueness, 
     *           where p[1] is the alpha parameter and 
     *           p[2] is the beta parameter
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException if p[1] = -1;
     */
    public SINProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	init(r0, p);
    }

    /**
     * Create an SINProjection object
     * @param p  projection parameters, where p[1] is the 
     *           SIN distance parameter, mu in units of 180/PI.
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     */
    public SINProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	init(0, p);
    }

    /**
     * sets up intermediate data
     * <pre>
     *   r0    reset to 180/pi if 0.
     *   w[0]   1/r0
     *   w[1]   alpha**2 + beta**2
     *   w[2]   2*(alpha**2 + beta**2)
     *   w[3]   2*(alpha**2 + beta**2 + 1)
     *   w[4]   alpha**2 + beta**2 - 1
     * where alpha=p[1] and beta=p[2]
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	if (p == null || p.length < 3) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "3 projection parameters");

	if (r0 == 0.0) r0 = R2D;
	this.r0 = r0;
	this.w = new double[5];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	w[0] = 1.0/r0;
	w[1] = p[1]*p[1] + p[2]*p[2];
	w[2] = 2.0*w[1];
	w[3] = w[2] + 2.0;
	w[4] = w[1] - 1.0;
    }	

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
    {
	double cthe, t, z;
	double[] out = new double[2];

	t = (90.0 - Math.abs(theta))*D2R;
	if (t < 1.0e-5) {
	    if (theta > 0.0) {
		z = -t*t/2.0;
	    } else {
		z = 2.0 - t*t/2.0;
	    }
	    cthe = t;
	} else {
	    z =  TrigD.sin(theta) - 1.0;
	    cthe = TrigD.cos(theta);
	}

	out[0] =  r0*(cthe*TrigD.sin(phi) + p[1]*z);
	out[1] = -r0*(cthe*TrigD.cos(phi) + p[2]*z);
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
	double tol = 1.0e-13;
	double a, b, c, d, r2, sth, sth1, sth2, sxy, x0, xp, y0, yp, z;
	double[] out = new double[2];

	// Compute intermediaries. 
	x0 = x*w[0];
	y0 = y*w[0];
	r2 = x0*x0 + y0*y0;

	if (w[1] == 0.0) {

	    // Orthographic projection. 
	    if (r2 != 0.0) {
		out[0] = TrigD.atan2(x0, -y0);
	    } else {
		out[0] = 0.0;
	    }

	    if (r2 < 0.5) {
		out[1] = TrigD.acos(Math.sqrt(r2));
	    } else if (r2 <= 1.0) {
		out[1] = TrigD.asin(Math.sqrt(1.0 - r2));
	    } else {
		throw new PixelBeyondProjectionException("SIN: (x,y) = (" + x +
						       ", " + y + ")");
	    }

	} else {

	    // "Synthesis" projection. 
	    if (r2 < 1.0e-10) {

		// Use small angle formula. 
		z = -r2/2.0;
		out[1] = 90.0 - R2D*Math.sqrt(r2/(1.0 - x0*p[1] + y0*p[2]));

	    } else {
		sxy = 2.0*(p[1]*x0 - p[2]*y0);

		a = w[3];
		b = -(sxy + w[2]);
		c = r2 + sxy + w[4];
		d = b*b - 2.0*a*c;

		// Check for a solution. 
		if (d < 0.0)  throw new 
			PixelBeyondProjectionException("SIN: (x,y) = (" + x +
						       ", " + y + ")");
		d = Math.sqrt(d);

		// Choose solution closest to pole. 
		sth1 = (-b + d)/a;
		sth2 = (-b - d)/a;
		sth = (sth1>sth2) ? sth1 : sth2;
		if (sth > 1.0) {
		    if (sth-1.0 < tol) {
			sth = 1.0;
		    } else {
			sth = (sth1<sth2) ? sth1 : sth2;
		    }
		}
		if (sth > 1.0 || sth < -1.0) 
		    throw new 
			PixelBeyondProjectionException("SIN: (x,y) = (" + x +
						       ", " + y + ")");

		out[1] = TrigD.asin(sth);
		z = sth - 1.0;
	    }

	    // Compute native coordinates. 
	    xp = -y0 - p[2]*z;
	    yp =  x0 - p[1]*z;
	    if (xp == 0.0 && yp == 0.0) {
		out[0] = 0.0;
	    } else {
		out[0]   = TrigD.atan2(yp,xp);
	    }

	}
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
	    double[] use = { 0.0, 0.0, 0.0 };
	    init(r0, use);
	} 
    }

    /**
     * set the projection parameters 
     * @param p array of projection parameters; if p is null, no update is
     *          to be made.
     */
    public void setProjParm(double[] p)
	throws ArrayIndexOutOfBoundsException
    {
	if (p == null) return;
	if (p.length < 3) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "3 projection parameters");
	init(r0, p);
    }
}



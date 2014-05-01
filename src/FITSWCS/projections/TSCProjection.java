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
*                       Renamed variables l <-> m in the quadcube 
*                         projections to accord with standard usage (and 
*                         revised WCS draft paper).
*                       check for pixel out of range for Orthographic 
*                         projection in rev() method: now throws 
*                         PixelBeyondProjectionException
*===========================================================================*/

package FITSWCS.projections;

import FITSWCS.*;
import FITSWCS.exceptions.*;

/**
 *   This class provides support for the tangential spherical cube
 *   projection (TSC) used by the FITS "World Coordinate System" (WCS) 
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
public class TSCProjection extends Projection {

    /**
     * Create an TSCProjection object
     */
    public TSCProjection() { init(0, new double[1]); }

    /**
     * Create an TSCProjection object
     * @param r0 sphere radius
     */
    public TSCProjection(double r0) { init(r0, new double[1]); }

    /**
     * Create an TSCProjection object
     * @param r0 sphere radius
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public TSCProjection(double r0, double[] p) 
    {
	init(r0, p);
    }

    /**
     * Create an TSCProjection object
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public TSCProjection(double[] p) 
    {
	init(0, p);
    }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   r0*(pi/4)
     *   w[1]   (4/pi)/r0
     * </pre>
     */
    private void init(double r0, double[] p) {
	w = new double[2];

	if (r0 == 0.0) {
	    r0 = R2D;
	    w[0] = 45.0;
	    w[1] = 1.0/45.0;
	} else {
	    w[0] = r0*PI/4.0;
	    w[1] = 1.0/w[0];
	}
	this.r0 = r0;
    }

    protected double tol = 1.0e-12;

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
	throws PixelBeyondProjectionException
    {
	int face;
	double costhe, l, m, n, rho, x0, xf, y0, yf;
	double[] out = new double[2];

	costhe = TrigD.cos(theta);
	l = costhe*TrigD.cos(phi);
	m = costhe*TrigD.sin(phi);
	n = TrigD.sin(theta);

	face = 0;
	rho  = n;
	if (l > rho) {
	    face = 1;
	    rho  = l;
	}
	if (m > rho) {
	    face = 2;
	    rho  = m;
	}
	if (-l > rho) {
	    face = 3;
	    rho  = -l;
	}
	if (-m > rho) {
	    face = 4;
	    rho  = -m;
	}
	if (-n > rho) {
	    face = 5;
	    rho  = -n;
	}

	xf = yf = x0 = y0 = 0.0;
	if (face == 0) {
	    xf =  m/rho;
	    yf = -l/rho;
	    x0 =  0.0;
	    y0 =  2.0;
	} else if (face == 1) {
	    xf =  m/rho;
	    yf =  n/rho;
	    x0 =  0.0;
	    y0 =  0.0;
	} else if (face == 2) {
	    xf = -l/rho;
	    yf =  n/rho;
	    x0 =  2.0;
	    y0 =  0.0;
	} else if (face == 3) {
	    xf = -m/rho;
	    yf =  n/rho;
	    x0 =  4.0;
	    y0 =  0.0;
	} else if (face == 4) {
	    xf =  l/rho;
	    yf =  n/rho;
	    x0 =  6.0;
	    y0 =  0.0;
	} else if (face == 5) {
	    xf =  m/rho;
	    yf =  l/rho;
	    x0 =  0.0;
	    y0 = -2.0;
	}

	if (Math.abs(xf) > 1.0) {
	    if (Math.abs(xf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "TSC: Solution not defined for phi, theta: " + phi + 
		    ", " + theta);

	    xf = (xf < 0) ? -1.0 : 1.0;
	}
	if (Math.abs(yf) > 1.0) {
	    if (Math.abs(yf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "TSC: Solution not defined for phi, theta: " + phi + 
		    ", " + theta);

	    yf = (yf < 0) ? -1.0 : 1.0;
	}

	out[0] = w[0]*(xf + x0);
	out[1] = w[0]*(yf + y0);
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
	int   face;
	boolean direct;
	double l, m, n,  xf, yf;
	double[] out = new double[2];

	xf = x*w[1];
	yf = y*w[1];

//	l = m = n = 0.0;
	// Determine the face.
	if (xf > 7.0) {
	    throw new PixelBeyondProjectionException("x = " + x);
	} else if (xf > 5.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);

	    // face = 4 
	    xf = xf - 6.0;
	    m  = -1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    l  = -m*xf;
	    n  = -m*yf;

	} else if (xf > 3.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);

	    // face = 3 
	    xf = xf - 4.0;
	    l  = -1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    m  =  l*xf;
	    n  = -l*yf;

	} else if (xf > 1.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);

	    // face = 2 
	    xf = xf - 2.0;
	    m  =  1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    l  = -m*xf;
	    n  =  m*yf;

	} else if (xf < -1.0) {
	    throw new PixelBeyondProjectionException("x = " + x);

	} else if (yf > 1.0) {
	    if (yf > 3.0) throw new PixelBeyondProjectionException("y = " + y);

	    // face = 0 
	    yf = yf - 2.0;
	    n  = 1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    l  = -n*yf;
	    m  =  n*xf;

	} else if (yf < -1.0) {
	    if (yf < -3.0) throw new PixelBeyondProjectionException("y = " + y);

	    // face = 5 
	    yf = yf + 2.0;
	    n  = -1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    l  = -n*yf;
	    m  = -n*xf;

	} else {

	    // face = 1 
	    l  =  1.0/Math.sqrt(1.0 + xf*xf + yf*yf);
	    m  =  l*xf;
	    n  =  l*yf;
	}

	if (l == 0.0 && m == 0.0) {
	    out[0] = 0.0;
	} else {
	    out[0] = TrigD.atan2(m, l);
	}
	out[1] = TrigD.asin(n);
	return out;
    }

    /**
     * set the sphere radius 
     */
    public void setR0(double r0) {
	this.r0 = (r0 == 0) ? R2D : r0;
	init(r0, p);
    }

    /**
     * set the projection parameters (which are not used);
     * @param p array of projection parameters; if p is null, no update is
     *          to be made.
     */
    public void setProjParm(double[] p)
    {
	if (p == null) return;
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);
    }
}



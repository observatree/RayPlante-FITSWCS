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
 *   This class provides support for the Bonne's 
 *   projection (BON) used by the FITS "World Coordinate System" (WCS) 
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
public class BONProjection extends Projection {

    protected GLSProjection gls=null;

    /**
     * Create an BONProjection object
     * @param r0     sphere radius; if 0, it defaults to 180/PI 
     * @param theta  Bonne conformal latitude, theta1, in degrees.
     */
    public BONProjection(double r0, double theta) 
    {
	double[] use = { 0, theta };
	init(r0, use);
    }

    /**
     * Create an BONProjection object
     * @param theta  Bonne conformal latitude, theta1, in degrees.
     */
    public BONProjection(double theta) 
    {
	this(0, theta);
    }

    /**
     * Create an BONProjection object
     * @param r0 sphere radius; if 0, it defaults to 180/PI 
     * @param p  projection parameters, where p[1] is the 
     *           Bonne conformal latitude, theta1, in degrees.
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     */
    public BONProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	init(r0, p);
    }

    /**
     * Create an BONProjection object
     * @param p  projection parameters, where p[1] is the 
     *           Bonne conformal latitude, theta1, in degrees.
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     */
    public BONProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	this(0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public BONProjection() { r0 = R2D; }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   r0*pi/180                   
     *   w[1]   Y0 = r0*cot(theta1) + theta1
     * where theta1 is p[1];
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException
    {
	double cxi, r, txi, xi;
	if (p == null || p.length < 2) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "2 projection parameters");

	if (p[1] == 0.0) {
	    gls = new GLSProjection(p);
	}
	else {
	    this.w = new double[3];
	    this.p = new double[p.length];
	    System.arraycopy(p, 0, this.p, 0, p.length);

	    if (r0 == 0.0) {
		r0 = R2D;
		w[1] = 1.0;
		w[2] = r0*TrigD.cos(p[1])/TrigD.sin(p[1]) + p[1];
	    } else {
		w[1] = r0*D2R;
		w[2] = r0*(TrigD.cos(p[1])/TrigD.sin(p[1]) + p[1]*D2R);
	    }
	    this.r0 = r0;
	}	
    }

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
    {
	double a, r;
	double[] out;
	if (p == null) throw new UnsetProjectionParameterException();

	if (gls != null) {

	    // Sanson-Flamsteed. 
	    out = gls.fwd(phi, theta);
	}
	else {
	    out = new double[2];

	    r = w[2] - theta*w[1];
	    a = r0*phi*TrigD.cos(theta)/r;

	    out[0] =        r*TrigD.sin(a);  // x
	    out[1] = w[2] - r*TrigD.cos(a);  // y
	}
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
    {
	double a, dy, costhe, r;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	if (gls != null) {

	    // Sanson-Flamsteed. 
	    out = gls.rev(x, y);
	}
	else {

	    dy = w[2] - y;
	    r = Math.sqrt(x*x + dy*dy);
	    if (p[1] < 0.0) r = -r;

	    if (r == 0.0) {
		a = 0.0;
	    } else {
		a = TrigD.atan2(x/r, dy/r);
	    }

	    out[1] = (w[2] - r)/w[1];
	    costhe = TrigD.cos(out[1]);
	    if (costhe == 0.0) {
		out[0] = 0.0;
	    } else {
		out[0] = a*(r/r0)/TrigD.cos(out[1]);
	    }
	}

	return out;
    }

    /**
     * same as rev(xy[0], xy[1]);
     */
    public double[] rev(double[] xy) { return rev(xy[0], xy[1]); }

    /**
     * set the sphere radius 
     */
    public void setR0(double r0) {
	this.r0 = (r0 == 0) ? R2D : r0;
	try {
	    if (p != null) init(r0, p);
	} catch (ArrayIndexOutOfBoundsException ex) {
	    p = null;
	} 
    }

    /**
     * set the projection parameters 
     * @param p array of projection parameters; a null value cause no updates
     *          to be made.
     * @exception ArrayIndexOutOfBoundsException if p.length < 2 (parameters
     *          left unchanged);
     * @exception BadProjectionParameterException if p[1] = -1 
     *          (parameters left in invalid state);
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



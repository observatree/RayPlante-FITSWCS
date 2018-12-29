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
 *   This class provides support for the conic equal area
 *   projection (COE) used by the FITS "World Coordinate System" (WCS) 
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
public class COEProjection extends Projection {

    protected double tol = 1.0e-12;

    /**
     * Create an COEProjection object with a given oblateness
     * @param sigma  the average of the latitudes of the standard parallels,
     *               in degrees
     * @param delta  half the difference in the latitudes of the standard 
     *               parallels, in degrees
     * @exception BadProjectionParameterException 
     */
    public COEProjection(double sigma, double delta) 
	throws BadProjectionParameterException
    {
	double[] use = { 0, sigma, delta };
	init(0, use);
    }

    /**
     * Create an COEProjection object
     * @param r0     sphere radius
     * @param sigma  the average of the latitudes of the standard parallels,
     *               in degrees
     * @param delta  half the difference in the latitudes of the standard 
     *               parallels, in degrees
     * @exception BadProjectionParameterException 
     */
    public COEProjection(double r0, double sigma, double delta) 
	throws BadProjectionParameterException
    {
	double[] use = { 0, sigma, delta };
	init(r0, use);
    }

    /**
     * Create an COEProjection object
     * @param r0 sphere radius
     * @param p  projection parameters, where p[1] the average and p[2] 
     *           is the difference in the standard latitudes
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException 
     */
    public COEProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(r0, p);
    }

    /**
     * Create an COEProjection object
     * @param p  projection parameters, where p[1] is mu, the distance of 
     *           point of projection from the centre of the generating 
     *           sphere, and p[2] is lambda, the radius of the cylinder of 
     *           projection, both in units of r0
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException 
     */
    public COEProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public COEProjection() { r0 = R2D; }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   C = (sin(theta1) + sin(theta2))/2
     *   w[1]   1/C
     *   w[2]   Y0 = chi*Math.sqrt(psi - 2C*TrigD.sin(sigma))
     *   w[3]   chi = r0/C
     *   w[4]   psi = 1 + sin(theta1)*sin(theta2)
     *   w[5]   2C
     *   w[6]   (1 + sin(theta1)*sin(theta2))*(r0/C)**2
     *   w[7]   C/(2*r0**2)
     *   w[8]   chi*Math.sqrt(psi + 2C)
     * where sigma is p[1] and delta is p[2], and theta1 = sigma + delta 
     * and theta2 = sigma - delta.
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	double theta1, theta2;

	if (p == null || p.length < 3) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "3 projection parameters");

	if (r0 == 0.0) r0 = R2D;
	this.r0 = r0;
	this.w = new double[9];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	theta1 = p[1] - p[2];
	theta2 = p[1] + p[2];

	w[0] = (TrigD.sin(theta1) + TrigD.sin(theta2))/2.0;
	if (w[0] == 0.0) throw new 
	    BadProjectionParameterException(
               "Projection parameters: sin(theta1) + sin(theta2) = 0");

	w[1] = 1.0/w[0];

	w[3] = r0/w[0];
	w[4] = 1.0 + TrigD.sin(theta1)*TrigD.sin(theta2);
	w[5] = 2.0*w[0];
	w[6] = w[3]*w[3]*w[4];
	w[7] = 1.0/(2.0*r0*w[3]);
	w[8] = w[3]*Math.sqrt(w[4] + w[5]);

	w[2] = w[3]*Math.sqrt(w[4] - w[5]*TrigD.sin(p[1]));
    }	

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
    {
	double a, r;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	a = phi*w[0];
	if (theta == -90.0) {
	    r = w[8];
	} else {
	    r = w[3]*Math.sqrt(w[4] - w[5]*TrigD.sin(theta));
	}

	out[0] =        r*TrigD.sin(a);       
	out[1] = w[2] - r*TrigD.cos(a);
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
	double a, dy, r, ww;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	dy = w[2] - y;
	r  = Math.sqrt(x*x + dy*dy);
	if (p[1] < 0.0) r = -r;

	if (r == 0.0) {
	    a = 0.0;
	} else {
	    a = TrigD.atan2(x/r, dy/r);
	}

	out[0] = a*w[1];                          

	if (Math.abs(r - w[8]) < tol) {
	    out[1] = -90.0;
	} else {
	    ww = (w[6] - r*r)*w[7];
	    if (Math.abs(ww) > 1.0) {
		if (Math.abs(ww-1.0) < tol) {
		    out[1] = 90.0;
		} else if (Math.abs(ww+1.0) < tol) {
		    out[1] = -90.0;
		} else {
		    throw new PixelBeyondProjectionException(
		 "COE: Calculation failed for (x,y) = (" + x + ", " + y + ")");
		}
	    } else {
		out[1] = TrigD.asin(ww);
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
	    p = null;
	} catch (BadProjectionParameterException ex) {
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
	if (p.length < 3) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "3 projection parameters");
	this.p = null;
	init(r0, p);
    }
}



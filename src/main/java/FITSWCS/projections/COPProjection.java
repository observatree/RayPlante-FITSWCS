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
 *   This class provides support for the conic perspective
 *   projection (COP) used by the FITS "World Coordinate System" (WCS) 
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
public class COPProjection extends Projection {

    /**
     * Create an COPProjection object with a given oblateness
     * @param sigma  the average of the latitudes of the standard parallels,
     *               in degrees
     * @param delta  half the difference in the latitudes of the standard 
     *               parallels, in degrees
     * @exception BadProjectionParameterException if mu = 0 or if 
     *               mu + lambda = 0;
     */
    public COPProjection(double sigma, double delta) 
	throws BadProjectionParameterException
    {
	double[] use = { 0, sigma, delta };
	init(0, use);
    }

    /**
     * Create an COPProjection object
     * @param r0     sphere radius
     * @param sigma  the average of the latitudes of the standard parallels,
     *               in degrees
     * @param delta  half the difference in the latitudes of the standard 
     *               parallels, in degrees
     * @exception BadProjectionParameterException if mu = 0 or if 
     *               mu + lambda = 0;
     */
    public COPProjection(double r0, double sigma, double delta) 
	throws BadProjectionParameterException
    {
	double[] use = { 0, sigma, delta };
	init(r0, use);
    }

    /**
     * Create an COPProjection object
     * @param r0 sphere radius
     * @param p  projection parameters, where p[1] the average and p[2] 
     *           is the difference in the standard latitudes
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0;
     */
    public COPProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(r0, p);
    }

    /**
     * Create an COPProjection object
     * @param p  projection parameters, where p[1] is mu, the distance of 
     *           point of projection from the centre of the generating 
     *           sphere, and p[2] is lambda, the radius of the cylinder of 
     *           projection, both in units of r0
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0;
     */
    public COPProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public COPProjection() { r0 = R2D; }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   C  = sin(sigma)
     *   w[1]   1/C
     *   w[2]   Y0 = r0*cos(delta)*cot(sigma)
     *   w[3]   r0*cos(delta)
     *   w[4]   1/(r0*cos(delta)
     *   w[5]   cot(sigma)
     * where sigma is p[1] and delta is p[2]
     * </pre>
     */
    private void init(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	if (p == null || p.length < 3) 
          throw new ArrayIndexOutOfBoundsException("Need at least " +
		 		                   "3 projection parameters");
	if (p[2] == 0.0) 
	  throw new BadProjectionParameterException("p[2] must be non-zero");
	if (p[1] == -p[2])
	  throw new BadProjectionParameterException("p[1] + p[2] = 0");

	if (r0 == 0.0) r0 = R2D;
	this.r0 = r0;
	this.w = new double[6];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	w[0] = TrigD.sin(p[1]);
	if (w[0] == 0.0) throw new 
	    BadProjectionParameterException("Bad value for sigma: " + p[1]);

	w[1] = 1.0/w[0];

	w[3] = r0*TrigD.cos(p[2]);
	if (w[0] == 0.0) throw new 
	    BadProjectionParameterException("Bad value for delta: " + p[2]);

	w[4] = 1.0/w[3];
	w[5] = 1.0/TrigD.tan(p[1]);
	w[2] = w[3]*w[5];
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

	a = w[0]*phi;
	r = w[2] - w[3]*TrigD.tan(theta-p[1]);

	out[0] = r*TrigD.sin(a);       
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
	double a, dy, r;
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
	out[1] = p[1] + TrigD.atan(w[5] - r*w[4]);
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
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0; 
     *               (parameters left in invalid state)
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



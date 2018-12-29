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
 *   This class provides support for the Hammer-Aitoff projection (AIT) 
 *   used by the FITS "World Coordinate System" (WCS) convention. <p>
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
public class AITProjection extends Projection {

    protected double tol = 1.0e-12;

    /**
     * Create an AITProjection object
     */
    public AITProjection() { init(0, new double[1]); }

    /**
     * Create an AITProjection object
     * @param r0 sphere radius
     */
    public AITProjection(double r0) { init(r0, new double[1]); }

    /**
     * Create an AITProjection object
     * @param r0 sphere radius
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public AITProjection(double r0, double[] p) 
    {
	init(r0, p);
    }

    /**
     * Create an AITProjection object
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public AITProjection(double[] p) 
    {
	init(0, p);
    }

    /**
     * sets up intermediate data
     * <pre>
     *   r0     reset to 180/pi if 0.
     *   w[0]   2*r0**2
     *   w[1]   1/(2*r0)**2
     *   w[2]   1/(4*r0)**2
     *   w[3]   1/(2*r0)
     * </pre>
     */
    private void init(double r0, double[] p) {
	w = new double[4];

	if (r0 == 0.0) r0 = R2D;

	w[0] = 2.0*r0*r0;
	w[1] = 1.0/(2.0*w[0]);
	w[2] = w[1]/4.0;
	w[3] = 1.0/(2.0*r0);
	this.r0 = r0;
    }

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
    {
	double costhe, ww;
	double[] out = new double[2];

	costhe = TrigD.cos(theta);
	ww = Math.sqrt(w[0]/(1.0 + costhe*TrigD.cos(phi/2.0)));
	out[0] = 2.0*ww*costhe*TrigD.sin(phi/2.0);
	out[1] = ww*TrigD.sin(theta);
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
	double s, u, xp, yp, z;
	double[] out = new double[2];

	u = 1.0 - x*x*w[2] - y*y*w[1];
	if (u < 0.0) throw new PixelBeyondProjectionException(
		    "AIT: Solution not defined for x,y: " + x + ", " + y);

	z = Math.sqrt(u);
	s = z*y/r0;
	if (s < -1.0 || s > 1.0) throw new PixelBeyondProjectionException(
		    "AIT: Solution not defined for x,y: " + x + ", " + y);

	xp = 2.0*z*z - 1.0;
	yp = z*x*w[3];
	if (xp == 0.0 && yp == 0.0) {
	    out[0] = 0.0;
	} else {
	    out[0] = 2.0*TrigD.atan2(yp, xp);
	}
	out[1] = TrigD.asin(s);

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

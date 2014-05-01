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
 *   This class provides support for the cylindrical perspective
 *   projection (CYP) used by the FITS "World Coordinate System" (WCS) 
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
public class CYPProjection extends Projection {

    /**
     * Create an CYPProjection object with a given oblateness
     * @param mu     Distance of point of projection from the centre of the
     *               generating sphere, in units of 180/PI.
     * @param lambda Radius of the cylinder of projection, in units of 180/PI
     * @exception BadProjectionParameterException if mu = 0 or if 
     *               mu + lambda = 0;
     */
    public CYPProjection(double mu, double lambda) 
	throws BadProjectionParameterException
    {
	if (mu == 0.0) 
	  throw new BadProjectionParameterException("mu must be non-zero");
	if (lambda == -mu)
	  throw new BadProjectionParameterException("mu + lambda = 0");
	
	double[] use = { 0, mu, lambda };
	init(0, use);
    }

    /**
     * Create an CYPProjection object
     * @param r0     sphere radius
     * @param mu     Distance of point of projection from the centre of the
     *               generating sphere, in units of 180/PI.
     * @param lambda Radius of the cylinder of projection, in units of 180/PI
     * @exception BadProjectionParameterException if mu = 0 or if 
     *               mu + lambda = 0;
     */
    public CYPProjection(double r0, double mu, double lambda) 
	throws BadProjectionParameterException
    {
	if (mu == 0.0) 
	  throw new BadProjectionParameterException("mu must be non-zero");
	if (lambda == -mu)
	  throw new BadProjectionParameterException("mu + lambda = 0");
	
	double[] use = { 0, mu, lambda };
	init(r0, use);
    }

    /**
     * Create an CYPProjection object
     * @param r0 sphere radius
     * @param p  projection parameters, where p[1] is mu, the distance of 
     *           point of projection from the centre of the generating 
     *           sphere, and p[2] is lambda, the radius of the cylinder of 
     *           projection, both in units of r0
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0;
     */
    public CYPProjection(double r0, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(r0, p);
    }

    /**
     * Create an empty Projection object; a call to setProjParm() must
     * be made to avoid throwing UnsetProjectionParameterException
     */
    public CYPProjection() { r0 = R2D; }

    /**
     * Create an CYPProjection object
     * @param p  projection parameters, where p[1] is mu, the distance of 
     *           point of projection from the centre of the generating 
     *           sphere, and p[2] is lambda, the radius of the cylinder of 
     *           projection, both in units of r0
     * @exception ArrayIndexOutOfBoundsException if p.length < 2;
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0;
     */
    public CYPProjection(double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException
    {
	init(0, p);
    }

    /**
     * sets up intermediate data
     * <pre>
     *   r0    reset to 180/pi if 0.
     *   w[0]   r0*lambda*(pi/180)  
     *   w[1]   (180/pi)/(r0*lambda)
     *   w[2]   r0*(mu + lambda)    
     *   w[3]   1/(r0*(mu + lambda))
     * where lambda is p[2], the radius of the cylinder of projection,
     * in units of r0
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
	this.w = new double[5];
	this.p = new double[p.length];
	System.arraycopy(p, 0, this.p, 0, p.length);

	w[0] = r0*p[2]*D2R;
	w[1] = 1.0/w[0];
	w[2] = r0*(p[1] + p[2]);
	w[3] = 1.0/w[2];
    }	

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
	throws PixelBeyondProjectionException
    {
	double s;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	s = p[1] + TrigD.cos(theta);
	if (s == 0.0) 
	    throw new 
		PixelBeyondProjectionException("CYP: theta out of bounds: " +
					       theta);

	out[0] = w[0]*phi;
	out[1] = w[2]*TrigD.sin(theta)/s;
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
	double eta;
	double[] out = new double[2];
	if (p == null) throw new UnsetProjectionParameterException();

	out[0] = x*w[1];
	eta    = y*w[3];
	out[1] = TrigD.atan2(eta,1.0) + 
	    TrigD.asin(eta*p[1]/Math.sqrt(eta*eta+1.0));

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
     * @exception ArrayIndexOutOfBoundsException if p.length < 3;
     * @exception BadProjectionParameterException if p[1] = 0 or if 
     *               p[1] + p[2] = 0;
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



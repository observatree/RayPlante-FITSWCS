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
 *   This class provides support for the COBE quadrilateralized spherical 
 *   cube projection (CSC) used by the FITS "World Coordinate System" (WCS) 
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
public class CSCProjection extends Projection {

    /**
     * Create an CSCProjection object
     */
    public CSCProjection() { init(0, new double[1]); }

    /**
     * Create an CSCProjection object
     * @param r0 sphere radius
     */
    public CSCProjection(double r0) { init(r0, new double[1]); }

    /**
     * Create an CSCProjection object
     * @param r0 sphere radius
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public CSCProjection(double r0, double[] p) 
    {
	init(r0, p);
    }

    /**
     * Create an CSCProjection object
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public CSCProjection(double[] p) 
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

    public final static float gstar  =  1.37484847732f; 
    public final static float mm     =  0.004869491981f;
    public final static float gamma  = -0.13161671474f; 
    public final static float omega1 = -0.159596235474f;
    public final static float d0  =  0.0759196200467f;  
    public final static float d1  = -0.0217762490699f;  
    public final static float c00 =  0.141189631152f;   
    public final static float c10 =  0.0809701286525f;  
    public final static float c01 = -0.281528535557f;   
    public final static float c11 =  0.15384112876f;    
    public final static float c20 = -0.178251207466f;   
    public final static float c02 =  0.106959469314f;   

    protected float tol = 1.0e-7f;

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). 
     * @return double[] a two-element array containing x,y
     */
    public double[] fwd(double phi, double theta) 
	throws PixelBeyondProjectionException
    {
	int face;
	double costhe, eta, l, m, n, rho, xi;
	float a, a2, a2b2, a4, ab, b, b2, b4, ca2, cb2, x0, xf, y0, yf;
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

	xi = eta = 0.0;
	x0 = y0 = 0.0f;
	if (face == 0) {
	    xi  =  m;
	    eta = -l;
	    x0  =  0.0f;
	    y0  =  2.0f;
	} else if (face == 1) {
	    xi  =  m;
	    eta =  n;
	    x0  =  0.0f;
	    y0  =  0.0f;
	} else if (face == 2) {
	    xi  = -l;
	    eta =  n;
	    x0  =  2.0f;
	    y0  =  0.0f;
	} else if (face == 3) {
	    xi  = -m;
	    eta =  n;
	    x0  =  4.0f;
	    y0  =  0.0f;
	} else if (face == 4) {
	    xi  =  l;
	    eta =  n;
	    x0  =  6.0f;
	    y0  =  0.0f;
	} else if (face == 5) {
	    xi  =  m;
	    eta =  l;
	    x0  =  0.0f;
	    y0  = -2.0f;
	}

	a = (float)  (xi/rho);
	b = (float) (eta/rho);

	a2 = a*a;
	b2 = b*b;
	ca2 = 1.0f - a2;
	cb2 = 1.0f - b2;

	// Avoid floating underflows. 
	ab   = Math.abs(a*b);
	a4   = (a2 > 1.0e-16) ? a2*a2 : 0.0f;
	b4   = (b2 > 1.0e-16) ? b2*b2 : 0.0f;
	a2b2 = (ab > 1.0e-16) ? a2*b2 : 0.0f;

	xf = a*(a2 + ca2*(gstar + b2*(gamma*ca2 + mm*a2 +
          cb2*(c00 + c10*a2 + c01*b2 + c11*a2b2 + c20*a4 + c02*b4)) +
          a2*(omega1 - ca2*(d0 + d1*a2))));
	yf = b*(b2 + cb2*(gstar + a2*(gamma*cb2 + mm*b2 +
          ca2*(c00 + c10*b2 + c01*a2 + c11*a2b2 + c20*b4 + c02*a4)) +
          b2*(omega1 - cb2*(d0 + d1*b2))));

	if (Math.abs(xf) > 1.0) {
	    if (Math.abs(xf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "CSC: Solution not defined for phi, theta: " + phi + 
		    ", " + theta);

	    xf = (xf < 0) ? -1.0f : 1.0f;
	}
	if (Math.abs(yf) > 1.0) {
	    if (Math.abs(yf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "CSC: Solution not defined for phi, theta: " + phi + 
		    ", " + theta);

	    yf = (yf < 0) ? -1.0f : 1.0f;
	}

	out[0] = w[0]*(x0 + xf);
	out[1] = w[0]*(y0 + yf);
	return out;
    }

    public static final float p00 = -0.27292696f;
    public static final float p10 = -0.07629969f;
    public static final float p20 = -0.22797056f;
    public static final float p30 =  0.54852384f;
    public static final float p40 = -0.62930065f;
    public static final float p50 =  0.25795794f;
    public static final float p60 =  0.02584375f;
    public static final float p01 = -0.02819452f;
    public static final float p11 = -0.01471565f;
    public static final float p21 =  0.48051509f;
    public static final float p31 = -1.74114454f;
    public static final float p41 =  1.71547508f;
    public static final float p51 = -0.53022337f;
    public static final float p02 =  0.27058160f;
    public static final float p12 = -0.56800938f;
    public static final float p22 =  0.30803317f;
    public static final float p32 =  0.98938102f;
    public static final float p42 = -0.83180469f;
    public static final float p03 = -0.60441560f;
    public static final float p13 =  1.50880086f;
    public static final float p23 = -0.93678576f;
    public static final float p33 =  0.08693841f;
    public static final float p04 =  0.93412077f;
    public static final float p14 = -1.41601920f;
    public static final float p24 =  0.33887446f;
    public static final float p05 = -0.63915306f;
    public static final float p15 =  0.52032238f;
    public static final float p06 =  0.14381585f;

    /**
     * Compute native spherical coordinates (phi,theta) from the 
     * (x,y) coordinates in the plane of projection. 
     * @return double[] a two-element array containing phi,theta
     */
    public double[] rev(double x, double y) 
	throws PixelBeyondProjectionException
    {
	int   face;
	double l, m, n;
	float     a, b, xf, xx, yf, yy, z0, z1, z2, z3, z4, z5, z6;
	double[] out = new double[2];

	xf = (float) (x*w[1]);
	yf = (float) (y*w[1]);

	// Determine the face.
	if (xf > 7.0) {
	    throw new PixelBeyondProjectionException("x = " + x);
	} else if (xf > 5.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 4;
	    xf = xf - 6.0f;
	} else if (xf > 3.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 3;
	    xf = xf - 4.0f;
	} else if (xf > 1.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 2;
	    xf = xf - 2.0f;
	} else if (xf < -1.0) {
	    throw new PixelBeyondProjectionException("x = " + x);
	} else if (yf > 1.0) {
	    if (yf > 3.0) throw new PixelBeyondProjectionException("y = " + y);
	    face = 0;
	    yf = yf - 2.0f;
	} else if (yf < -1.0) {
	    if (yf < -3.0) throw new PixelBeyondProjectionException("y = " + y);
	    face = 5;
	    yf = yf + 2.0f;
	} else {
	    face = 1;
	}

	xx  =  xf*xf;
	yy  =  yf*yf;

	z0 = p00 + 
	       xx*(p10 + xx*(p20 + xx*(p30 + xx*(p40 + xx*(p50 + xx*(p60))))));
	z1 = p01 + xx*(p11 + xx*(p21 + xx*(p31 + xx*(p41 + xx*(p51)))));
	z2 = p02 + xx*(p12 + xx*(p22 + xx*(p32 + xx*(p42))));
	z3 = p03 + xx*(p13 + xx*(p23 + xx*(p33)));
	z4 = p04 + xx*(p14 + xx*(p24));
	z5 = p05 + xx*(p15);
	z6 = p06;

	a = z0 + yy*(z1 + yy*(z2 + yy*(z3 + yy*(z4 + yy*(z5 + yy*z6)))));
	a = xf + xf*(1.0f - xx)*a;

	z0 = p00 + 
	       yy*(p10 + yy*(p20 + yy*(p30 + yy*(p40 + yy*(p50 + yy*(p60))))));
	z1 = p01 + yy*(p11 + yy*(p21 + yy*(p31 + yy*(p41 + yy*(p51)))));
	z2 = p02 + yy*(p12 + yy*(p22 + yy*(p32 + yy*(p42))));
	z3 = p03 + yy*(p13 + yy*(p23 + yy*(p33)));
	z4 = p04 + yy*(p14 + yy*(p24));
	z5 = p05 + yy*(p15);
	z6 = p06;

	b = z0 + xx*(z1 + xx*(z2 + xx*(z3 + xx*(z4 + xx*(z5 + xx*z6)))));
	b = yf + yf*(1.0f - yy)*b;

	l = m = n = 0.0;
	if (face == 0) {
	    n =  1.0/Math.sqrt(a*a + b*b + 1.0);
	    l = -b*n;
	    m =  a*n;
	} else if (face == 1) {
	    l =  1.0/Math.sqrt(a*a + b*b + 1.0);
	    m =  a*l;
	    n =  b*l;
	} else if (face == 2) {
	    m =  1.0/Math.sqrt(a*a + b*b + 1.0);
	    l = -a*m;
	    n =  b*m;
	} else if (face == 3) {
	    l = -1.0/Math.sqrt(a*a + b*b + 1.0);
	    m =  a*l;
	    n = -b*l;
	} else if (face == 4) {
	    m = -1.0/Math.sqrt(a*a + b*b + 1.0);
	    l = -a*m;
	    n = -b*m;
	} else if (face == 5) {
	    n = -1.0/Math.sqrt(a*a + b*b + 1.0);
	    l = -b*n;
	    m = -a*n;
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



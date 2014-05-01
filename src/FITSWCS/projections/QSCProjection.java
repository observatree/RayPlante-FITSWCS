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
*                         projection in rev() method
*===========================================================================*/

package FITSWCS.projections;

import FITSWCS.*;
import FITSWCS.exceptions.*;

/**
 *   This class provides support for the quadrilateralized spherical cube
 *   projection (QSC) used by the FITS "World Coordinate System" (WCS) 
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
public class QSCProjection extends Projection {

    /**
     * Create an QSCProjection object
     */
    public QSCProjection() { init(0, new double[1]); }

    /**
     * Create an QSCProjection object
     * @param r0 sphere radius
     */
    public QSCProjection(double r0) { init(r0, new double[1]); }

    /**
     * Create an QSCProjection object
     * @param r0 sphere radius
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public QSCProjection(double r0, double[] p) 
    {
	init(r0, p);
    }

    /**
     * Create an QSCProjection object
     * @param p  projection parameters.  None are actually used by 
     *           this projection; this constructor is provided to 
     *           parallel other subclasses of the Projection class
     */
    public QSCProjection(double[] p) 
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
	double chi, costhe, eta, l, m, n, p, psi, rho, rhu, t, x0, 
	       xf, xi, y0, yf;
	double[] out = new double[2];

	if (Math.abs(theta) == 90.0) {
	    out[0] = 0.0;
	    out[1] = (theta < 0.0) ? -Math.abs(2.0*w[0]) : Math.abs(2.0*w[0]);
	    return out;
	}

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

	rhu = 1.0 - rho;

	xi = eta = x0 = y0 = 0.0;
	if (face == 0) {
	    xi  =  m;
	    eta = -l;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = (90.0 - theta)*D2R;
		rhu = t*t/2.0;
	    }
	    x0  =  0.0;
	    y0  =  2.0;
	} else if (face == 1) {
	    xi  =  m;
	    eta =  n;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = theta*D2R;
		p = Math.IEEEremainder(phi,360.0);
		if (p < -180.0) p += 360.0;
		if (p >  180.0) p -= 360.0;
		p *= D2R;
		rhu = (p*p + t*t)/2.0;
	    }
	    x0  =  0.0;
	    y0  =  0.0;
	} else if (face == 2) {
	    xi  = -l;
	    eta =  n;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = theta*D2R;
		p = Math.IEEEremainder(phi,360.0);
		if (p < -180.0) p += 360.0;
		p = (90.0 - p)*D2R;
		rhu = (p*p + t*t)/2.0;
	    }
	    x0  =  2.0;
	    y0  =  0.0;
	} else if (face == 3) {
	    xi  = -m;
	    eta =  n;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = theta*D2R;
		p = Math.IEEEremainder(phi,360.0);
		if (p < 0.0) p += 360.0;
		p = (180.0 - p)*D2R;
		rhu = (p*p + t*t)/2.0;
	    }
	    x0  =  4.0;
	    y0  =  0.0;
	} else if (face == 4) {
	    xi  =  l;
	    eta =  n;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = theta*D2R;
		p = Math.IEEEremainder(phi,360.0);
		if (p > 180.0) p -= 360.0;
		p *= (90.0 + p)*D2R;
		rhu = (p*p + t*t)/2.0;
	    }
	    x0  =  6;
	    y0  =  0.0;
	} else if (face == 5) {
	    xi  =  m;
	    eta =  l;
	    if (rhu < 1.0e-8) {

		// Small angle formula. 
		t = (90.0 + theta)*D2R;
		rhu = t*t/2.0;
	    }
	    x0  =  0.0;
	    y0  = -2;
	}

	xf = yf = 0.0;
	if (xi == 0.0 && eta == 0.0) {
	    xf  = 0.0;
	    yf  = 0.0;
	} else if (-xi >= Math.abs(eta)) {
	    psi = eta/xi;
	    chi = 1.0 + psi*psi;
	    xf  = -Math.sqrt(rhu/(1.0-1.0/Math.sqrt(1.0+chi)));
	    yf  = (xf/15.0)*(TrigD.atan(psi) - 
			     TrigD.asin(psi/Math.sqrt(chi+chi)));
	} else if (xi >= Math.abs(eta)) {
	    psi = eta/xi;
	    chi = 1.0 + psi*psi;
	    xf  =  Math.sqrt(rhu/(1.0-1.0/Math.sqrt(1.0+chi)));
	    yf  = (xf/15.0)*(TrigD.atan(psi) - 
			     TrigD.asin(psi/Math.sqrt(chi+chi)));
	} else if (-eta > Math.abs(xi)) {
	    psi = xi/eta;
	    chi = 1.0 + psi*psi;
	    yf  = -Math.sqrt(rhu/(1.0-1.0/Math.sqrt(1.0+chi)));
	    xf  = (yf/15.0)*(TrigD.atan(psi) - 
			     TrigD.asin(psi/Math.sqrt(chi+chi)));
	} else if (eta > Math.abs(xi)) {
	    psi = xi/eta;
	    chi = 1.0 + psi*psi;
	    yf  =  Math.sqrt(rhu/(1.0-1.0/Math.sqrt(1.0+chi)));
	    xf  = (yf/15.0)*(TrigD.atan(psi) - 
			     TrigD.asin(psi/Math.sqrt(chi+chi)));
	}

	if (Math.abs(xf) > 1.0) {
	    if (Math.abs(xf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "QSC: Solution not defined for phi, theta: " + phi + 
		    ", " + theta);

	    xf = (xf < 0) ? -1.0 : 1.0;
	}
	if (Math.abs(yf) > 1.0) {
	    if (Math.abs(yf) > 1.0+tol) 
		throw new PixelBeyondProjectionException(
		    "QSC: Solution not defined for phi, theta: " + phi + 
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
	double chi, l, m, n, psi, rho, rhu, xf, yf, ww;
	double[] out = new double[2];

	xf = x*w[1];
	yf = y*w[1];

	// Determine the face.
	if (xf > 7.0) {
	    throw new PixelBeyondProjectionException("x = " + x);
	} else if (xf > 5.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 4;
	    xf = xf - 6.0;
	} else if (xf > 3.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 3;
	    xf = xf - 4.0;
	} else if (xf > 1.0) {
	    if (Math.abs(yf) > 1.0) 
		throw new PixelBeyondProjectionException("y = " + y);
	    face = 2;
	    xf = xf - 2.0;
	} else if (xf < -1.0) {
	    throw new PixelBeyondProjectionException("x = " + x);
	} else if (yf > 1.0) {
	    if (yf > 3.0) throw new PixelBeyondProjectionException("y = " + y);
	    face = 0;
	    yf = yf - 2.0;
	} else if (yf < -1.0) {
	    if (yf < -3.0) throw new PixelBeyondProjectionException("y = " + y);
	    face = 5;
	    yf = yf + 2.0;
	} else {
	    face = 1;
	}

	direct = (Math.abs(xf) > Math.abs(yf));
	if (direct) {
	    if (xf == 0.0) {
		psi = 0.0;
		chi = 1.0;
		rho = 1.0;
		rhu = 0.0;
	    } else {
		ww = 15.0*yf/xf;
		psi = TrigD.sin(ww)/(TrigD.cos(ww) - SQRT2INV);
		chi = 1.0 + psi*psi;
		rhu = xf*xf*(1.0 - 1.0/Math.sqrt(1.0 + chi));
		rho = 1.0 - rhu;
	    }
	} else {
	    if (yf == 0.0) {
		psi = 0.0;
		chi = 1.0;
		rho = 1.0;
		rhu = 0.0;
	    } else {
		ww = 15.0*xf/yf;
		psi = TrigD.sin(ww)/(TrigD.cos(ww) - SQRT2INV);
		chi = 1.0 + psi*psi;
		rhu = yf*yf*(1.0 - 1.0/Math.sqrt(1.0 + chi));
		rho = 1.0 - rhu;
	    }
	}

	if (rho < -1.0) {
	    if (rho < -1.0-tol) 
		throw new PixelBeyondProjectionException(
		    "QSC: Solution not defined for x, y: " + x + ", " + y);

	    rho = -1.0;
	    rhu =  2.0;
	    ww   =  0.0;
	} else {
	    ww = Math.sqrt(rhu*(2.0-rhu)/chi);
	}

	l = m = n = 0.0;
	if (face == 0) {
	    n = rho;
	    if (direct) {
		m = ww;
		if (xf < 0.0) m = -m;
		l = -m*psi;
	    } else {
		l = ww;
		if (yf > 0.0) l = -l;
		m = -l*psi;
	    }
	} else if (face == 1) {
	    l = rho;
	    if (direct) {
		m = ww;
		if (xf < 0.0) m = -m;
		n = m*psi;
	    } else {
		n = ww;
		if (yf < 0.0) n = -n;
		m = n*psi;
	    }
	} else if (face == 2) {
	    m =  rho;
	    if (direct) {
		l = ww;
		if (xf > 0.0) l = -l;
		n = -l*psi;
	    } else {
		n = ww;
		if (yf < 0.0) n = -n;
		l = -n*psi;
	    }
	} else if (face == 3) {
	    l = -rho;
	    if (direct) {
		m = ww;
		if (xf > 0.0) m = -m;
		n = -m*psi;
	    } else {
		n = ww;
		if (yf < 0.0) n = -n;
		m = -n*psi;
	    }
	} else if (face == 4) {
	    m = -rho;
	    if (direct) {
		l = ww;
		if (xf < 0.0) l = -l;
		n = l*psi;
	    } else {
		n = ww;
		if (yf < 0.0) n = -n;
		l = n*psi;
	    }
	} else if (face == 5) {
	    n = -rho;
	    if (direct) {
		m = ww;
		if (xf < 0.0) m = -m;
		l = m*psi;
	    } else {
		l = ww;
		if (yf < 0.0) l = -l;
		m = l*psi;
	    }
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



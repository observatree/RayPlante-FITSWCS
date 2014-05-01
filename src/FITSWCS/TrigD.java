/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS;

/**
 *   This class provides trigonometric or inverse trigonometric
 *   functions which take or return angular arguments in decimal degrees. 
 *   <p>
 *
 *   These routines explicitly handle angles which are a multiple of 
 *   90 degrees returning an exact result.
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
public class TrigD {

    public final static double PI = Math.PI;
    public final static double d2r = PI / 180.0;
    public final static double r2d = 180.0 / PI;
    private final static double tol = 1e-10;

    public final static double cos(double angle) {

	double resid;

	resid = Math.abs(Math.IEEEremainder(angle,360.0));
	if (resid == 0.0) {
	    return 1.0;
	} else if (resid == 90.0) {
	    return 0.0;
	} else if (resid == 180.0) {
	    return -1.0;
	} else if (resid == 270.0) {
	    return 0.0;
	}

	return Math.cos(angle*d2r);
    }

    public final static double sin(double angle) {

	double resid;

	resid = Math.IEEEremainder(angle-90.0,360.0);
	if (resid == 0.0) {
	    return 1.0;
	} else if (resid == 90.0) {
	    return 0.0;
	} else if (resid == 180.0) {
	    return -1.0;
	} else if (resid == 270.0) {
	    return 0.0;
	}

	return Math.sin(angle*d2r);
    }

    public final static double tan(double angle) {

	double resid;

	resid = Math.IEEEremainder(angle,360.0);
	if (resid == 0.0 || Math.abs(resid) == 180.0) {
	    return 0.0;
	} else if (resid == 45.0 || resid == 225.0) {
	    return 1.0;
	} else if (resid == -135.0 || resid == -315.0) {
	    return -1.0;
	}

	return Math.tan(angle*d2r);
    }

    public final static double acos(double v) {

	if (v >= 1.0) {
	    if (v-1.0 <  tol) return 0.0;
	} else if (v == 0.0) {
	    return 90.0;
	} else if (v <= -1.0) {
	    if (v+1.0 > -tol) return 180.0;
	}

	return Math.acos(v)*r2d;
    }

    public final static double asin(double v) {

	if (v <= -1.0) {
	    if (v+1.0 > -tol) return -90.0;
	} else if (v == 0.0) {
	    return 0.0;
	} else if (v >= 1.0) {
	    if (v-1.0 <  tol) return 90.0;
	}

	return Math.asin(v)*r2d;
    }

    public final static double atan(double v) {

	if (v == -1.0) {
	    return -45.0;
	} else if (v == 0.0) {
	    return 0.0;
	} else if (v == 1.0) {
	    return 45.0;
	}

	return Math.atan(v)*r2d;
    }

    public final static double atan2(double y, double x) {

	if (y == 0.0) {
	    if (x >= 0.0) {
		return 0.0;
	    } else if (x < 0.0) {
		return 180.0;
	    }
	} else if (x == 0.0) {
	    if (y > 0.0) {
		return 90.0;
	    } else if (y < 0.0) {
		return -90.0;
	    }
	}

	return Math.atan2(y,x)*r2d;
    }
}

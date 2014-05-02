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

import FITSWCS.exceptions.*;

/**
 *   This class provides support for spherical map projections
 *   used by the FITS "World Coordinate System" (WCS) convention. <p>
 *
 *   The FITSWCS package was translated from the WCSLIB C library
 *   (V2.3).  This original library was written in support for coordinate 
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
 *   argument lists should make it clear what is intended. <p>
 *   
 *   <b> Accuracy </b><p>
 *   
 *   Closure to a precision of at least 1.0-10 degree of longitude and latitude
 *   has been verified for typical projection parameters on the 1 degree grid
 *   of native longitude and latitude (to within 5 degrees of any latitude
 *   where the projection may diverge). <p>
 *
 *   Notwithstanding this, absolutely no claim is made for the accuracy or
 *   reliability of these routines.  They are supplied as is, with no warranty
 *   of fitness for any purpose.
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
public abstract class Projection implements ProjectionType {

    /**
     * The radius of the generating sphere for the projection, a linear
     * scaling parameter.  If this is zero, it will be reset to the default
     * value of 180/pi (the value for FITS WCS).
     */
    protected double r0;

    /**
     * the ten projection parameters corresponding to the PROJPn keywords 
     * in FITS, so p[0] is PROJP0, and p[9] is PROJP9.  Many projections 
     * use p[1] (PROJP1) and some also use p[2] (PROJP2).  ZPN is the 
     * only projection which uses any of the others.
     */
    protected double[] p=null;

    /**
     * intermediate values derived from the projection parameters
     */
    protected double[] w;

    /**
     * an intermedate datum
     */
    protected int n;

    public final static double PI = Math.PI;
    public final static double D2R = PI / 180.0;
    public final static double R2D = 180.0 / PI;
    public final static double SQRT2 = Math.sqrt(2);
    public final static double SQRT2INV = 1.0 / SQRT2;

    /**
     * Compute (x,y) coordinates in the plane of projection from native 
     * spherical coordinates (phi,theta). <p>
     *
     * The values of phi and theta (the native longitude and latitude)
     * normally lie in the range [-180,180] for phi, and [-90,90] for theta.
     * However, all forward projections will accept any value of phi and will
     * not normalize it. <p>
     * 
     * Although many of the forward projections will accept values of theta
     * outside the range [-90,90] such latitudes are not meaningful and should
     * normally be marked as an error.  However, in the interests of
     * efficiency, the forward projection routines do not check for this,
     * although they do check for any invalid values of theta within the
     * range [-90,90].
     *
     * @return double[] a two-element array containing x,y
     */
    public abstract double[] fwd(double phi, double theta)
	throws PixelBeyondProjectionException;

    /**
     * same as fwd(phitheta[0], phitheta[1])
     */
    public double[] fwd(double[] phitheta)
	throws PixelBeyondProjectionException 
    { 
	return fwd(phitheta[0], phitheta[1]);
    }

    /**
     * Compute native spherical coordinates (phi,theta) from the 
     * (x,y) coordinates in the plane of projection. <p>
     *
     * Error checking on the projected coordinates (x,y) is limited to that
     * required to ascertain whether a solution exists.  Where a solution does
     * exist no check is made that the value of phi and theta obtained lie
     * within the ranges [-180,180] for phi, and [-90,90] for theta.
     *
     * @return double[] a two-element array containing phi,theta
     */
    public abstract double[] rev(double x, double y)
	throws PixelBeyondProjectionException;

    /**
     * same as rev(xy[0], xy[1])
     */
    public double[] rev(double[] xy)
	throws PixelBeyondProjectionException { return rev(xy[0], xy[1]); }

    /**
     * return the value of r0
     */
    public double getR0() { return r0; }

    /**
     * set the value of r0
     */
    public abstract void setR0(double r0);

    /**
     * set with new projection parameters.
     */
    public abstract void setProjParm(double[] p)
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException;

    /** 
     * return a copy of the projection parameters
     */
    public double[] getProjParm() { 
	double[] out = new double[p.length];
	System.arraycopy(p, 0, out, 0, p.length);
	return out;
    }

    public final static Projection getProjection(String pcode, 
						 double[] projparm) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException,
	       UnsupportedProjectionException
    {
	int i;
	Projection out;
	String pkg = "FITSWCS.projections.";

	for(i=1; 
	    i < NTYPES && 
		! pcode.equalsIgnoreCase(code[i]); 
	    i++);
	if (i >= NTYPES) throw new UnsupportedProjectionException(pcode);

	String classname = new String(pkg + code[i] + "Projection");
	try {
	    out = (Projection) Class.forName(classname).newInstance();
	}
	catch (ClassNotFoundException ex) {
	    throw new UnsupportedProjectionException(pcode + 
						     "(class not found)");
	}
	catch (IllegalAccessException ex) {
	    throw new UnsupportedProjectionException(pcode + 
						     "(unable to instantiate)");
	}
	catch (InstantiationException ex) {
	    throw new UnsupportedProjectionException(pcode + 
						     "(unable to instantiate)");
	}

	try {
	    out.setProjParm(projparm);
	} catch (ArrayIndexOutOfBoundsException ex) {
	    throw new ArrayIndexOutOfBoundsException(
		code[i] + ex.getMessage());
	} catch (BadProjectionParameterException ex) {
	    throw new BadProjectionParameterException(
		code[i] + ex.getMessage());
	}

	return out;
    }

    /**
     * return an integer representing the projection type which can be 
     * compared with this.ARC, SIN, TAN, etc.  -1 is returned if the 
     * name is not recognized.
     */
    public final static int projectionType(String name) {
	int i;
	for(i=0; i < NTYPES && 
		! name.equalsIgnoreCase(code[i]); 
	    i++);
	if (i == NTYPES) i = -1;
	return i;
    }

}



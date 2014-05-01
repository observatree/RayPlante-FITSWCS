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
 *   This class provides support for spherical coordinate
 *   transformations used by the FITS "World Coordinate System" (WCS) 
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
 *   argument lists should make it clear what is intended.
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
public class SphericalTransform {

    static final double tol = 1.0e-5;

    protected double[] euler;

    /**
     * Create a SphericalTransform object (storing euler angles internally)
     *    assuming the reference point is at the native pole and 
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     */
    public SphericalTransform(double reflong, double reflat)
	throws BadReferenceParameterException
    {
	euler = getEuler(reflong, reflat, 999.0, true, 90.0);
    }

    /**
     * Create a SphericalTransform object (storing euler angles internally)
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param pcode    the projection code of applied to the celestial system
     */
    public SphericalTransform(double reflong, double reflat, String pcode) 
	throws BadReferenceParameterException, UnsupportedProjectionException
    {
	double theta0 = getTheta0(pcode);
	euler = getEuler(reflong, reflat, 999.0, true, theta0);
    }

    /**
     * Create a SphericalTransform object (storing euler angles internally)
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param theta0   native latitude of reference point
     */
    public SphericalTransform(double reflong, double reflat, double theta0) 
	throws BadReferenceParameterException
    {
	euler = getEuler(reflong, reflat, 999.0, true, theta0);
    }

    /**
     * Create a SphericalTransform object (storing euler angles internally)
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param poleref  either the native longitude or native latitude of 
     *                 the celestial coordinate system, in degrees; a 
     *                 value of 999.0 indicates that this should be 
     *                 set to a default value.
     * @param islongpole true if poleref refers to the native longitude;
     *                 otherwise, its the native latitude.
     * @param theta0   native latitude of reference point
     */
    public SphericalTransform(double reflong, double reflat,
			      double poleref, boolean islongpole,
			      double theta0)
	throws BadReferenceParameterException
    {
	euler = getEuler(reflong, reflat, poleref, islongpole, theta0);
    }

    /**
     * Create a SphericalTransform object (storing euler angles internally)
     * @param ref[4] an array containing the reference point parameters:
     *        <dl> <dt> ref[0,1] <dd> celestial longitude & latitude 
     *                                (usually right ascension) of the 
     *                                reference point of the projection, 
     *                                in degrees
     *        <dl> <dt> ref[2,3] <dd> native longitude & latitude of the
     *                 pole of the celestial coordinate system, in degrees; 
     *                 either or both may have a value of 999.0, indicating 
     *                 that a default value should be calulated.
     *        </dl>
     * @param pcode    the projection code of applied to the celestial system
     */
    public SphericalTransform(double[] ref, String pcode) 
	throws BadReferenceParameterException, UnsupportedProjectionException
    {
	double poleref;
	boolean islongpole = true;
	double theta0 = getTheta0(pcode);

	if (ref[2] == 999.0) {
	    islongpole = false;
	    poleref = ref[3]; 
	}
	else {
	    poleref = ref[2];
	}
	
	euler = getEuler(ref[0], ref[1], poleref, islongpole, theta0);
    }

    /**
     * Create a SphericalTransform object, using the given (5) euler 
     * angles
     */
    public SphericalTransform(double[] euler) 
	throws ArrayIndexOutOfBoundsException
    {
	if (euler.length < 5) 
	    throw new ArrayIndexOutOfBoundsException(
		"input Euler array must have at least 5 elements");
	this.euler = new double[euler.length];
	System.arraycopy(euler, 0, this.euler, 0, euler.length);
    }

    /**
     * return a copy of the euler angles in use by this object
     */
    public double[] getEuler() { 
	double[] out = new double[euler.length];
	System.arraycopy(euler, 0, out, 0, euler.length);
	return out;
    }

    /**
     * Do a forward transformation, returning the result, phi and theta, 
     * as a double
     * @param lng spherical longitude in degrees
     * @param lat spherical latitude in degrees
     * @return phi and theta, the longitude and latitude in the native 
     *            coordinate system of the projection in degrees, returned 
     *            as a two element array
     */
    public double[] fwd(double lng, double lat) {
	return SphericalTransform.fwd(lng, lat, euler);
    }

    /**
     * same as fwd(lnglat[0], lnglat[1])
     */
    public double[] fwd(double[] lnglat) {
	return SphericalTransform.fwd(lnglat[0], lnglat[1], euler);
    }

    /**
     * Do a reverse transformation
     * @param phi Longitude in the native coordinate system of the 
     *            projection, in degrees.
     * @param theta Latitude in the native coordinate system of the 
     *            projection, in degrees.
     * @return spherical longitude and latitude in degrees, returned as a 
     *            two element array
     */
    public double[] rev(double phi, double theta) {
	return SphericalTransform.rev(phi, theta, euler);
    }

    /**
     * same as rev(phitheta[0], phitheta[1])
     */
    public double[] rev(double[] phitheta) {
	return SphericalTransform.rev(phitheta[0], phitheta[1], euler);
    }

    /**
     * Do a forward transformation, returning the result, phi and theta, 
     * as a double
     * @param lng spherical longitude in degrees
     * @param lat spherical latitude in degrees
     * @param eul[5] Euler angles for the transformation:
     *               <pre>
     *                  0: Celestial longitude of the native pole, in
     *                     degrees.
     *                  1: Celestial colatitude of the native pole, or
     *                     native colatitude of the celestial pole, in
     *                     degrees.
     *                  2: Native longitude of the celestial pole, in
     *                     degrees.
     *                  3: cos(eul[1])
     *                  4: sin(eul[1])
     *               </pre>
     * @return phi and theta, the longitude and latitude in the native 
     *            coordinate system of the projection in degrees, returned 
     *            as a two element array
     * @exception ArrayIndexOutOfBoundsException if eul has less than 5
     *               elements
     */
    public static double[] fwd(double lng, double lat, double[] eul) {
	double coslat, coslng, dlng, dphi, sinlat, sinlng, x, y, z;
	double phi, theta;
			   
	coslat = TrigD.cos(lat);
	sinlat = TrigD.sin(lat);

	dlng = lng - eul[0];
	coslng = TrigD.cos(dlng);
	sinlng = TrigD.sin(dlng);

	// Compute the native longitude.
	x = sinlat*eul[4] - coslat*eul[3]*coslng;
	if (Math.abs(x) < tol) {

	    // Rearange formula to reduce roundoff errors.
	    x = -TrigD.cos(lat+eul[1]) + coslat*eul[3]*(1.0 - coslng);
	}
	y = -coslat*sinlng;
	if (x != 0.0 || y != 0.0) {
	    dphi = TrigD.atan2(y,x);
	} else {

	    // Change of origin of longitude.
	    dphi = dlng - 180.0;
	}
	phi = eul[2] + dphi;

	// Normalize the native longitude.
	if (phi > 180.0) {
	    phi -= 360.0;
	} else if (phi < -180.0) {
	    phi += 360.0;
	}

	// Compute the native latitude. 
	if (Math.IEEEremainder(dlng,180.0) == 0.0) {
	    theta = lat + coslng*eul[1];
	    if (theta >  90.0) theta =  180.0 - theta;
	    if (theta < -90.0) theta = -180.0 - theta;
	} else {
	    z = sinlat*eul[3] + coslat*eul[4]*coslng;
	    if (Math.abs(z) > 0.99) {

		// Use an alternative formula for greater numerical accuracy. 
		double tmp = TrigD.acos(Math.sqrt(x*x+y*y));
		theta = (z < 0.0) ? -Math.abs(tmp) : Math.abs(tmp);
	    } else {
		theta = TrigD.asin(z);
	    }
	}

	// return result as a two element array
	double[] out = { phi, theta };
	return out;
    }

    /**
     * Do a reverse transformation
     * @param phi Longitude in the native coordinate system of the 
     *            projection, in degrees.
     * @param theta Latitude in the native coordinate system of the 
     *            projection, in degrees.
     * @param eul[5] Euler angles for the transformation:
     *               <pre>
     *                  0: Celestial longitude of the native pole, in
     *                     degrees.
     *                  1: Celestial colatitude of the native pole, or
     *                     native colatitude of the celestial pole, in
     *                     degrees.
     *                  2: Native longitude of the celestial pole, in
     *                     degrees.
     *                  3: cos(eul[1])
     *                  4: sin(eul[1])
     *               </pre>
     * @return spherical longitude and latitude in degrees, returned as a 
     *            two element array
     * @exception ArrayIndexOutOfBoundsException if eul has less than 5
     *               elements
     */
    public static double[] rev(double phi, double theta, double[] eul) {

	double cosphi, costhe, dlng, dphi, sinphi, sinthe, x, y, z;
	double lng, lat;

	costhe = TrigD.cos(theta);
	sinthe = TrigD.sin(theta);

	dphi = phi - eul[2];
	cosphi = TrigD.cos(dphi);
	sinphi = TrigD.sin(dphi);

	// Compute the celestial longitude. 
	x = sinthe*eul[4] - costhe*eul[3]*cosphi;
	if (Math.abs(x) < tol) {

	    // Rearrange formula to reduce roundoff errors. 
	    x = - TrigD.cos(theta+eul[1]) + 
		costhe*eul[3]*(1.0 - cosphi);
	}
	y = -costhe*sinphi;
	if (x != 0.0 || y != 0.0) {
	    dlng = TrigD.atan2(y, x);
	} else {

	    // Change of origin of longitude. 
	    dlng = dphi + 180.0;
	}
	lng = eul[0] + dlng;

	// Normalize the celestial longitude. 
	if (eul[0] >= 0.0) {
	    if (lng < 0.0) lng += 360.0;
	} else {
	    if (lng > 0.0) lng -= 360.0;
	}

	if (lng > 360.0) {
	    lng -= 360.0;
	} else if (lng < -360.0) {
	    lng += 360.0;
	}

	// Compute the celestial latitude. 
	if (Math.IEEEremainder(dphi,180.0) == 0.0) {
	    lat = theta + cosphi*eul[1];
	    if (lat >  90.0) lat =  180.0 - lat;
	    if (lat < -90.0) lat = -180.0 - lat;
	} else {
	    z = sinthe*eul[3] + costhe*eul[4]*cosphi;
	    if (Math.abs(z) > 0.99) {

		// Use an alternative formula for greater numerical accuracy. 
		double tmp = TrigD.acos(Math.sqrt(x*x+y*y));
		lat = (z < 0.0) ? -Math.abs(tmp) : Math.abs(tmp);
	    } else {
		lat = TrigD.asin(z);
	    }
	}

	// return result as a two element array
	double[] out = { lng, lat };
	return out;
    }

    /**
     * Compute the euler angles for a given set of reference angles
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param poleref  either the native longitude or native latitude of 
     *                 the celestial coordinate system, in degrees; a 
     *                 value of 999.0 indicates that this should be 
     *                 set to a default value.
     * @param islongpole true if poleref refers to the native longitude;
     *                 otherwise, its the native latitude.
     * @param theta0   native latitude of reference point
     */
    public final static double[] getEuler(double reflong, double reflat,
					  double poleref, boolean islongpole,
					  double theta0)
	throws BadReferenceParameterException
    {
	boolean dophip;
	double tol = 1.0e-10;
	double clat0, cphip, cthe0, slat0, sphip, sthe0;
	double latp, latp1, latp2;
	double u, v, x, y, z;
	double tmp;
	int i;
	double[] euler = new double[5];

	double[] ref = { reflong, reflat, poleref, 999.0 };
	if (! islongpole) { 
	    tmp = ref[2];
	    ref[2] = ref[3];
	    ref[3] = tmp;
	}

	// Set default for native longitude of the celestial pole? 
	dophip = (ref[2] == 999.0);

	// Compute celestial coordinates of the native pole. 
	if (theta0 == 90.0) {

	    // Reference point is at the native pole. 
	    if (dophip) {

		// Set default for longitude of the celestial pole. 
		ref[2] = 180.0;
	    }

	    latp = ref[1];
	    ref[3] = latp;

	    euler[0] = ref[0];
	    euler[1] = 90.0 - latp;

	} else {

	    // Reference point away from the native pole. 
	    // Set default for longitude of the celestial pole. 
	    if (dophip) {
		ref[2] = (ref[1] < theta0) ? 180.0 : 0.0;
	    }

	    clat0 = TrigD.cos(ref[1]);
	    slat0 = TrigD.sin(ref[1]);
	    cphip = TrigD.cos(ref[2]);
	    sphip = TrigD.sin(ref[2]);
	    cthe0 = TrigD.cos(theta0);
	    sthe0 = TrigD.sin(theta0);

	    x = cthe0*cphip;
	    y = sthe0;
	    z = Math.sqrt(x*x + y*y);
	    if (z == 0.0) {
		if (slat0 != 0.0) 
		    throw new BadReferenceParameterException(
			"longpole=" + ref[2] + 
			" is incompatible with reference latitude=" + ref[1]);

		// latp determined by LATPOLE in this case. 
		latp = ref[3];
	    } else {
		if (Math.abs(slat0/z) > 1.0) 
		    throw new BadReferenceParameterException(
			"reference latitude=" + ref[1] + 
			" is incompatible with longpole=" + ref[2]);

		u = TrigD.atan2(y,x);
		v = TrigD.cos(slat0/z);

		latp1 = u + v;
		if (latp1 > 180.0) {
		    latp1 -= 360.0;
		} else if (latp1 < -180.0) {
		    latp1 += 360.0;
		}

		latp2 = u - v;
		if (latp2 > 180.0) {
		    latp2 -= 360.0;
		} else if (latp2 < -180.0) {
		    latp2 += 360.0;
		}

		if (Math.abs(ref[3]-latp1) < Math.abs(ref[3]-latp2)) {
		    if (Math.abs(latp1) < 90.0+tol) {
			latp = latp1;
		    } else {
			latp = latp2;
		    }
		} else {
		    if (Math.abs(latp2) < 90.0+tol) {
			latp = latp2;
		    } else {
			latp = latp1;
		    }
		}

		ref[3] = latp;
	    }

	    euler[1] = 90.0 - latp;

	    z = TrigD.cos(latp)*clat0;
	    if (Math.abs(z) < tol) {
		if (Math.abs(clat0) < tol) {

		    // Celestial pole at the reference point. 
		    euler[0] = ref[0];
		    euler[1] = 90.0 - theta0;
		} else if (latp > 0.0) {

		    // Celestial pole at the native north pole.
		    euler[0] = ref[0] + ref[2] - 180.0;
		    euler[1] = 0.0;
		} else if (latp < 0.0) {

		    // Celestial pole at the native south pole. 
		    euler[0] = ref[0] - ref[2];
		    euler[1] = 180.0;
		}
	    } else {
		x = (sthe0 - TrigD.sin(latp)*slat0)/z;
		y =  sphip*cthe0/clat0;
		if (x == 0.0 && y == 0.0) 
		    throw new BadReferenceParameterException(
			"Unable to calculate euler parameters for given " +
			"reference parameters");

		euler[0] = ref[0] - TrigD.atan2(y,x);
	    }

	    // Make euler[0] the same sign as ref[0]. 
	    if (ref[0] >= 0.0) {
		if (euler[0] < 0.0) euler[0] += 360.0;
	    } else {
		if (euler[0] > 0.0) euler[0] -= 360.0;
	    }
	}

	euler[2] = ref[2];
	euler[3] = TrigD.cos(euler[1]);
	euler[4] = TrigD.sin(euler[1]);

	// Check for ill-conditioned parameters. 
	if (Math.abs(latp) > 90.0+tol) 
	    throw new BadReferenceParameterException(
		"Ill-conditioned reference paramters");

	return euler;
    }

    /**
     * Compute the euler angles for a given set of reference angles, 
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param theta0   native latitude of reference point
     */
    public final static double[] getEuler(double reflong, double reflat, 
					  double theta0) 
	throws BadReferenceParameterException
    {
	return getEuler(reflong, reflat, 999.0, true, theta0);
    }

    /**
     * Compute the euler angles for a given set of reference angles, 
     *    assuming the reference point is at the native pole and 
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     */
    public final static double[] getEuler(double reflong, double reflat)
	throws BadReferenceParameterException
    {
	return getEuler(reflong, reflat, 999.0, true, 90.0);
    }

    /**
     * Compute the euler angles for a given set of reference angles, 
     *    taking defaults for the longitude & latitude of the pole
     *    of the reference system.
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param pcode    the projection code of applied to the celestial system
     */
    public final static double[] getEuler(double reflong, double reflat, 
					  String pcode) 
	throws BadReferenceParameterException, UnsupportedProjectionException
    {
	double theta0 = getTheta0(pcode);
	return getEuler(reflong, reflat, 999.0, true, theta0);
    }

    /**
     * Compute the euler angles for a given set of reference angles
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param poleref  either the native longitude or native latitude of the
     *                 pole of the celestial coordinate system, in degrees; 
     *                 a value of 999.0 indicates that this should be 
     *                 set to a default value.
     * @param islongpole true if poleref refers to the native longitude;
     *                 otherwise, its the native latitude.
     * @param pcode    the projection code of applied to the celestial system
     */
    public final static double[] getEuler(double reflong, double reflat,
					  double poleref, boolean islongpole,
					  String pcode)
	throws BadReferenceParameterException, UnsupportedProjectionException
    {
	double theta0 = getTheta0(pcode);
	return getEuler(reflong, reflat, poleref, islongpole, theta0);
    }

    /**
     * determine the proper value of theta0 (used by getEuler() and the 
     * SphericalTransform constructor) for a given projection code.
     */
    public static double getTheta0(String pcode) 
	throws UnsupportedProjectionException
    {
	int i;

	for(i=0;
	    i < ProjectionType.NTYPES && 
	      ! pcode.equalsIgnoreCase(ProjectionType.code[i]);
	    i++);
	if (i >= ProjectionType.NTYPES) 
		throw new UnsupportedProjectionException(
		    "Unrecognized Projection code: " + pcode);

	return ((i <= ProjectionType.AIR) ? 90.0 : 0.0);
    }
}

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
import FITSWCS.projections.*;

/**
 *   This class provides support for celestial coordinate transformations 
 *   (encapsulating spherical coordinate and projection transformations)  
 *   used by the FITS "World Coordinate System" (WCS) convention. <p>
 *
 *   The FITSWCS package was translated from the WCSLIB C library
 *   This original library was written in support for coordinate 
 *   systems used by astronomical data stored in FITS format.  For more 
 *   information on these coordinate systems, refer to the paper by Greisen 
 *   and Calabretta at:
 *   <blockquote>
 *       ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z 
 *   </blockquote>
 *
 *   <a name="refdat"><b> Coordinate System Reference Parameters </b></a><p>
 *
 *   This class can be instantiated by passing to the constructor an 
 *   array of four double values representing the coordinate system reference 
 *   data.  The first two elements should be set to the celestial longitude
 *   and latitude (usually right ascension and declination) of the
 *   reference point of the projection. <p>
 *   
 *   The two values are the native longitude and latitude of
 *   the pole of the celestial coordinate system and correspond to the
 *   FITS keywords LONGPOLE and LATPOLE. <p>
 *   
 *   LONGPOLE defaults to 0 degrees if the celestial latitude of the
 *   reference point of the projection is greater than the native
 *   latitude, otherwise 180 degrees.  (This is the condition for the
 *   celestial latitude to increase in the same direction as the native
 *   latitude at the reference point.)  ref[2] may be set to 999.0 to
 *   indicate that the correct default should be substituted. <p>
 *   
 *   In some circumstances the latitude of the native pole may be
 *   determined by the first three values only to within a sign and
 *   LATPOLE is used to choose between the two solutions.  LATPOLE is
 *   set in ref[3] and the solution closest to this value is used to
 *   reset ref[3].  It is therefore legitimate, for example, to set
 *   ref[3] to 999.0 to choose the more northerly solution - the default
 *   if the LATPOLE card is omitted from the FITS header.  For the
 *   special case where the reference point of the projection is at
 *   native latitude zero, its celestial latitude is zero, and
 *   LONGPOLE = +/- 90 then the native latitude of the pole is not
 *   determined by the first three reference values and LATPOLE
 *   specifies it completely. <p>
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
public class CelestialTransform {

    /**
     * The projection type code (used only for error messages)
     */
    protected String pcode;

    /**
     * a double array containing the coordinate system reference parameters.
     * See the above section on <a href="#refdat">reference parameters</a> for
     * more details.
     */
    protected double[] ref;

    /**
     * Euler angles and associated intermediaries derived from the
     * coordinate reference values.
     */
    protected double[] euler;

    /**
     * The Projection transform object to use
     */
    protected Projection prj;

    /**
     * Construct a CelestialTransform object
     * @param pcode  3-character code indicating desired projection
     * @param ref    4-element array containing coordinate reference 
     *               parameters (see <a href="#refdat">above</a> for 
     *               more details).
     * @param r0     sphere radius, if 0 defaults to 180/PI.
     * @param p      array containing projection parameters (up to 
     *               10 elements may be used, depending on pcode).
     * @exception ArrayIndexOutOfBoundsException if ref.length < 4 or if 
     *               p.length < number needed for particular pcode.
     * @exception BadProjectionParameterException if p contains one or more
     *               bad values for the given pcode
     * @exception UnsupportedProjectionException if pcode is refers to an 
     *               unrecognized projection type.
     * @exception BadReferenceParameterException if ref contains one or more
     *               bad values for the given pcode
     */
    public CelestialTransform(String pcode, double[] ref, double r0, double[] p)
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException,
	       UnsupportedProjectionException, BadReferenceParameterException
    {
	init(pcode, ref, r0, p);
    }

    /**
     * Construct a CelestialTransform object
     * @param pcode  3-character code indicating desired projection
     * @param ref    4-element array containing coordinate reference 
     *               parameters (see <a href="#refdat">above</a> for 
     *               more details).
     * @param p      array containing projection parameters (up to 
     *               10 elements may be used, depending on pcode).
     * @exception ArrayIndexOutOfBoundsException if ref.length < 4 or if 
     *               p.length < number needed for particular pcode.
     * @exception BadProjectionParameterException if p contains one or more
     *               bad values for the given pcode
     * @exception UnsupportedProjectionException if pcode is refers to an 
     *               unrecognized projection type.
     * @exception BadReferenceParameterException if ref contains one or more
     *               bad values for the given pcode
     */
    public CelestialTransform(String pcode, double[] ref, double[] p) 
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException,
	       UnsupportedProjectionException, BadReferenceParameterException
    {
	init(pcode, ref, 0, p);
    }

    /**
     * Construct a CelestialTransform object taking the reference parameters
     * as four seperate arguments (reflong, reflat, longpole, latpole).
     * See the above section on <a href="#refdat">reference parameters</a> for
     * more details.
     * @param pcode    3-character code indicating desired projection
     * @param reflong  celestial longitude (usually right ascension) of 
     *                 the reference point of the projection, in degrees
     * @param reflat   celestial latitude (usually declination) of 
     *                 the reference point of the projection, in degrees
     * @param longpole native longitude of the celestial coordinate system,
     *                 in degrees; a value of 999.0 indicates that this 
     *                 should be set to a default value (either 0 or 180, 
     *                 so that celestial and native latitudes increase in 
     *                 the same direction).
     * @param latpole  native latiitude of the celestial coordinate system,
     *                 in degrees; a value of 999.0 indicates that this 
     *                 should be set to a default value (derived from 
     *                 reflong, reflat, & longpole).
     * @param p        array containing projection parameters (up to 
     *                 10 elements may be used, depending on pcode).
     * @exception ArrayIndexOutOfBoundsException if p.length < number 
     *                 needed for particular pcode.
     * @exception BadProjectionParameterException if p contains one or more
     *                 bad values for the given pcode
     * @exception UnsupportedProjectionException if pcode is refers to an 
     *                 unrecognized projection type.
     * @exception BadReferenceParameterException if any of reflong, reflat, 
     *                 longpole, or latpole contains a bad value for the 
     *                 given pcode
     */
    public CelestialTransform(String pcode, double reflong, double reflat,
			      double longpole, double latpole, double[] p)
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException,
	       UnsupportedProjectionException, BadReferenceParameterException
    {
	double[] use = { reflong, reflat, longpole, latpole };
	init(pcode, use, 0, p);
    }

    private void init(String pcode, double[] refdat, double r0, double[] p)
	throws ArrayIndexOutOfBoundsException, BadProjectionParameterException,
	       UnsupportedProjectionException, BadReferenceParameterException
    {
	boolean dophip;
	double tol = 1.0e-10;
	double clat0, cphip, cthe0, theta0, slat0, sphip, sthe0;
	double latp, latp1, latp2;
	double u, v, x, y, z;
	int i;

	// remember the pcode
	this.pcode = pcode;

	// make a copy of coordinate system reference parameters
	if (refdat.length < 4) throw new 
	    ArrayIndexOutOfBoundsException(
		"Need at least 4 elements in reference data array");

	ref = new double[refdat.length];
	System.arraycopy(refdat, 0, ref, 0, ref.length);

	if (p == null) { 
	    p = new double[1]; 
	    p[0] = 0.0;
	}

	// Determine the proper projection to apply
	try {
	    if (pcode.equalsIgnoreCase("AZP")) {
		prj = new AZPProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("TAN")) {
		prj = new TANProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("SIN")) {
		prj = new SINProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("STG")) {
		prj = new STGProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("ARC")) {
		prj = new ARCProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("ZPN")) {
		prj = new ZPNProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("ZEA")) {
		prj = new ZEAProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("AIR")) {
		prj = new AIRProjection(r0, p);
		theta0 = 90.0;
	    } else if (pcode.equalsIgnoreCase("CYP")) {
		prj = new CYPProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("CAR")) {
		prj = new CARProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("MER")) {
		prj = new MERProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("CEA")) {
		prj = new CEAProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("COP")) {
		prj = new COPProjection(r0, p);
		theta0 = p[1];
	    } else if (pcode.equalsIgnoreCase("COD")) {
		prj = new CODProjection(r0, p);
		theta0 = p[1];
	    } else if (pcode.equalsIgnoreCase("COE")) {
		prj = new COEProjection(r0, p);
		theta0 = p[1];
	    } else if (pcode.equalsIgnoreCase("COO")) {
		prj = new COOProjection(r0, p);
		theta0 = p[1];
	    } else if (pcode.equalsIgnoreCase("BON")) {
		prj = new BONProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("PCO")) {
		prj = new PCOProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("GLS")) {
		prj = new GLSProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("PAR")) {
		prj = new PARProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("AIT")) {
		prj = new AITProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("MOL")) {
		prj = new MOLProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("CSC")) {
		prj = new CSCProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("QSC")) {
		prj = new QSCProjection(r0, p);
		theta0 = 0.0;
	    } else if (pcode.equalsIgnoreCase("TSC")) {
		prj = new TSCProjection(r0, p);
		theta0 = 0.0;
	    } else {

		// Unrecognized projection code. 
		throw new UnsupportedProjectionException(
		    "Unrecognized Projection code: " + pcode);
	    }
	}
	catch(BadProjectionParameterException ex) {
	    throw new BadProjectionParameterException(pcode + ": "  +
						      ex.getMessage());
	}
	catch(ArrayIndexOutOfBoundsException ex) {
	    throw new ArrayIndexOutOfBoundsException(pcode  + ": "  +
						      ex.getMessage());
	}

	// Set default for native longitude of the celestial pole? 
	dophip = (ref[2] == 999.0);

	euler = new double[5];

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
			"Bad celestial position for this projection: " +
			ref[0] + ", " + ref[1] + " " + pcode);

		// latp determined by LATPOLE in this case. 
		latp = ref[3];
	    } else {
		if (Math.abs(slat0/z) > 1.0) 
		    throw new BadReferenceParameterException(
			"Bad celestial position for this projection: " +
			ref[0] + ", " + ref[1] + " " + pcode);

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
    }

    /**
     * Compute (x,y) coordinates in the plane of projection from celestial
     * coordinates (lng,lat).
     * @exception InvalidCelestialTransformException if the given coordinates
     *              are invalid for this system
     */
    public double[] fwd(double lng, double lat) 
	throws InvalidCelestialTransformException
    {
	int    err;
	double[] out, phitheta;

	// Compute native coordinates. 
	phitheta = SphericalTransform.fwd(lng, lat, euler);

	// Apply forward projection. 
	try {
	    out = prj.fwd(phitheta[0], phitheta[1]);
	}
	catch (PixelBeyondProjectionException ex) {
	    throw new InvalidCelestialCoordException(pcode, lng, lat);
	}
	return out;
    }

    public double[] rev(double x, double y) 
	throws InvalidCelestialTransformException
    {
	int    err;
	double[] out, phitheta;

	// Apply reverse projection. 
	try {
	    phitheta = prj.rev(x, y);
	}
	catch (PixelBeyondProjectionException ex) {
	    throw new InvalidMapCoordException(pcode, x, y);
	}

	// Compute native coordinates. 
	out = SphericalTransform.rev(phitheta[0], phitheta[1], euler);
	return out;
    }

    /**
     * return a copy of the projection code
     */
    public String getProjectionCode() { return new String(pcode); }

    /** 
     * return the projection object currently in use
     */
    public Projection getProjection() { return prj; }

    /**
     * return a copy of the reference system parameters 
     */
    public double[] getRefParm() {
	double[] out = new double[ref.length];
	System.arraycopy(ref, 0, out, 0, ref.length);
	return out;
    }

    /**
     * return a copy of the Euler angles currently in use
     */
    public double[] getEuler() {
	double[] out = new double[euler.length];
	System.arraycopy(euler, 0, out, 0, euler.length);
	return out;
    }
}

	

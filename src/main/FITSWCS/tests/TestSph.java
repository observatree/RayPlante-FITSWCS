/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS.tests;

import FITSWCS.*;
import FITSWCS.exceptions.*;

/**
 *   This class verifies the SphericalTransform class for closure errors.<p>
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
public class TestSph {

    public static void main(String args[]) {

	int err, j, lat, lng;
	double[] eul = new double[5];
//	double coslat, lng1, lng2, lat1, lat2, phi, theta, zeta;
	double coslat=-90;
	double lat1, lat2;
	double phi, theta, zeta;
	double lng1, lng2;
	double DEGTORAD = Math.PI / 180.0;
	double tol = 1.0e-12;
	double[] out;
 
	System.out.println("Testing closure of WCSLIB coordinate " + 
			   "transformation routines");
	System.out.println("-------------------------------------" + 
			   "-----------------------");

	// Set reference angles. 
	eul[0] =  90.0;
	eul[1] =  30.0;
	eul[2] = -90.0;

	System.out.println("\nCelestial longitude and latitude of the " +
			   "native pole, and native");
	System.out.print("longitude of the celestial pole (degrees): ");
        System.out.format("%-10.4f", eul[0])
                  .format("%-10.4f", eul[1])
                  .format("%-10.4f", eul[2]).println();

	eul[3] = Math.cos(eul[1] * DEGTORAD);
	eul[4] = Math.sin(eul[1] * DEGTORAD);

	System.out.format("Closure tolerance: %-7.1e degrees of arc.\n\n", tol);

	for (lat = 90; lat >= -90; lat--) {
	    lat1 = lat;
	    coslat = Math.cos(lat1 * DEGTORAD);

	    for (lng = -180; lng <= 180; lng++) {
		lng1 = lng;

		out = SphericalTransform.fwd(lng1, lat1, eul);
		phi = out[0]; 
		theta = out[1];
		out = SphericalTransform.rev(phi, theta, eul);
		lng2 = out[0];
		lat2 = out[1];

		if (Math.abs(lat2-lat1) > tol || 
		    (Math.abs(lng2-lng1)-360.0)*coslat > tol) {

                System.out.printf("Unclosed: lng1 = %-20.15f  lat1 = %-20.15f",
                                  lng1, lat1).println();
                System.out.printf("           phi = %-20.15f theta = %-20.15f",
                                  phi, theta).println();
                System.out.printf("          lng2 = %-20.15f  lat2 = %-20.15f",
                                  lng2, lat2).println();
		}
	    }
	}

	// Test closure at points close to the pole. 
	for (j = -1; j <= 1; j += 2) {
	    zeta = 1.0;
	    lng1 = -180.0;

	    for (lat = 0; lat < 12; lat++) {
		lat1 = j*(90.0 - zeta);

		out = SphericalTransform.fwd(lng1, lat1, eul);
		phi = out[0];
		theta = out[1];

		out = SphericalTransform.rev(phi, theta, eul); 
		lng2 = out[0];
		lat2 = out[1];

		if (Math.abs(lat2-lat1) > tol || 
		    (Math.abs(lng2-lng1)-360.0)*coslat > tol) {
		
                System.out.printf("Unclosed: lng1 = %-20.15f  lat1 = %-20.15f",
                                  lng1, lat1).println();
                System.out.printf("           phi = %-20.15f theta = %-20.15f",
                                  phi, theta).println();
                System.out.printf("          lng2 = %-20.15f  lat2 = %-20.15f",
                                  lng2, lat2).println();
		}

		zeta /= 10.0;
		lng1 += 30.0;
	    }
	}

	System.out.println("Testing complete.");
    }
}

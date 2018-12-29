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
import java.util.BitSet;

/**
 *   This class verifies the CelestialTransform class for closure errors.<p>
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
public class TestCel {

    public static void main(String args[]) {

	double tol = 1.0e-10;
	double[] p = new double[10];
	double[] r = new double[4];
	CelestialTransform cel;

	System.out.println("Testing closure of WCSLIB coordinate " + 
			   "transformation routines");
	System.out.println("-------------------------------------" + 
			   "-----------------------");

	BitSet doproj = whichProj(args);

	for (int j = 0; j < 10; p[j++] = 0.0);
	r[0] = 54.5;     // sample celestial reference longitude
	r[1] = 32.125;   // sample celestial reference latitude
	r[2] = 999.0;    // use default for LONGPOLE
	r[3] = 999.0;    // use default for LATPOLE

	if (doproj.get(AZP)) {
	    
	    // AZP: zenithal/azimuthal perspective.
	    p[1] = 2.0;
	    try {
		cel = new CelestialTransform("AZP", r, p);
		runTest(cel, 90, -30, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }

	} 
	if (doproj.get(TAN)) {

	    // TAN: gnomonic.  
            // Note: projection goes unstable far from reference position;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.0;
	    try {
		cel = new CelestialTransform("TAN", ref, p);
		runTest(cel, 90, 5, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(SIN)) {

	    // SIN: orthographic/synthesis
	    p[1] = 0.3;
	    p[2] = 1.5;
	    try {
		cel = new CelestialTransform("SIN", r, p);
		runTest(cel, 90, 60, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(STG)) {

	    // STG: stereographic.
            // Note: projection goes unstable far from reference position;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.0;
	    try {
		cel = new CelestialTransform("STG", ref, p);
		runTest(cel, 90, -85, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(ARC)) {

	    // ARC: zenithal/azimuthal equidistant.
	    try {
		cel = new CelestialTransform("ARC", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(ZPN)) {

	    // ZPN: zenithal/azimuthal equidistant.
	    p[0] =  0.00000;
	    p[1] =  0.95000;
	    p[2] = -0.02500;
	    p[3] = -0.15833;
	    p[4] =  0.00208;
	    p[5] =  0.00792;
	    p[6] = -0.00007;
	    p[7] = -0.00019;
	    p[8] =  0.00000;
	    p[9] =  0.00000;
	    try {
		cel = new CelestialTransform("ZPN", r, p);
		runTest(cel, 90, 10, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(ZEA)) {

	    // ZEA: zenithal/azimuthal equal area. 
	    try {
		cel = new CelestialTransform("ZEA", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(AIR)) {
	    
	    // AIR: zenithal/azimuthal perspective.
            // Note: projection goes unstable far from reference position;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.0;
	    p[1] = 45.0;
	    try {
		cel = new CelestialTransform("AIR", ref, p);
		runTest(cel, 90, -85, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(CYP)) {
	    
	    // CYP: cylindrical perspective.
	    p[1] = 3.0;
	    p[2] = 0.8;
	    try {
		cel = new CelestialTransform("CYP", r, p);
		runTest(cel, 90, -90, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(CAR)) {

	    // CAR: Cartesian. 
	    try {
		cel = new CelestialTransform("CAR", r, p);
		runTest(cel, 90, -90, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(MER)) {

	    // MER: Cartesian. 
	    try {
		cel = new CelestialTransform("MER", r, p);
		runTest(cel, 85, -85, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(CEA)) {
	    
	    // CEA: cylindrical equal area.
	    p[1] = 0.75;
	    try {
		cel = new CelestialTransform("CEA", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(COP)) {
	    
	    // COP: conic perspective.
            // Note: projection goes unstable far from reference position;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89;

	    p[1] = 60.0;
	    p[2] = 15.0;
	    try {
		cel = new CelestialTransform("COP", ref, p);
		runTest(cel, 90, 0, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println("Failed to test COP: "+ex.getMessage());
	    }
	}
	if (doproj.get(COD)) {
	    
	    // COD: conic equidistant
            // Note: projection goes unstable far from reference position
            // and opposite reference longitude;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.5;

	    p[1] = -60.0;
	    p[2] = 15.0;
	    try {
		cel = new CelestialTransform("COD", ref, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(COE)) {
	    
	    // COE: conic equal area.
            // Note: projection goes unstable near reference pole
            // and opposite reference longitude;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.5;

	    p[1] = 60.0;
	    p[2] = -15.0;
	    try {
		cel = new CelestialTransform("COE", ref, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println("Failed to test COE: "+ex.getMessage());
	    }
	}
	if (doproj.get(COO)) {
	    
	    // COO: conic orthomorphic.
            // Note: projection goes unstable far from reference position;
            // Thus, we choose a ref. position near +90.  
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.5;

	    p[1] = -60.0;
	    p[2] = -15.0;
	    try {
		cel = new CelestialTransform("COO", ref, p);
		runTest(cel, 85, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(BON)) {

	    // BON: Bonne's projection.
	    p[1] = 30.0;
	    try {
		cel = new CelestialTransform("BON", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(PCO)) {

	    // PCO: polyconic.
            // Note: projection goes unstable near pole and the 
            // reference longitude
            double[] ref = new double[4];
            System.arraycopy(r, 0, ref, 0, r.length);
            ref[1] = 89.0;

	    try {
		cel = new CelestialTransform("PCO", ref, p);
		runTest(cel, 88, -88, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(GLS)) {

	    // GLS: Sanson-Flamsteed (global sinusoid). 
	    try {
		cel = new CelestialTransform("GLS", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(PAR)) {

	    // PAR: parabolic. 
	    try {
		cel = new CelestialTransform("PAR", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(AIT)) {

	    // AIT: Hammer-Aitoff.
	    try {
		cel = new CelestialTransform("AIT", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(MOL)) {

	    // MOL: Mollweide's projection. 
	    try {
		cel = new CelestialTransform("MOL", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(CSC)) {

	    // CSC: COBE quadrilateralized spherical cube.
	    try {
		cel = new CelestialTransform("QSC", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(QSC)) {

	    // QSC: quadrilateralized spherical cube.
	    try {
		cel = new CelestialTransform("QSC", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(TSC)) {

	    // TSC: tangential spherical cube.
	    try {
		cel = new CelestialTransform("TSC", r, p);
		runTest(cel, 90, -89, tol);
	    }
	    catch (FITSWCSException ex) {
		System.out.println(ex.getMessage());
	    }
	}  
    }

    protected static void runTest(CelestialTransform cel,
				  int north, int south, double tol) 
    {
	int err, lat, lng, j;
	double dlat, dlatmx, dlng, dlngmx, dr, drmax, lat1, lat2, lng1, lng2;
	double r, theta, x, x1, x2, y, y1, y2;
	double[] out;
	String pcode = cel.getProjectionCode();
	Projection prj = cel.getProjection();
	double ux, uy;

	System.out.printf("Testing %s native latitudes %d to %d, ",
                          pcode, north, south)
                  .printf("closure tolerance %8.1e deg.\n", tol);

	drmax = 0.0;

	for (lat = north; lat >= south; lat--) {
	    lat1 = (double)lat;

	    for (lng = -180; lng <= 180; lng++) {
		lng1 = (double)lng;

		try {
		    out = cel.fwd(lng1, lat1);
		}
		catch (InvalidCelestialTransformException ex) {
		    System.out.printf("Error: lng1 = %20.15f  lat1 = %20.15f\n",
                                      lng1, lat1);
		    System.out.println("       " + ex.getMessage());
		    continue;
		}
		x1 = out[0];
		y1 = out[1];

		try {
		    out = cel.rev(x1, y1);
		}
		catch (InvalidCelestialTransformException ex) {
		    System.out.printf("Error: lng1 = %20.15f  lat1 = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("          x = %20.15f     y = %20.15f\n",
                                      x1, y1);
		    System.out.println("       " + ex.getMessage());
		    continue;
		}
		lng2 = out[0];
		lat2 = out[1];

		try {
		    out = cel.fwd(lng2, lat2);
		}
		catch (InvalidCelestialTransformException ex) {
		    System.out.printf("Error: lng1 = %20.15f  lat1 = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("         x1 = %20.15f    y1 = %20.15f\n",
                                      x1, y1);
		    System.out.printf("Error: lng2 = %20.15f  lat2 = %20.15f\n",
                                      lng2, lat2);
		    System.out.println("       " + ex.getMessage());
		    continue;
		}
		x2 = out[0];
		y2 = out[1];

		dr = Math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
		if (dr > drmax) drmax = dr;
		if (dr > tol) {
		    System.out.printf("  lng = %20.15f  lat = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("  %s:    x1 = %20.15f    y1 = %20.15f\n",
                                      pcode, x1, y1);
		    System.out.printf("     : lng2 = %20.15f  lat2 = %20.15f\n",
                                      lng2, lat2);
		    System.out.printf("     :   x2 = %20.15f    y2 = %20.15f\n",
                                      x2, y2);
		}
	    }
	}

	System.out.printf("  Maximum residual (map): dR: %12.6e\n", drmax);
    }

    protected static final int nprojs = 26;
    protected static final int NON = 0;
    protected static final int AZP = 1;
    protected static final int TAN = 2;
    protected static final int SIN = 3;
    protected static final int STG = 4;
    protected static final int ARC = 5;
    protected static final int ZPN = 6;
    protected static final int ZEA = 7;
    protected static final int AIR = 8;
    protected static final int CYP = 9;
    protected static final int CAR = 10;
    protected static final int MER = 11;
    protected static final int CEA = 12;
    protected static final int COP = 13;
    protected static final int COD = 14;
    protected static final int COE = 15;
    protected static final int COO = 16;
    protected static final int BON = 17;
    protected static final int PCO = 18;
    protected static final int GLS = 19;
    protected static final int PAR = 20;
    protected static final int AIT = 21;
    protected static final int MOL = 22;
    protected static final int CSC = 23;
    protected static final int QSC = 24;
    protected static final int TSC = 25;

    static BitSet whichProj(String[] args) {
	BitSet out = new BitSet(nprojs);
	int i, j;
	String[] pcode = { "NON", "AZP", "TAN", "SIN", "STG", "ARC", 
			   "ZPN", "ZEA", "AIR", "CYP", "CAR", 
			   "MER", "CEA", "COP", "COD", "COE", 
			   "COO", "BON", "PCO", "GLS", "PAR", 
			   "AIT", "MOL", "CSC", "QSC", "TSC"  }; 

	if (args == null || args.length == 0) {

	    // do all projections
	    for(i=0; i < nprojs; i++) out.set(i);
	}
	else {
	    for(j=0; j < args.length; j++) {
		args[j].trim();
		for(i=0; i < nprojs; i++) {
		    if (args[j].equalsIgnoreCase(pcode[i])) {
			out.set(i);
			continue;
		    }
		}
	    }
	}
		
	return out;
    }
}

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
import FITSWCS.projections.*;
import java.util.BitSet;
import Acme.Fmt;

/**
 *   This class verifies the Projection class for closure errors.<p>
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
public class TestProj {

    public static void main(String args[]) {

	double tol = 1.0e-10;
	Projection prj;
	double[] p = new double[10];

	System.out.println("Testing closure of WCSLIB coordinate " + 
			   "transformation routines");
	System.out.println("-------------------------------------" + 
			   "-----------------------");

	for (int j = 0; j < 10; p[j++] = 0.0);

	BitSet doproj = whichProj(args);

	if (doproj.get(AZP)) {
	    
	    // AZP: zenithal/azimuthal perspective.
	    p[1] = 2.0;
	    try {
		prj = new AZPProjection(p);
//  	        prj = new AZPProjection(p[1]);
		runTest("AZP", prj, 90, -30, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }

	} 
	if (doproj.get(TAN)) {

	    // TAN: gnomonic.
//	    prj = new TANProjection(p);
	    prj = new TANProjection();
	    runTest("TAN", prj, 90, 5, tol);

	} 
	if (doproj.get(SIN)) {

	    // SIN: orthographic/synthesis
	    p[1] = 0.3;
	    p[2] = 1.5;
	    prj = new SINProjection(p);
//	    prj = new SINProjection();
	    runTest("SIN", prj, 90, 60, tol);

	} 
	if (doproj.get(STG)) {

	    // STG: stereographic.
	    prj = new STGProjection(p);
//	    prj = new STGProjection();
	    runTest("STG", prj, 90, -85, tol);

	} 
	if (doproj.get(ARC)) {

	    // ARC: zenithal/azimuthal equidistant.
	    prj = new ARCProjection(p);
//	    prj = new ARCProjection();
	    runTest("ARC", prj, 90, -90, tol);

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
		prj = new ZPNProjection(p);
		runTest("ZPN", prj, 90, 10, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(ZEA)) {

	    // ZEA: zenithal/azimuthal equal area. 
	    prj = new ZEAProjection(p);
//	    prj = new ZEAProjection();
	    runTest("ZEA", prj, 90, -90, tol);

	}
	if (doproj.get(AIR)) {
	    
	    // AIR: Airy's zenithal projection. 
	    p[1] = 45.0;
	    try {
		prj = new AIRProjection(p);
//  	        prj = new AIRProjection(p[1]);
		runTest("AIR", prj, 90, -85, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(CYP)) {
	    
	    // CYP: cylindrical perspective.
	    p[1] = 3.0;
	    p[2] = 0.8;
	    try {
		prj = new CYPProjection(p);
//  	        prj = new CYPProjection(p[1], p[2]);
		runTest("CYP", prj, 90, -90, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(CAR)) {

	    // CAR: Cartesian. 
//	    prj = new CARProjection(p);
	    prj = new CARProjection();
	    runTest("CAR", prj, 90, -90, tol);

	}
	if (doproj.get(MER)) {

	    // MER: Cartesian. 
//	    prj = new MERProjection(p);
	    prj = new MERProjection();
	    runTest("MER", prj, 85, -85, tol);

	}
	if (doproj.get(CEA)) {
	    
	    // CEA: cylindrical equal area.
	    p[1] = 0.75;
	    try {
//		prj = new CEAProjection(p);
  	        prj = new CEAProjection(p[1]);
		runTest("CEA", prj, 90, -90, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(COP)) {
	    
	    // COP: conic perspective.
	    p[1] = 60.0;
	    p[2] = 15.0;
	    try {
//		prj = new COPProjection(p);
  	        prj = new COPProjection(p[1], p[2]);
		runTest("COP", prj, 90, -25, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(COD)) {
	    
	    // COD: conic equidistant
	    p[1] = -60.0;
	    p[2] = 15.0;
	    try {
//		prj = new CODProjection(p);
  	        prj = new CODProjection(p[1], p[2]);
		runTest("COD", prj, 90, -90, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	} 
	if (doproj.get(COE)) {
	    
	    // COE: conic equal area.
	    p[1] = 60.0;
	    p[2] = -15.0;
	    try {
//		prj = new COEProjection(p);
  	        prj = new COEProjection(p[1], p[2]);
		runTest("COE", prj, 90, -90, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(COO)) {
	    
	    // COO: conic orthomorphic.
	    p[1] = -60.0;
	    p[2] = -15.0;
	    try {
//		prj = new COOProjection(p);
  	        prj = new COOProjection(p[1], p[2]);
		runTest("COO", prj, 85, -90, tol);
	    }
	    catch (BadProjectionParameterException ex) {
		System.out.println(ex.getMessage());
	    }
	}
	if (doproj.get(BON)) {

	    // BON: Bonne's projection.
	    p[1] = 30.0;
	    prj = new BONProjection(p);
//	    prj = new BONProjection(p[1]);
//	    p[1] = 0.0;
//	    prj = new BONProjection(p);
	    runTest("BON", prj, 90, -90, tol);

	}
	if (doproj.get(PCO)) {

	    // PCO: polyconic.
	    prj = new PCOProjection(p);
//	    prj = new PCOProjection();
	    runTest("PCO", prj, 90, -90, tol);

	}
	if (doproj.get(GLS)) {

	    // GLS: Sanson-Flamsteed (global sinusoid). 
	    prj = new GLSProjection(p);
//	    prj = new GLSProjection();
	    runTest("GLS", prj, 90, -90, tol);

	}
	if (doproj.get(PAR)) {

	    // PAR: parabolic. 
	    prj = new PARProjection(p);
//	    prj = new PARProjection();
	    runTest("PAR", prj, 90, -90, tol);

	}
	if (doproj.get(AIT)) {

	    // AIT: Hammer-Aitoff.
	    prj = new AITProjection(p);
//	    prj = new AITProjection();
	    runTest("AIT", prj, 90, -90, tol);

	}
	if (doproj.get(MOL)) {

	    // MOL: Mollweide's projection. 
	    prj = new MOLProjection(p);
//	    prj = new MOLProjection();
	    runTest("MOL", prj, 90, -90, tol);

	}
	if (doproj.get(CSC)) {

	    // CSC: COBE quadrilateralized spherical cube.
	    prj = new CSCProjection(p);
//	    prj = new CSCProjection();
	    runTest("CSC", prj, 90, -90, 4.0e-2);

	}
	if (doproj.get(QSC)) {

	    // QSC: quadrilateralized spherical cube.
	    prj = new QSCProjection(p);
//	    prj = new QSCProjection();
	    runTest("QSC", prj, 90, -90, tol);

	}
	if (doproj.get(TSC)) {

	    // TSC: tangential spherical cube.
	    prj = new TSCProjection(p);
//	    prj = new TSCProjection();
	    runTest("TSC", prj, 90, -90, tol);

	}  
    }

    protected static void runTest(String pcode, Projection prj, 
				  int north, int south, double tol) 
    {
	int err, lat, lng, j;
	double dlat, dlatmx, dlng, dlngmx, dr, drmax, lat1, lat2, lng1, lng2;
	double r, theta, x, x1, x2, y, y1, y2;
	double[] out;

	System.out.printf("Testing %s latitudes %s to %s, ",
                          pcode, north, south)
                  .printf("closure tolerance %8.1e deg.\n", tol);

	dlngmx = 0.0;
	dlatmx = 0.0;

	for (lat = north; lat >= south; lat--) {
	    lat1 = (double)lat;

	    for (lng = -180; lng <= 180; lng++) {
		lng1 = (double)lng;

		try {
		    out = prj.fwd(lng1, lat1);
		}
		catch (PixelBeyondProjectionException ex) {
		    System.out.printf("Error: lng1 = %20.15f  lat = %20.15f\n",
                                       lng1, lat1);
		    System.out.println("       " + ex.getMessage());
		    continue;
		}

		x = out[0];
		y = out[1];
		try {
		    out = prj.rev(x, y);
		}
		catch (PixelBeyondProjectionException ex) {
		    System.out.printf("Error: lng1 = %20.15f  lat1 = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("          x = %20.15f     y = %20.15f\n",
                                      x, y);
		    System.out.println("       " + ex.getMessage());
		    continue;
		}
		lng2 = out[0];
		lat2 = out[1];

		dlng = Math.abs(lng2-lng1);
		if (dlng > 180.0) dlng = Math.abs(dlng-360.0);
		if (Math.abs(lat) != 90 && dlng > dlngmx) dlngmx = dlng;
		dlat = Math.abs(lat2-lat1);
		if (dlat > dlatmx) dlatmx = dlat;

		if (dlat > tol) {
		    System.out.printf("  lng = %20.15f  lat = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("  %s:     x = %20.15f     y = %20.15f\n",
                                      pcode, x, y);
		    System.out.printf("     : lng2 = %20.15f  lat2 = %20.15f\n",
                                      lng2, lat2);
		} 
		else if (Math.abs(lat) != 90) {
		    if (dlng > tol) {
		    System.out.printf("  lng = %20.15f  lat = %20.15f\n",
                                      lng1, lat1);
		    System.out.printf("  %s:     x = %20.15f     y = %20.15f\n",
                                      pcode, x, y);
		    System.out.printf("     : lng2 = %20.15f  lat2 = %20.15f\n",
                                      lng2, lat2);
		    } 
		}
	    }
	}

	System.out.println("  Maximum residual (sky): lng: " + 
			   Fmt.fmt(dlngmx, 20, 15) + "   lat: " +
			   Fmt.fmt(dlatmx, 20, 15));

	// Test closure at a point close to the reference point. 
	r = 10.0;
	theta = -195.0;

	drmax = 0.0;

	for (j = 1; j <= 12; j++) {
	    r /= 10.0;
	    theta += 15.0;

	    x1 = r*TrigD.cos(theta);
	    y1 = r*TrigD.sin(theta);

	    try {
		out = prj.rev(x1, y1);
	    }
	    catch (PixelBeyondProjectionException ex) {
		System.out.printf("Error: (r,th)=(%8.4f, %8.4f) " +
                                  "x1 = %20.15f  y1 = %20.15f\n",
                                  r, theta, x1, y1);
		System.out.println("       " + ex.getMessage());
		continue;
	    }
	    lng1 = out[0];
	    lat1 = out[1];

	    try {
		out = prj.fwd(lng1, lat1);
	    }
	    catch (PixelBeyondProjectionException ex) {
		System.out.printf("Error: (r,th)=(%8.4f, %8.4f) " +
                                  "x1 = %20.15f  y1 = %20.15f\n",
                                  r, theta, x1, y1);
		System.out.printf("       lng1 = %20.15f  lat1 = %20.15f\n",
                                  lng1, lat1); 
		System.out.println("       " + ex.getMessage());
		continue;
	    }
	    x2 = out[0];
	    y2 = out[1];

	    dr = Math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	    if (dr > drmax) drmax = dr;
	    if (dr > tol) {
                System.out.printf("  %s:  x1 = %20.15f   y1 = %20.15f\n",
                                  pcode, x1, y1);
                System.out.printf("     : lng = %20.15f  lat = %20.15f\n",
                                  lng1, lat1);
                System.out.printf("     :  x2 = %20.15f   y2 = %20.15f\n",
                                  x2, y2);
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

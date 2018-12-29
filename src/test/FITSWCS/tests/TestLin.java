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

import java.util.Formatter;

/**
 *   This class verifies the LinearTransform class for closure errors.<p>
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
public class TestLin {

    public static void main(String args[]) {

	double[] crpix =  {256.0, 256.0,  64.0, 128.0,   1.0};
	double[][] pc = {{  1.0,   0.5,   0.0,   0.0,   0.0},
			 {  0.5,   1.0,   0.0,   0.0,   0.0},
			 {  0.0,   0.0,   1.0,   0.0,   0.0},
			 {  0.0,   0.0,   0.0,   1.0,   0.0},
			 {  0.0,   0.0,   0.0,   0.0,   1.0}};
	double[] cdelt = {  1.2,   2.3,   3.4,   4.5,   5.6};
	double[] pix =   {303.0, 265.0, 112.4, 144.5,  28.2};
	double[] img;

	int j;
	StringBuffer out;
        Formatter fmtr;
	LinearTransform lin;

	try {
	    lin = new LinearTransform(5, crpix, pc, cdelt);
	} catch (SingularMatrixException ex) {
	    throw new InternalError(ex.getMessage());
	}

	System.out.println("Testing closure of WCSLIB coordinate " + 
			   "transformation routines");
	System.out.println("-------------------------------------" + 
			   "-----------------------");

	out = new StringBuffer("pix:");
        fmtr = new Formatter(out);
	for (j=0; j < 5; j++) {
	    fmtr.format("%14.8f", pix[j]);
	}
	System.out.println(out);

	img = lin.rev(pix);
	out = new StringBuffer("img:");
        fmtr = new Formatter(out);
	try {
	    for (j=0; j < 5; j++) {
                fmtr.format("%14.8f", img[j]);
	    }
	} catch (ArrayIndexOutOfBoundsException ex) {
	    System.out.println("Missing elements in img");
	}
	System.out.println(out);

	pix = lin.fwd(img);
	out = new StringBuffer("pix:");
        fmtr = new Formatter(out);
	try {
	    for (j=0; j < 5; j++) {
                fmtr.format("%14.8f", pix[j]);
	    }
	} catch (ArrayIndexOutOfBoundsException ex) {
	    System.out.println("Missing elements in pix");
	}
	System.out.println(out);

	img = lin.rev(pix);
	out = new StringBuffer("img:");
        fmtr = new Formatter(out);
	try {
	    for (j=0; j < 5; j++) {
                fmtr.format("%14.8f", img[j]);
	    }
	} catch (ArrayIndexOutOfBoundsException ex) {
	    System.out.println("Missing elements in img");
	}
	System.out.println(out);

    }
}

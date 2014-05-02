/*===========================================================================
 * 
 * Adapted for JUnit from FITSWCS/tests/TestLin.java, a translation of tlin.c
 * from WCSlib (Copyright (C) 1995,1996 Mark Calabretta).
 *===========================================================================*/

package FITSWCS.tests;

import FITSWCS.LinearTransform;
import FITSWCS.exceptions.SingularMatrixException;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import static org.junit.Assert.*;

public class LinClosureTest {

    static final double[] crpix =  {256.0, 256.0,  64.0, 128.0,   1.0};
    static final double[][] pc = {{  1.0,   0.5,   0.0,   0.0,   0.0},
                                  {  0.5,   1.0,   0.0,   0.0,   0.0},
                                  {  0.0,   0.0,   1.0,   0.0,   0.0},
                                  {  0.0,   0.0,   0.0,   1.0,   0.0},
                                  {  0.0,   0.0,   0.0,   0.0,   1.0}};
    static final double[] cdelt = {  1.2,   2.3,   3.4,   4.5,   5.6};
    static final double[] pix =   {303.0, 265.0, 112.4, 144.5,  28.2};
    static final double[] img =   { 61.8, 74.75, 164.56, 74.25, 152.32 };
    static final double tol = 0.5e-8;

    LinearTransform lin;

    @Before
    public void setup() {
        lin = null;
    }

    @Test
    public void testCtor() {
	try {
	    lin = new LinearTransform(5, crpix, pc, cdelt);
	} catch (SingularMatrixException ex) {
	    fail(ex.getMessage());
	}
    }

    @Test
    public void testClosure() {
        testCtor();

        double[] trx = lin.rev(pix);
        assertArrayEquals(img, trx, tol);

        trx = lin.fwd(trx);
        assertArrayEquals(pix, trx, tol);

        trx = lin.rev(trx);
        assertArrayEquals(img, trx, tol);

        trx = lin.fwd(trx);
        assertArrayEquals(pix, trx, tol);
    }

    public static void main(String args[]) {
      org.junit.runner.JUnitCore.main("dalserver.conf.NTFITSHeaderKeywordsTest");
    }
}
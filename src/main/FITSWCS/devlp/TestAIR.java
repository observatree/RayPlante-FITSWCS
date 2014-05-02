package FITSWCS;

import java.util.BitSet;
import Acme.Fmt;

public class TestAIR {

    public static void main(String args[]) {

	int j;
	int    err, lat, lng;
	double dlat, dlatmx, dlng, dlngmx, dr, drmax, lat1, lat2, lng1, lng2;
	double r, theta, x, x1, x2, y, y1, y2;
	String pcode = "AIR";
	int north=90, south=-85;
	double tol = 1.0e-10;
//	double[] p = new double[10];
	double[] out = new double[2];
	Projection prj;

	System.out.println("Testing " + pcode + " latitudes " + north +
			   " to " + south + ", closure tolerance " + 
			   Fmt.fmt(tol, 8, 1) + " deg.");

	try {
	    prj = new AIRProjection((double)45.0);
	} catch (BadProjectionParameterException ex) {
	    throw new InternalError(ex.getMessage());
	}

	dlngmx = 0.0;
	dlatmx = 0.0;

	/* Test closure at a point close to the reference point. */
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
		System.out.println("Error: (r,th)=(" + r + ", " + theta +
				   ") x1 =" + 
				   Fmt.fmt(x1, 20, 15) +
				   "  y1 =" + 
				   Fmt.fmt(y1, 20, 15));
		System.out.println("       " + ex.getMessage());
		continue;
	    }
	    lng1 = out[0];
	    lat1 = out[1];

	    try {
		out = prj.fwd(lng1, lat1);
	    }
	    catch (PixelBeyondProjectionException ex) {
		System.out.println("Error: (r,th)=(" + r + ", " + theta +
				   ") x1 =" + 
				   Fmt.fmt(x1, 20, 15) +
				   "  y1 =" + 
				   Fmt.fmt(y1, 20, 15));
		System.out.println("       lng1 =" + 
				   Fmt.fmt(lng1, 20, 15) +
				   "  lat =" + 
				   Fmt.fmt(lat1, 20, 15));
		System.out.println("       " + ex.getMessage());
		continue;
	    }
	    x2 = out[0];
	    y2 = out[1];

	    dr = Math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	    if (dr > drmax) drmax = dr;
	    if (dr > tol) {
		    System.out.println("  " + pcode + ": x1 =" +
				       Fmt.fmt(x1, 20, 15) +
				       "  y1 =" + 
				       Fmt.fmt(y1, 20, 15));
		    System.out.println("     : lng =" +
				       Fmt.fmt(lng1, 20, 15) +
				       "  lat =" + 
				       Fmt.fmt(lat1, 20, 15));
		    System.out.println("     : x2 =" +
				       Fmt.fmt(x2, 20, 15) +
				       "  y2 =" + 
				       Fmt.fmt(y2, 20, 15));
	    }
	}

	System.out.println("  Maximum residual (map): dR: " + 
			   Fmt.fmt(drmax, 10, 31));

    }
}

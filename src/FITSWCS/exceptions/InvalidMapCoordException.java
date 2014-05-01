/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS.exceptions;

/**
 * an exception thrown when an invalid position is given to a method of a
 * CelestialTransform object.
 */
public class InvalidMapCoordException 
    extends InvalidCelestialTransformException 
{
    public InvalidMapCoordException(String s) { super(s); }
    public InvalidMapCoordException() { super(); }

    public InvalidMapCoordException(String pcode, double x, double y) 
    {
	super("Invalid map position (x, y) for " + pcode + 
	      "projection: " + x + ", " + y);
    }
}

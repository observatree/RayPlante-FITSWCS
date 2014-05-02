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
 * A given pixel lies beyond the defined spherical projection or the 
 *   conversion diverges at that position.
 */
public class PixelBeyondProjectionException extends ProjectionException {
    public PixelBeyondProjectionException(String s) { super(s); }
    public PixelBeyondProjectionException() { 
	super("Conversion of input pixel diverges in this projection."); }
}

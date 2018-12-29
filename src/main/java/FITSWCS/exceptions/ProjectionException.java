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
 * an exception associated with use of a Projection object.
 */
public class ProjectionException extends FITSWCSException {
    public ProjectionException(String s) { super(s); }
    public ProjectionException() { super(); }
}

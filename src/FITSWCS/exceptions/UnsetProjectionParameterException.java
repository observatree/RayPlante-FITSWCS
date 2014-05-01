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

/*
 * an exception indicating that a call to Projection object's fwd() or 
 * rev() method was made without first setting the projection parameters
 * (either in the Projection object's constructor or with the setProjParm()
 * method).  
 */
public class UnsetProjectionParameterException extends RuntimeException {
    public UnsetProjectionParameterException(String s) {
        super(s);
    }

    public UnsetProjectionParameterException() {
        super("Projection parameters not set.");
    }
}

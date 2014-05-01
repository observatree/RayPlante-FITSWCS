/*============================================================================
*
*   FITSWCS - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995,1996 Mark Calabretta
*   Translated into Java(TM) from WCSLIB (C impl) 8/1996
*   by Raymond L. Plante, copyright (c) 1996
*
*   $Id: wcstrig.c,v 2.1 1996/05/07 20:05:10 mcalabre Exp $
*===========================================================================*/

package FITSWCS;

import FITSWCS.exceptions.*;
import java.util.*;

/**
 * an interface containing codes and indexes for various projection
 * types
 */
public interface ProjectionType {

    public static final int NTYPES = 26;

    public static final int NON =  0;
    public static final int AZP =  1;
    public static final int TAN =  2;
    public static final int SIN =  3;
    public static final int STG =  4;
    public static final int ARC =  5;
    public static final int ZPN =  6;
    public static final int ZEA =  7;
    public static final int AIR =  8;
    public static final int CYP =  9;
    public static final int CAR = 10;
    public static final int MER = 11;
    public static final int CEA = 12;
    public static final int COP = 13;
    public static final int COD = 14;
    public static final int COE = 15;
    public static final int COO = 16;
    public static final int BON = 17;
    public static final int PCO = 18;
    public static final int GLS = 19;
    public static final int PAR = 20;
    public static final int AIT = 21;
    public static final int MOL = 22;
    public static final int CSC = 23;
    public static final int QSC = 24;
    public static final int TSC = 25;

    public String[] code = { "NON", "AZP", "TAN", "SIN", "STG", "ARC", 
			     "ZPN", "ZEA", "AIR", "CYP", "CAR", 
			     "MER", "CEA", "COP", "COD", "COE", 
			     "COO", "BON", "PCO", "GLS", "PAR", 
			     "AIT", "MOL", "CSC", "QSC", "TSC"  }; 
}

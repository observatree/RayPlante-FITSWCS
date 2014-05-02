// Util - assorted static utility routines
//
// Copyright (C) 1996 by Jef Poskanzer <jef@acme.com>.  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// Visit the ACME Labs Java page for up-to-date versions of this and other
// fine Java utilities: http://www.acme.com/java/

package Acme;

import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.applet.*;
import java.net.*;

/// Assorted static utility routines.
// <P>
// Whenever I come up with a static routine that might be of general use,
// I put it here.  So far the class includes:
// <UL>
// <LI> some string routines that were left out of java.lang.String
// <LI> a color-spec parser
// <LI> a full-precision replacement for Double.toString()
// <LI> a fixed version of java.io.InputStream's byte-array read routine
// <LI> a standard broken-image icon
// <LI> a thick-line drawing routine
// </UL>
// and lots more.
// <P>
// <A HREF="/resources/classes/Acme/Util.java">Fetch the software.</A><BR>
// <A HREF="/resources/classes/Acme.tar.Z">Fetch the entire Acme package.</A>

public class Util
    {

    /// Returns a date string formatted in Unix ls style - if it's within
    // six months of now, Mmm dd hh:ss, else Mmm dd  yyyy.
    @SuppressWarnings("deprecation")
    public static String lsDateStr( Date date )
        {
        long dateTime = date.getTime();
	if ( dateTime == -1L )
	    return "------------";
        long nowTime = (new Date()).getTime();
        String[] months = {
            "Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
        String part1 = months[date.getMonth()] + Fmt.fmt( date.getDate(), 3 );
        if ( Math.abs( nowTime - dateTime ) < 183L * 24L * 60L * 60L * 1000L )
            return part1 + Fmt.fmt( date.getHours(), 3 ) + ":" +
                Fmt.fmt( date.getMinutes(), 2, Fmt.ZF );
        else
            return part1 + Fmt.fmt( date.getYear() + 1900, 6 );
        }


    /// Returns the length of the initial segment of str which consists
    // entirely of characters from charSet.
    public static int strSpan( String str, String charSet )
	{
	int i;
	for ( i = 0; i < str.length(); ++i )
	    if ( charSet.indexOf( str.charAt( i ) ) == -1 )
		break;
	return i;
	}

    /// Returns the length of the initial segment of str which consists
    // entirely of characters NOT from charSet.
    public static int strCSpan( String str, String charSet )
	{
	int i;
	for ( i = 0; i < str.length(); ++i )
	    if ( charSet.indexOf( str.charAt( i ) ) != -1 )
		break;
	return i;
	}

    /// Returns the length of the initial segment of str1 and str2 that matches.
    public static int sameLength( String str1, String str2 )
	{
	int i;
	for ( i = 0;
	      i < str1.length() && i < str2.length() &&
		str1.charAt( i ) == str2.charAt( i );
	      ++i )
	    ;
	return i;
	}

    /// Returns the number of times the given character appears in the string.
    public static int charCount( String str, char c )
	{
	int n = 0;
	for ( int i = 0; i < str.length(); ++i )
	    if ( str.charAt( i ) == c )
		++n;
	return n;
	}


    /// Turns a String into an array of Strings, by using StringTokenizer
    // to split it up at whitespace.
    public static String[] splitStr( String str )
	{
	StringTokenizer st = new StringTokenizer( str );
	int n = st.countTokens();
	String[] strs = new String[n];
	for ( int i = 0; i < n; ++i )
	    strs[i] = st.nextToken();
	return strs;
	}


    /// Turns an array of Strings into a single String, with the components
    // separated by spaces.
    public static String flattenStrarr( String[] strs )
	{
	StringBuffer sb = new StringBuffer();
	for ( int i = 0; i < strs.length; ++i )
	    {
	    if ( i > 0 )
		sb.append( ' ' );
	    sb.append( strs[i] );
	    }
	return sb.toString();
	}


    /// Returns the number a raised to the power of b.  Long version
    // of Math.pow().  Throws ArithmeticException if b is negative.
    public static long pow( long a, long b )
	throws ArithmeticException
	{
	if ( b < 0 )
	    throw new ArithmeticException();
	long r = 1;
	while ( b != 0 )
	    {
	    if ( odd( b ) )
		r *= a;
	    b >>>= 1;
	    a *= a;
	    }
	return r;
	}


    /// Improved version of Double.toString(), returns more decimal places.
    public static String doubleToString( double d )
	{
	// As of JDK 1.0, Double.toString() is a native method that just
	// does an sprintf with a %g.  At least on some systems, this returns
	// only six decimal places.  This replacement version gives the full
	// sixteen digits.

	// Handle special numbers first, to avoid complications.
	if ( Double.isNaN( d ) )
	    return "NaN";
	if ( d == Double.NEGATIVE_INFINITY )
	    return "-Inf";
	if ( d == Double.POSITIVE_INFINITY )
	    return "Inf";

	// Grab the sign, and then make the number positive for simplicity.
	boolean negative = false;
	if ( d < 0.0D )
	    {
	    negative = true;
	    d = -d;
	    }

	// Get the native version of the unsigned value, as a template.
	String unsStr = Double.toString( d );

	// Dissect out the exponent.
	String mantStr, expStr;
	int exp;
	int eInd = unsStr.indexOf( 'e' );
	if ( eInd == -1 )
	    {
	    mantStr = unsStr;
	    expStr = "";
	    exp = 0;
	    }
	else
	    {
	    mantStr = unsStr.substring( 0, eInd );
	    expStr = unsStr.substring( eInd + 1 );
	    if ( expStr.startsWith( "+" ) )
		exp = Integer.parseInt( expStr.substring( 1 ) );
	    else
		exp = Integer.parseInt( expStr );
	    }

	// Dissect out the number part.
	String numStr;
	int dotInd = mantStr.indexOf( '.' );
	if ( dotInd == -1 )
	    numStr = mantStr;
	else
	    numStr = mantStr.substring( 0, dotInd );
	long num;
	if ( numStr.length() == 0 )
	    num = 0;
	else
	    num = Integer.parseInt( numStr );

	// Build the new mantissa.
	StringBuffer newMantBuf = new StringBuffer( numStr + "." );
	double p = Math.pow( 10, exp );
	double frac = d - num * p;
	String digits = "0123456789";
	int nDigits = 16 - numStr.length();	// about 16 digits in a double
	for ( int i = 0; i < nDigits; ++i )
	    {
	    p /= 10.0D;
	    int dig = (int) ( frac / p );
	    if ( dig < 0 ) dig = 0;
	    if ( dig > 9 ) dig = 9;
	    newMantBuf.append( digits.charAt( dig ) );
	    frac -= dig * p;
	    }

	if ( (int) ( frac / p + 0.5D ) == 1 )
	    {
	    // Round up.
	    boolean roundMore = true;
	    for ( int i = newMantBuf.length() - 1; i >= 0; --i )
		{
		int dig = digits.indexOf( newMantBuf.charAt( i ) );
		if ( dig == -1 )
		    continue;
		++dig;
		if ( dig == 10 )
		    {
		    newMantBuf.setCharAt( i, '0' );
		    continue;
		    }
		newMantBuf.setCharAt( i, digits.charAt( dig ) );
		roundMore = false;
		break;
		}
	    if ( roundMore )
		{
		// If this happens, we need to prepend a 1.  But I haven't
		// found a test case yet, so I'm leaving it out for now.
		// But if you get this message, please let me know!
		newMantBuf.append( "ROUNDMORE" );
		}
	    }

	// Chop any trailing zeros.
	int len = newMantBuf.length();
	while ( newMantBuf.charAt( len - 1 ) == '0' )
	    newMantBuf.setLength( --len );
	// And chop a trailing dot, if any.
	if ( newMantBuf.charAt( len - 1 ) == '.' )
	    newMantBuf.setLength( --len );

	// Done.
	return ( negative ? "-" : "" ) +
	       newMantBuf +
	       ( expStr.length() != 0 ? ( "e" + expStr ) : "" );
	}


    private final static int bWidth = 20;
    private final static int bHeight = 19;
    private static int[] bPixels = {	// color model is AARRGGBB
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0x00ffffff,
	0x00ffffff, 0x00ffffff, 0x00ffffff, 0x00ffffff, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,
	0x00ffffff, 0x00ffffff, 0x00ffffff, 0x00ffffff, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xffffffff,
	0xff000000, 0x00ffffff, 0x00ffffff, 0x00ffffff, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000,
	0xff00ff00, 0xff00ff00, 0xff00ff00, 0xff00ff00, 0xff000000,
	0xff000000, 0xffff00ff, 0xffff00ff, 0xff000000, 0xffffffff,
	0xffffffff, 0xff000000, 0x00ffffff, 0x00ffffff, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff00ff00,
	0xff00ff00, 0xff000000, 0xff000000, 0xff00ff00, 0xff00ff00,
	0xff000000, 0xff000000, 0xffff00ff, 0xff000000, 0xffffffff,
	0xffffffff, 0xffffffff, 0xff000000, 0x00ffffff, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff00ff00,
	0xff000000, 0xff000000, 0xffff00ff, 0xff000000, 0xff00ff00,
	0xff00ff00, 0xff000000, 0xffff00ff, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0x00ffffff,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000,
	0xff000000, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff00ff00,
	0xff00ff00, 0xff000000, 0xffff00ff, 0xffff00ff, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xff000000, 0xff00ff00, 0xff00ff00,
	0xff000000, 0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xff000000, 0xff00ff00, 0xff00ff00, 0xff000000,
	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xff000000, 0xff00ff00, 0xff00ff00, 0xff000000, 0xff000000,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xff000000, 0xff00ff00, 0xff000000, 0xff000000, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xff000000, 0xff000000, 0xff000000, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xff000000, 0xff000000, 0xff000000, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xff000000, 0xff00ff00, 0xff00ff00, 0xff000000, 0xff000000,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xff000000, 0xff00ff00, 0xff000000, 0xff000000, 0xff000000,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xff000000, 0xff000000, 0xff000000, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff, 0xffff00ff,
	0xffff00ff, 0xffff00ff, 0xffff00ff, 0xff000000, 0xff000000,

	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,

	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	0xff000000, 0xff000000, 0xff000000, 0xff000000, 0xff000000,
	};

    /// Draw a broken-image image.
    @SuppressWarnings("deprecation")
    public static void brokenImage( Graphics graphics, Component comp )
	{
	Image bimg = comp.createImage(
	    new MemoryImageSource( bWidth, bHeight, bPixels, 0, bWidth ) );
	Dimension d = comp.size();
	graphics.setColor( comp.getBackground() );
	graphics.fillRect( 0, 0, d.width, d.height );
	graphics.drawImage(
	    bimg, ( d.width - bWidth ) / 2, ( d.height - bHeight ) / 2, null );
	}


    /// Parse a color string into a Color.  The color can be specified
    // by name as one of:
    // <BLOCKQUOTE>
    // black blue cyan darkGray gray green lightGray
    // magenta orange pink red white yellow
    // </BLOCKQUOTE>
    // Or, as an #rrggbb hex number, like in Netscape.
    public static Color parseColor( String str )
	{
	if ( str.startsWith( "#" ) )
	    {
	    try
		{
		int h = Integer.parseInt( str.substring( 1 ), 16 );
		return new Color(
		    ( h >>> 16 ) & 0xff, ( h >>> 8 ) & 0xff, h & 0xff );
		}
	    catch ( NumberFormatException e )
		{
		return null;
		}
	    }
	Color color;
	color = parseNamedColor( str );
	if ( color != null )
	    return color;
	if ( str.substring( 0, 4 ).equalsIgnoreCase( "dark" ) )
	    {
	    color = parseNamedColor( str.substring( 4 ) );
	    if ( color != null )
		return color.darker();
	    }
	if ( str.substring( 0, 5 ).equalsIgnoreCase( "light" ) )
	    {
	    color = parseNamedColor( str.substring( 5 ) );
	    if ( color != null )
		return color.brighter();
	    }
	if ( str.substring( 0, 6 ).equalsIgnoreCase( "bright" ) )
	    {
	    color = parseNamedColor( str.substring( 6 ) );
	    if ( color != null )
		return color.brighter();
	    }
	return null;
	}

    private static Color parseNamedColor( String str )
	{
	if ( str.equalsIgnoreCase( "black" ) )
	    return Color.black;
	if ( str.equalsIgnoreCase( "blue" ) )
	    return Color.blue;
	if ( str.equalsIgnoreCase( "cyan" ) )
	    return Color.cyan;
	if ( str.equalsIgnoreCase( "darkGray" ) )
	    return Color.darkGray;
	if ( str.equalsIgnoreCase( "gray" ) )
	    return Color.gray;
	if ( str.equalsIgnoreCase( "green" ) )
	    return Color.green;
	if ( str.equalsIgnoreCase( "lightGray" ) )
	    return Color.lightGray;
	if ( str.equalsIgnoreCase( "magenta" ) )
	    return Color.magenta;
	if ( str.equalsIgnoreCase( "orange" ) )
	    return Color.orange;
	if ( str.equalsIgnoreCase( "pink" ) )
	    return Color.pink;
	if ( str.equalsIgnoreCase( "red" ) )
	    return Color.red;
	if ( str.equalsIgnoreCase( "white" ) )
	    return Color.white;
	if ( str.equalsIgnoreCase( "yellow" ) )
	    return Color.yellow;
	return null;
	}


    /// Handle the standard BGCOLOR parameter.  Call as:
    // <BLOCKQUOTE>
    // Acme.Util.handleBgcolor( this );
    // </BLOCKQUOTE>
    // at the start of your init() method.
    public static void handleBgcolor( Applet applet )
	{
	String param = applet.getParameter( "bgcolor" );
	if ( param != null )
	    {
	    Color color = parseColor( param );
	    if ( color != null )
		applet.setBackground( color );
	    }
	}


    /// Test is a number is even.
    public static boolean even( long n )
	{
	return ( n & 1 ) == 0;
	}

    /// Test is a number is odd.
    public static boolean odd( long n )
	{
	return ( n & 1 ) != 0;
	}


    /// Count the number of 1-bits in a byte.
    public static int countOnes( byte n )
	{
	return countOnes( n & 0xffL );
	}

    /// Count the number of 1-bits in an int.
    public static int countOnes( int n )
	{
	return countOnes( n & 0xffffffffL );
	}

    /// Count the number of 1-bits in a long.
    public static int countOnes( long n )
	{
	// There are faster ways to do this, all the way up to looking
	// up bytes in a 256-element table.  But this is not too bad.
	int count = 0;
	while ( n != 0 )
	    {
	    if ( odd( n ) )
		++count;
	    n >>>= 1;
	    }
	return count;
	}


    /// A fixed version of java.io.InputStream.read(byte[], int, int).  The
    // standard version catches and ignores IOExceptions from below.
    // This version sends them on to the caller.
    public static int fixedRead( InputStream in, byte[] b, int off, int len ) throws IOException
        {
        if ( len <= 0 )
            return 0;
        int c = in.read();
        if ( c == -1 )
            return -1;
        if ( b != null )
            b[off] = (byte) c;
        int i;
        for ( i = 1; i < len ; ++i )
            {
            c = in.read();
            if ( c == -1 )
                break;
            if ( b != null )
                b[off + i] = (byte) c;
            }
        return i;
        }


    private static final int SPLINE_THRESH = 3;

    /// Draw a three-point spline.
    public static void drawSpline( Graphics graphics, int x1, int y1, int x2, int y2, int x3, int y3 )
	{
	int xa, ya, xb, yb, xc, yc, xp, yp;

	xa = ( x1 + x2 ) / 2;
	ya = ( y1 + y2 ) / 2;
	xc = ( x2 + x3 ) / 2;
	yc = ( y2 + y3 ) / 2;
	xb = ( xa + xc ) / 2;
	yb = ( ya + yc ) / 2;

	xp = ( x1 + xb ) / 2;
	yp = ( y1 + yb ) / 2;
	if ( Math.abs( xa - xp ) + Math.abs( ya - yp ) > SPLINE_THRESH )
	    drawSpline( graphics, x1, y1, xa, ya, xb, yb );
	else
	    graphics.drawLine( x1, y1, xb, yb );

	xp = ( x3 + xb ) / 2;
	yp = ( y3 + yb ) / 2;
	if ( Math.abs( xc - xp ) + Math.abs( yc - yp ) > SPLINE_THRESH )
	    drawSpline( graphics, xb, yb, xc, yc, x3, y3 );
	else
	    graphics.drawLine( xb, yb, x3, y3 );
	}


    private static final int DDA_SCALE = 8192;

    /// Draw a thick line.
    public static void drawThickLine( Graphics graphics, int x1, int y1, int x2, int y2, int linewidth )
	{
	// Draw the starting point filled.
	graphics.fillOval(
	    x1 - linewidth / 2, y1 - linewidth / 2, linewidth, linewidth );

	// Short-circuit zero-length lines.
	if ( x1 == x2 && y1 == y2 )
	    return;

	/* Draw, using a simple DDA. */
	if ( Math.abs( x2 - x1 ) > Math.abs( y2 - y1 ) )
	    {
	    // Loop over X domain.
	    int dy, srow;
	    int dx, col, row, prevrow;

	    if ( x2 > x1 )
		dx = 1;
	    else
		dx = -1;
	    dy = ( y2 - y1 ) * DDA_SCALE / Math.abs( x2 - x1 );
	    prevrow = row = y1;
	    srow = row * DDA_SCALE + DDA_SCALE / 2;
	    col = x1;
	    for (;;)
		{
		if ( row != prevrow )
		    {
		    graphics.drawOval(
			col - linewidth / 2, prevrow - linewidth / 2,
			linewidth, linewidth );
		    prevrow = row;
		    }
		graphics.drawOval(
		    col - linewidth / 2, row - linewidth / 2,
		    linewidth, linewidth );
		if ( col == x2 )
		    break;
		srow += dy;
		row = srow / DDA_SCALE;
		col += dx;
		}
	    }
	else
	    {
	    // Loop over Y domain.
	    int dx, scol;
	    int dy, col, row, prevcol;

	    if ( y2 > y1 )
		dy = 1;
	    else
		dy = -1;
	    dx = ( x2 - x1 ) * DDA_SCALE / Math.abs( y2 - y1 );
	    row = y1;
	    prevcol = col = x1;
	    scol = col * DDA_SCALE + DDA_SCALE / 2;
	    for ( ; ; )
		{
		if ( col != prevcol )
		    {
		    graphics.drawOval(
			prevcol - linewidth / 2, row - linewidth / 2,
			linewidth, linewidth );
		    prevcol = col;
		    }
		graphics.drawOval(
		    col - linewidth / 2, row - linewidth / 2,
		    linewidth, linewidth );
		if ( row == y2 )
		    break;
		row += dy;
		scol += dx;
		col = scol / DDA_SCALE;
		}
	    }
	}


    /// Make a URL with no ref part and no query string.  Also, if it's
    // a directory then make sure there's a trailing slash.
    public static URL plainUrl( URL context, String urlStr )
	throws MalformedURLException
	{
	URL url = new URL( context, urlStr );
	String fileStr = url.getFile();
	int i = fileStr.indexOf( '?' );
	if ( i != -1 )
	    fileStr = fileStr.substring( 0, i );
	url = new URL(
	    url.getProtocol(), url.getHost(), url.getPort(), fileStr );
	if ( ( ! fileStr.endsWith( "/" ) ) &&
	     urlStrIsDir( url.toExternalForm() ) )
	    {
	    fileStr = fileStr + "/";
	    url = new URL(
		url.getProtocol(), url.getHost(), url.getPort(), fileStr );
	    }
	return url;
	}

    /// Make a URL with no ref part and no query string.  Also, if it's
    // a directory then make sure there's a trailing slash.
    public static URL plainUrl( String urlStr )
	throws MalformedURLException
	{
        return plainUrl( null, urlStr );
        }

    /// Figure out the base URL for a given URL.  What this means is
    // if the URL points to a directory, you get that directory; if the
    // URL points to a file, you get the directory the file is in.
    public static String baseUrlStr( String urlStr )
	{
	if ( urlStr.endsWith( "/" ) )
	    return urlStr;
	if ( urlStrIsDir( urlStr ) )
	    return urlStr + "/";
	return urlStr.substring( 0, urlStr.lastIndexOf( '/' ) + 1 );
	}

    /// Makes sure if a URL is a directory, it ends with a slash.
    public static String fixDirUrlStr( String urlStr )
	{
	if ( urlStr.endsWith( "/" ) )
	    return urlStr;
	if ( urlStrIsDir( urlStr ) )
	    return urlStr + "/";
	return urlStr;
	}

    /// Figures out whether a URL points to a directory or not.
    // Web servers are lenient and accept directory-URLs without
    // the trailing slash.  What they actually do is return a
    // redirect to the same URL with the trailing slash appended.
    // Unfortunately, Java doesn't let us see that such a redirect
    // happened.  Instead we have to figure out it's a directory
    // indirectly and heuristically.
    public static boolean urlStrIsDir( String urlStr )
	{
	// If it ends with a slash, it's probably a directory.
	if ( urlStr.endsWith( "/" ) )
	    return true;

	// If the last component has a dot, it's probably not a directory.
	int lastSlash = urlStr.lastIndexOf( '/' );
	int lastPeriod = urlStr.lastIndexOf( '.' );
	if ( lastPeriod != -1 && ( lastSlash == -1 || lastPeriod > lastSlash ) )
	    return false;

	// Otherwise, append a slash and try to connect.  This is
	// fairly expensive.
	String urlStrWithSlash = urlStr + "/";
	try
	    {
	    URL url = new URL( urlStrWithSlash );
	    InputStream f = url.openStream();
	    f.close();
	    // Worked fine - it's probably a directory.
	    return true;
	    }
	catch ( Exception e )
	    {
	    // Got an error - must not be a directory.
	    return false;
	    }
	}


    // Figures out whether a URL is absolute or not.
    public static boolean urlStrIsAbsolute( String urlStr )
	{
	if ( urlStr.startsWith( "/" ) || urlStr.indexOf( ":/" ) != -1 )
	    return true;
	// Should handle :8000/ and such too.
	return false;
	}

    // Returns an equivalent URL string that is guaranteed to be absolute.
    public static String absoluteUrlStr( String urlStr, URL contextUrl ) throws MalformedURLException
	{
	URL url = new URL( contextUrl, urlStr );
	return url.toExternalForm();
	}


    /// Check if an array contains a given element.
    public static boolean arraycontains( Object[] array, Object element )
	{
	for ( int i = 0; i < array.length; ++i )
	    if ( array[i].equals( element ) )
		return true;
	return false;
	}


    /// Run a program on the host system.
    // <P>
    // This routine runs the specified command, waits for it to
    // finish, and returns the exit status.
    // This is like the Unix system() routine.  Unlike the Unix version,
    // though, stdout and stderr get thrown away unless you redirect them.
    public static int system( String cmd )
	{
	try
	    {
	    Runtime runtime = Runtime.getRuntime();
	    String[] shCmd = new String[3];
	    shCmd[0] = "/bin/sh";
	    shCmd[1] = "-c";
	    shCmd[2] = cmd;
	    Process process = runtime.exec( shCmd );
	    return process.waitFor();
	    }
	catch ( IOException e )
	    {
	    return -1;
	    }
	catch ( InterruptedException e )
	    {
	    return -1;
	    }
	}

    }

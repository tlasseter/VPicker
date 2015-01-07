package com.Sts.Utilities;

import com.magician.fonts.*;

import javax.media.opengl.*;

/**
 * <p>Title: Vertical Font</p>
 * <p>Description: Develop fonts for vertical axis</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

public class StsVerticalFont extends GLBitmapFont
{
    private GLBitmapFont font;
    private int firstChar;
    private int numChars;
    private GLBitmapChar[] fontChars;

    static final boolean debug = false;

    public StsVerticalFont(GLBitmapFont font, GL gl)
    {
		super(gl);
        GLBitmapChar c;

        try
        {
            this.font = font;
            firstChar = font.getFirstChar();
            setFirstChar(firstChar);
            numChars = font.getNumChars();
            setNumChars(numChars);

            fontChars = new GLBitmapChar[numChars];

            GLBitmapChar[] horizontalChars = font.getFontData();
            for(int n = 0; n < numChars; n++)
            {
                c = horizontalChars[n];
                if(c == null) continue;
                int width = c.width;
                int height = c.height;
                float x = c.xorig;
                float y = c.yorig;
                float advance = c.advance;
                byte[] bitmap = c.bitmap;

                int nBits = width*height;
                byte[] verticalBitmap = swapBitmap(c.bitmap, width, height);
                fontChars[n] = new GLBitmapChar(c.height, c.width, (float)(c.height-c.yorig), c.xorig, c.advance, verticalBitmap);

                if(debug)
                {
                    printChar(c, n, "Horizontal");
                    printChar(fontChars[n], n, "Vertical");
                }
            }
            setFontData(fontChars);
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalFont.constructor() failed.",
                e, StsException.WARNING);
        }
    }

    /** a GLBitmapChar is width by height bytes.  Each row consists of 8 bits per byte
     *  padded out to byte boundary, e.g., if the width is 8 and the height is 12,
     *  there are 12 rows with 2 bytes per row.
     *  swapBitmap reverses rows and columns, so each new row has the same number of
     *  bits as the original column, padded out to byte boundary. So in our example
     *  above, we would now have 5 rows consisting of 2 bytes per row.
     *  Because the character is rotated counterclockwise, the row is the reverse
     *  of the original column with the bits beginning at the top of the original column
     *  padded out as needed below the original column.
     */
    protected byte[] swapBitmap(byte[] inBytes, int width, int height)
    {
        int nBits = width*height;
        if(nBits == 0) return inBytes;  // zero-length character

        int nBytesPerRow = width/8;
        if(nBytesPerRow*8 < width) nBytesPerRow++;
        int nBytesPerCol = height/8;
        if(nBytesPerCol*8 < height) nBytesPerCol++;
        int nBitsPerRow = nBytesPerRow*8;
        int nBitsPerCol = nBytesPerCol*8;

        // nRows = nBitsPerCol including padding out to byte boundary
        // nCols = nBitsPerRow including padding out to byte boundary
        int[][] bits = new int[nBitsPerCol][nBitsPerRow];
        int nByte = 0;
        for(int row = 0; row < height; row++)
        {
            bits[row] = getRowBits(inBytes, nByte, nBytesPerRow);
            nByte += nBytesPerRow;
        }

        byte[] swapBytes = new byte[nBytesPerCol*width];
        int nSwapByte = 0;
        int[] rowBits = new int[nBitsPerCol];
        for(int row = 0; row < width; row++)
        {
            int swapCol = 0;
            for(int col = height-1; col >= 0; col--)
                rowBits[swapCol++] = bits[col][row];
            setRowBits(rowBits, nBytesPerCol, swapBytes, nSwapByte);
            nSwapByte += nBytesPerCol;
        }
        return swapBytes;
    }

    private int[] getRowBits(byte[] inBytes, int firstByte, int nBytesPerRow)
    {
        int[] rowBits = new int[nBytesPerRow*8];

        int bb = firstByte;
        int nBit = 0;
        for(int b = 0; b < nBytesPerRow; b++)
        {
            int i = (int)inBytes[bb++] & 0xFF;
            if(debug)
                System.out.println("integer: " + i + " binaryString: " + Integer.toBinaryString(i));
            for(int s = 7; s >= 0; s--)
                rowBits[nBit++] = i >> s & 1;

            if(debug)
                System.out.println("inRowByteBits: " + rowBits[nBit-8] + rowBits[nBit-7] + rowBits[nBit-6] + rowBits[nBit-5] +
                                                     + rowBits[nBit-4] + rowBits[nBit-3] + rowBits[nBit-2] + rowBits[nBit-1]);
        }

        return rowBits;
    }

    public void setRowBits(int[] rowBits, int nRowBytes, byte[] rowBytes, int nByte)
    {
        int nBit = 0;
        for(int b = 0; b < nRowBytes; b++)
        {
            int i = 0;
            if(debug)
                System.out.println("outRowByteBits: " + rowBits[nBit] + rowBits[nBit+1] + rowBits[nBit+2] + rowBits[nBit+3] +
                                                      + rowBits[nBit+4] + rowBits[nBit+5] + rowBits[nBit+6] + rowBits[nBit+7]);
            for(int s = 7; s >= 0; s--)
                i = rowBits[nBit++] << s | i;
            rowBytes[nByte++] = (byte)i;
            if(debug) System.out.println("rowByte: " + i);
        }
    }

    /** Overrides drawCharacter() in superclass GLBitmapFont with only difference
     *  that advance is vertical rather than horizontal
     */
    public void drawCharacter(GL gl, char charToDraw)
    {
//        System.err.println( "drawCharacter()" );

        int[] swapbytes = new int[1];
        int[] lsbfirst = new int[1];
        int[] rowlength = new int[1];
        int[] skiprows = new int[1];
        int[] skippixels = new int[1];
        int[] alignment = new int[1];

//        System.err.println( "GLBitmapFont::charToDraw: " + (int)charToDraw );

        /** Check char is valid */
        if ( (int)charToDraw < firstChar ||
             (int)charToDraw >= firstChar + numChars ||
             fontChars[charToDraw - firstChar] == null ) {
            System.err.println( "Char outwith range" );
          }

        /** Save current modes */
        gl.glGetIntegerv( GL.GL_UNPACK_SWAP_BYTES, swapbytes, 0 );
        gl.glGetIntegerv( GL.GL_UNPACK_LSB_FIRST, lsbfirst, 0 );
        gl.glGetIntegerv( GL.GL_UNPACK_ROW_LENGTH, rowlength, 0 );
        gl.glGetIntegerv( GL.GL_UNPACK_SKIP_ROWS, skiprows, 0 );
        gl.glGetIntegerv( GL.GL_UNPACK_SKIP_PIXELS, skippixels, 0 );
        gl.glGetIntegerv( GL.GL_UNPACK_ALIGNMENT, alignment, 0 );

        /* Little endian machines (DEC Alpha for example) could
           benefit from setting GL_UNPACK_LSB_FIRST to GL_TRUE
           instead of GL_FALSE, but this would require changing the
           generated bitmaps too. */
        gl.glPixelStorei( GL.GL_UNPACK_SWAP_BYTES, GL.GL_FALSE );
        gl.glPixelStorei( GL.GL_UNPACK_LSB_FIRST, GL.GL_FALSE );
        gl.glPixelStorei( GL.GL_UNPACK_ROW_LENGTH, 0 );
        gl.glPixelStorei( GL.GL_UNPACK_SKIP_ROWS, 0 );
        gl.glPixelStorei( GL.GL_UNPACK_SKIP_PIXELS, 0 );
        gl.glPixelStorei( GL.GL_UNPACK_ALIGNMENT, 1 );
        gl.glBitmap( fontChars[charToDraw - firstChar].width,
                      fontChars[charToDraw - firstChar].height,
                      fontChars[charToDraw - firstChar].xorig,
                      fontChars[charToDraw - firstChar].yorig,
                      0,
                      fontChars[charToDraw - firstChar].advance,
                      fontChars[charToDraw - firstChar].bitmap, 0 );

        /* Restore saved modes. */
        gl.glPixelStorei( GL.GL_UNPACK_SWAP_BYTES, swapbytes[0] );
        gl.glPixelStorei( GL.GL_UNPACK_LSB_FIRST, lsbfirst[0] );
        gl.glPixelStorei( GL.GL_UNPACK_ROW_LENGTH, rowlength[0] );
        gl.glPixelStorei( GL.GL_UNPACK_SKIP_ROWS, skiprows[0] );
        gl.glPixelStorei( GL.GL_UNPACK_SKIP_PIXELS, skippixels[0] );
        gl.glPixelStorei( GL.GL_UNPACK_ALIGNMENT, alignment[0] );
    }

    /**
     * Draws the specified string in the current font
     * <P>
     * @param stringToDraw The string to draw
     */
    public final void drawString(GL gl, String stringToDraw )
    {

        /** Convert the string into an array of characters */
        char[] tmpArray = stringToDraw.toCharArray();

        /** Iterate through the string and call drawCharacter() */
        for ( int i = 0 ; i < tmpArray.length ; i++ )
            drawCharacter(gl, tmpArray[i] );
    }

    /**
     * Draws the specified string in the current font
     * <P>
     * @param stringToDraw The string to draw
     */
    public final void drawString(GL gl,  char[] stringToDraw )
    {
        for ( int i = 0 ; i < stringToDraw.length ; i++ )
            drawCharacter(gl, stringToDraw[i] );
    }

    /**
     * Draws the specified string in the current font
     * <P>
     * @param stringToDraw The string to draw
     */
    public final void drawString(GL gl, byte[] stringToDraw )
    {
        for ( int i = 0 ; i < stringToDraw.length ; i++ )
            drawCharacter(gl, (char)stringToDraw[i] );
    }

    private void printChar(GLBitmapChar c, int n, String orientation)
    {
        System.out.println("Character " + n + orientation + "\n");

        int width = c.width;
        int height = c.height;
        float x = c.xorig;
        float y = c.yorig;
        float advance = c.advance;
        byte[] bitmap = c.bitmap;

        int nBits = width*height;
        if(nBits == 0)
        {
            System.out.println("  zero length character");
            return;
        }
        int nBytesPerRow = width/8;
        if(nBytesPerRow*8 < width) nBytesPerRow++;
        int nBitsPerRow = nBytesPerRow*8;

        int[] rowBits = new int[nBitsPerRow];

        for(int row = height-1; row >= 0; row--)
        {
            int nByte = row*nBytesPerRow;
            rowBits = getRowBits(bitmap, nByte, nBytesPerRow);
            StringBuffer stringBuffer = new StringBuffer(nBitsPerRow+1);
            for(int b = 0; b < nBitsPerRow; b++)
            {
                if(rowBits[b] == 1) stringBuffer.append('X');
                else stringBuffer.append('.');
            }
            String string = stringBuffer.toString();
            System.out.println(string);
        }
        System.out.println("\n");
    }
}


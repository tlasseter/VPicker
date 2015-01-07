package com.Sts.Utilities;

import com.Sts.DBTypes.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 27, 2008
 * Time: 8:23:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsBitmap
{
    private int nRows;
    private int nCols;
    private byte[] bytes;
    private float pixelXOrigin;
    private float pixelYOrigin;

    public static final byte SMALL = 0;
    public static final byte BIG = 1;

    // WELL HEAD BITMAPS
    static public byte DRY = 0;
    static public byte GAS = 1;
    static public byte OIL = 2;
    static public byte GASOIL = 3;
    static public byte INJECTOR = 4;
    static public byte OBSERVATION = 5;
    static public byte[] wellTypes = {DRY, GAS, OIL, GASOIL, INJECTOR, OBSERVATION};
    static String[] wellTypeStrings = {"Dry", "Gas", "Oil", "Gas & Oil", "Injector", "Observation" };
    static public StsBitmap dryWellBitmap = new StsBitmap( new String[]
                                                                      { "00010000",
                                                                        "00111000",
                                                                        "01000100",
                                                                        "11000011",
                                                                        "01000010",
                                                                        "01000100",
                                                                        "00111000",
                                                                        "00010000"});
    static public StsBitmap gasWellBitmap = new StsBitmap( new String[]
                                                                     { "00010000",
                                                                       "10111010",
                                                                       "01000100",
                                                                       "11000011",
                                                                       "01000010",
                                                                       "01000100",
                                                                       "10111010",
                                                                       "00010000"});
    static public StsBitmap observationWellBitmap = new StsBitmap( new String[]
                                                                      { "00000000",
                                                                        "00111000",
                                                                        "11000111",
                                                                        "01000010",
                                                                        "01000010",
                                                                        "11000111",
                                                                        "00111000",
                                                                        "00000000"});
    static public StsBitmap oilWellBitmap = new StsBitmap( new String[]
                                                                      { "00000000",
                                                                        "00111000",
                                                                        "01111100",
                                                                        "01111110",
                                                                        "01111110",
                                                                        "01111100",
                                                                        "00111000",
                                                                        "00000000"});
    static public StsBitmap oilAndGasWellBitmap = new StsBitmap( new String[]
                                                                      { "00010000",
                                                                        "10111010",
                                                                        "01111100",
                                                                        "11111111",
                                                                        "01111110",
                                                                        "01111100",
                                                                        "10111010",
                                                                        "00010000"});
    static public StsBitmap injectorBitmap = new StsBitmap( new String[]
                                                                      { "00011000",
                                                                        "00011000",
                                                                        "00111100",
                                                                        "00111100",
                                                                        "01111110",
                                                                        "01111110",
                                                                        "11111111",
                                                                        "11111111"});
    static public StsBitmap[] wellBitmaps = {dryWellBitmap, gasWellBitmap, oilWellBitmap, oilAndGasWellBitmap, injectorBitmap, observationWellBitmap};

    // EQUIPMENT BITMAPS
    static public byte PACKER = StsEquipmentMarker.PACKER;
    static public byte CASING = StsEquipmentMarker.CASING;
    static public byte SENSOR = StsEquipmentMarker.SENSOR;
    static public byte[] equipTypes = StsEquipmentMarker.subTypes;
    static String[] equipTypesStrings = StsEquipmentMarker.subTypeStrings;
    static public StsBitmap packerBitmap = new StsBitmap( new String[]
                                                                 { "00011000",
                                                                   "00111100",
                                                                   "00111100",
                                                                   "01111110",
                                                                   "01111110",
                                                                   "01111110",
                                                                   "01111110",
                                                                   "11111111"});
    static public StsBitmap casingBitmap = new StsBitmap( new String[]
                                                                     { "00000000",
                                                                       "00011000",
                                                                       "01111110",
                                                                       "11111111",
                                                                       "11111111",
                                                                       "01111110",
                                                                       "00011000",
                                                                       "00000000"});
    static public StsBitmap sensorBitmap = new StsBitmap( new String[]
                                                                     { "00000001",
                                                                       "00000111",
                                                                       "00001111",
                                                                       "11111111",
                                                                       "11111111",
                                                                       "00001111",
                                                                       "00000111",
                                                                       "00000001"});
    static public StsBitmap[] equipBitmaps = {packerBitmap, casingBitmap, sensorBitmap};

    // GENERAL SHAPE BITMAPS
    static public StsBitmap triangleBitmap = new StsBitmap( new String[]
                            { "0000000100000000",
                              "0000001110000000",
                              "0000011111000000",
                              "0000111111100000",
                              "0001111111110000",
                              "0011111111111000",
                              "0111111111111100",
                              "1111111111111110" }
                            );
    static public StsBitmap smallTriangleBitmap = new StsBitmap( new String[]
                                                                       { "00010000",
                                                                         "00111000",
                                                                         "01111100",
                                                                         "11111110"}
                                                                       );
    static public StsBitmap fontBitmapF = new StsBitmap( new String[]
                            { "1111111100000000",
                              "1111111100000000",
                              "1100000000000000",
                              "1100000000000000",
                              "1100000000000000",
                              "1111111000000000",
                              "1111111000000000",
                              "1100000000000000",
                              "1100000000000000",
                              "1100000000000000",
                              "1100000000000000",
                              "1100000000000000" }
                            );

       static public StsBitmap diamondBitmap = new StsBitmap( new String[]
                                                                       { "0000000100000000",
                                                                         "0000001110000000",
                                                                         "0000011111000000",
                                                                         "0000111111100000",
                                                                         "0001111111110000",
                                                                         "0011111111111000",
                                                                         "0001111111110000",
                                                                         "0000111111100000",
                                                                         "0000011111000000",
                                                                         "0000001110000000" }
                                                                       );

       static public StsBitmap smallDiamondBitmap = new StsBitmap( new String[]
                                                                         { "00010000",
                                                                           "00111000",
                                                                           "01111100",
                                                                           "11111110",
                                                                           "01111100",
                                                                           "00111000",
                                                                           "00010000",
                                                                           "00000000"}
                                                                         );

       static public StsBitmap starBitmap = new StsBitmap( new String[]
                                                                       { "0000000100000000",
                                                                         "0000001110000000",
                                                                         "0000011111000000",
                                                                         "0000111111100000",
                                                                         "1111111111111111",
                                                                         "0111111111111110",
                                                                         "0011111111111100",
                                                                         "0001111111111000",
                                                                         "0001111101111100",
                                                                         "0011110000111100",
                                                                         "0111000000001110",
                                                                         "1100000000000011" }
                                                                       );
       static public StsBitmap smallStarBitmap = new StsBitmap( new String[]
                                                                              { "00000000",
    		   																	"01000100",
                                                                                "00101000",
                                                                                "01111100",
                                                                                "00101000",
                                                                                "01000100",
                                                                                "00000000",
                                                                                "00000000"}
                                                                              );
       static public StsBitmap circleBitmap = new StsBitmap( new String[]
                                                                      { "0000011111100000",
                                                                        "0001111111111000",
                                                                        "0011111111111100",
                                                                        "0011111111111100",
                                                                        "0011111111111100",
                                                                        "0011111111111100",
                                                                        "0001111111111000",
                                                                        "0000011111100000"
                                                                        }
                                                                      );

       static public StsBitmap smallCircleBitmap = new StsBitmap( new String[]
                                                                           { "00111000",
                                                                             "01111100",
                                                                             "01111100",
                                                                             "00111000"
                                                                             }
                                                                           );

    static public final boolean debug = false;

    public StsBitmap(String[] stringRowBits)
    {
        bytes = bitMap2dToBytes(stringRowBits);
    }
    static void drawCircle(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallCircleBitmap.glBitmap(gl, xyz, color);
    	else
    		circleBitmap.glBitmap(gl, xyz, color);
    }
    static void drawCircle2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallCircleBitmap.glBitmap2d(gl, xyz, color);
    	else
    		circleBitmap.glBitmap2d(gl, xyz, color);
    }
    static void drawTriangle(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallTriangleBitmap.glBitmap(gl, xyz, color);
    	else
    		triangleBitmap.glBitmap(gl, xyz, color);
    }
    static void drawTriangle2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallTriangleBitmap.glBitmap2d(gl, xyz, color);
    	else
    		triangleBitmap.glBitmap2d(gl, xyz, color);
    }
    static void drawStar(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallStarBitmap.glBitmap(gl, xyz, color);
    	else
    		starBitmap.glBitmap(gl, xyz, color);
    }
    static void drawStar2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallStarBitmap.glBitmap2d(gl, xyz, color);
    	else
    		starBitmap.glBitmap2d(gl, xyz, color);
    }
    static void drawDiamond(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallDiamondBitmap.glBitmap(gl, xyz, color);
    	else
    		diamondBitmap.glBitmap(gl, xyz, color);
    }
    static void drawDiamond2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type == SMALL)
    		smallDiamondBitmap.glBitmap2d(gl, xyz, color);
    	else
            diamondBitmap.glBitmap2d(gl, xyz, color);
    }
    static void drawF(GL gl, float[] xyz, StsColor color, byte type)
    {
        fontBitmapF.glBitmap(gl, xyz, color);
    }
    static void drawBitmap(GL gl, float[] xyz, StsBitmap bitmap, StsColor color)
    {
    	bitmap.glBitmap(gl, xyz, color);
    }
    static boolean drawEquipmentBitmap(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type >= equipTypes.length)
    		return false;

    	equipBitmaps[type].glBitmap(gl, xyz, color);
    	return true;
    }
    static boolean drawEquipmentBitmap2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type >= equipTypes.length)
    		return false;

    	equipBitmaps[type].glBitmap2d(gl, xyz, color);
    	return true;
    }
    static boolean drawEquipmentBitmap2d(GL gl, int x, int y, StsColor color, byte type)
    {
    	if(type >= equipTypes.length)
    		return false;

    	equipBitmaps[type].glBitmap2d(gl, x, y, color);
    	return true;
    }
    static boolean drawWellBitmap(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type >= wellTypes.length)
    		return false;

    	wellBitmaps[type].glBitmap(gl, xyz, color);
    	return true;
    }
    static boolean drawWellBitmap2d(GL gl, float[] xyz, StsColor color, byte type)
    {
    	if(type >= wellTypes.length)
    		return false;

    	wellBitmaps[type].glBitmap2d(gl, xyz, color);
    	return true;
    }
    public void glBitmap(GL gl, float[] xyz, StsColor color)
    {
        color.setGLColor(gl);
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
        gl.glRasterPos3f(xyz[0], xyz[1], xyz[2]);
        gl.glBitmap(nCols, nRows, pixelXOrigin, pixelYOrigin, 0.0f, 0.0f, bytes, 0);
        gl.glFlush();
    }

    public void glBitmap2d(GL gl, float[] xyz, StsColor color)
    {
        color.setGLColor(gl);
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
        gl.glRasterPos2f(xyz[0], xyz[1]);
        gl.glBitmap(nCols, nRows, pixelXOrigin, pixelYOrigin, 0.0f, 0.0f, bytes, 0);
        gl.glFlush();
    }

    public void glBitmap2d(GL gl, int x, int y, StsColor color)
    {
        color.setGLColor(gl);
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
        gl.glRasterPos2i(x, y);
        gl.glBitmap(nCols, nRows, pixelXOrigin, pixelYOrigin, 0.0f, 0.0f, bytes, 0);
        // jbw this is wrong never use 4   gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 4);
        gl.glFlush();
    }

    private  byte[] bitMap2dToBytes(String[] stringBits)
    {
        nRows = stringBits.length;
        nCols = stringBits[0].length();
        pixelXOrigin = nRows/2;
        pixelYOrigin = nCols/2;
        int nRowBytes = nCols/8;
        int nBytes = nRowBytes*nRows;
        byte[] bytes = new byte[nBytes];
        int nByte = 0;
        for(int row = nRows-1; row >= 0; row--)
        {
            char[] chars = stringBits[row].toCharArray();
            int nRowBit = 0;
            if(debug) System.out.println("row " + row + ": " + stringBits[row]);
            for(int rowByte = 0; rowByte < nRowBytes; rowByte++)
           {
                int i = 0;
                for(int s = 7; s >= 0; s--)
                {
                    char c = chars[nRowBit++];
                    int ii;
                    if(c == '0')
                        ii = 0;
                    else
                        ii = 1;
                    i = ii << s | i;
                }
                bytes[nByte++] = (byte)i;
                if(debug) System.out.println("rowByte: " + i);
            }
        }
        return bytes;
    }
}

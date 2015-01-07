package com.Sts.MVC.View3d;

import com.Sts.DBTypes.StsSeismicVelocityModel;
import com.Sts.Interfaces.*;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Types.*;
import com.Sts.Utilities.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import java.nio.*;

/**
 * Textures consists of n by n square texels.  The center of the first texel
 * is 1/(2*n) and the center of the last is 1-1/(2*n).  Each tile butts against
 * the edge of the adjoining tile, i.e., they share texels along that common edge.
 */

public class StsTextureTile
{
    int nTotalRows, nTotalCols;
    int nTotalPoints;
    public int rowMin, rowMax, colMin, colMax;
    public int croppedRowMin, croppedRowMax, croppedColMin, croppedColMax;
    public int nRows, nCols;
    StsTextureTiles textureTiles;
    int nBackgroundRows, nBackgroundCols;
    public double minRowTexCoor, minColTexCoor;
    public double maxRowTexCoor, maxColTexCoor;
    public double dRowTexCoor, dColTexCoor;
    public double[][] xyzPlane = new double[4][];
    int texture = 0;
    ByteBuffer tileData;
    int displayListNum = 0;
    int displayListNum2d = 0;
    boolean axesFlipped = false;
    int nTile = -1;

    static private final int XDIR = StsCursor3d.XDIR;
    static private final int YDIR = StsCursor3d.YDIR;
    static private final int ZDIR = StsCursor3d.ZDIR;

    static private boolean runTimer = false;
    static private StsTimer timer = null;
    static public boolean debug = false;

    public StsTextureTile(int nTotalRows, int nTotalCols, int row, int col, StsTextureTiles textureTiles, int maxTextureSize, int nTile)
    {
        this.nTotalRows = nTotalRows;
        this.nTotalCols = nTotalCols;
        this.nTotalPoints = nTotalRows * nTotalCols;
        this.rowMin = row;
        this.colMin = col;
        this.nTile = nTile;
        this.textureTiles = textureTiles;

        if(runTimer && timer == null) timer = new StsTimer();

        rowMax = Math.min(row + maxTextureSize - 1, nTotalRows - 1);
        colMax = Math.min(col + maxTextureSize - 1, nTotalCols - 1);

        initializeCroppedRange();

        nRows = rowMax - rowMin + 1;
        nCols = colMax - colMin + 1;

        nBackgroundRows = StsMath.nextBaseTwoInt(nRows);
        nBackgroundCols = StsMath.nextBaseTwoInt(nCols);

        double rowBorder = 1.0 / (2.0 * nBackgroundRows);
        double colBorder = 1.0 / (2.0 * nBackgroundCols);

        minRowTexCoor = rowBorder;
        minColTexCoor = colBorder;
        maxRowTexCoor = (double) nRows / nBackgroundRows - rowBorder;
        maxColTexCoor = (double) nCols / nBackgroundCols - colBorder;
        dRowTexCoor = (double) (maxRowTexCoor - minRowTexCoor) / (nRows - 1);
        dColTexCoor = (double) (maxColTexCoor - minColTexCoor) / (nCols - 1);

        int dirNo = textureTiles.dir;
        if(dirNo < 0) return;

        double minRowF = (double) rowMin / (nTotalRows - 1);
        double minColF = (double) colMin / (nTotalCols - 1);
        double maxRowF = (double) rowMax / (nTotalRows - 1);
        double maxColF = (double) colMax / (nTotalCols - 1);

        double xTileMin = 0.0f, xTileMax = 0.0f;
        double yTileMin = 0.0f, yTileMax = 0.0f;
        double zTileMin = 0.0f, zTileMax = 0.0f;

        StsRotatedGridBoundingBox boundingBox = textureTiles.boundingBox;
        switch(dirNo)
        {
            case StsCursor3d.XDIR:
                yTileMin = boundingBox.yMin + minRowF * (boundingBox.yMax - boundingBox.yMin);
                yTileMax = boundingBox.yMin + maxRowF * (boundingBox.yMax - boundingBox.yMin);
                zTileMin = boundingBox.zMin + minColF * (boundingBox.zMax - boundingBox.zMin);
                zTileMax = boundingBox.zMin + maxColF * (boundingBox.zMax - boundingBox.zMin);
                xyzPlane[0] = new double[]
                        {0.0f, yTileMin, zTileMin};
                xyzPlane[1] = new double[]
                        {0.0f, yTileMin, zTileMax};
                xyzPlane[2] = new double[]
                        {0.0f, yTileMax, zTileMax};
                xyzPlane[3] = new double[]
                        {0.0f, yTileMax, zTileMin};
                break;
            case StsCursor3d.YDIR:
                xTileMin = boundingBox.xMin + minRowF * (boundingBox.xMax - boundingBox.xMin);
                xTileMax = boundingBox.xMin + maxRowF * (boundingBox.xMax - boundingBox.xMin);
                zTileMin = boundingBox.zMin + minColF * (boundingBox.zMax - boundingBox.zMin);
                zTileMax = boundingBox.zMin + maxColF * (boundingBox.zMax - boundingBox.zMin);
                xyzPlane[0] = new double[]
                        {xTileMin, 0.0f, zTileMin};
                xyzPlane[1] = new double[]
                        {xTileMin, 0.0f, zTileMax};
                xyzPlane[2] = new double[]
                        {xTileMax, 0.0f, zTileMax};
                xyzPlane[3] = new double[]
                        {xTileMax, 0.0f, zTileMin};
                break;
            case StsCursor3d.ZDIR:
                xTileMin = boundingBox.xMin + minColF * (boundingBox.xMax - boundingBox.xMin);
                xTileMax = boundingBox.xMin + maxColF * (boundingBox.xMax - boundingBox.xMin);
                yTileMin = boundingBox.yMin + minRowF * (boundingBox.yMax - boundingBox.yMin);
                yTileMax = boundingBox.yMin + maxRowF * (boundingBox.yMax - boundingBox.yMin);
                xyzPlane[0] = new double[]
                        {xTileMin, yTileMin, 0.0f};
                xyzPlane[1] = new double[]
                        {xTileMax, yTileMin, 0.0f};
                xyzPlane[2] = new double[]
                        {xTileMax, yTileMax, 0.0f};
                xyzPlane[3] = new double[]
                        {xTileMin, yTileMax, 0.0f};
                break;
        }

        if(debug)
        {
            System.out.println("StsTextureTile.constructor() called.");
            System.out.println("    Geometry for tile in direction: " + dirNo);
            System.out.println("        xTileMin: " + xTileMin + " xTileMax: " + xTileMax);
            System.out.println("        yTileMin: " + yTileMin + " yTileMax: " + yTileMax);
            System.out.println("        zTileMin: " + zTileMin + " zTileMax: " + zTileMax);
            System.out.println("        maxRowTexCoor: " + maxRowTexCoor + " maxColTexCoor: " + maxColTexCoor);
        }
    }

    public StsTextureTile(int nTotalRows, int nTotalCols, int row, int col, StsTextureTiles textureTiles, int maxTextureSize, float[][] axisRanges, int nTile)
    {
        this.nTotalRows = nTotalRows;
        this.nTotalCols = nTotalCols;
        this.nTotalPoints = nTotalRows * nTotalCols;
        this.rowMin = row;
        this.colMin = col;
        this.textureTiles = textureTiles;
        this.nTile = nTile;

        if(runTimer && timer == null) timer = new StsTimer();

        rowMax = Math.min(row + maxTextureSize - 1, nTotalRows - 1);
        colMax = Math.min(col + maxTextureSize - 1, nTotalCols - 1);

        initializeCroppedRange();

        nRows = rowMax - rowMin + 1;
        nCols = colMax - colMin + 1;

        nBackgroundRows = StsMath.nextBaseTwoInt(nRows);
        nBackgroundCols = StsMath.nextBaseTwoInt(nCols);

        double rowBorder = 1.0 / (2.0 * nBackgroundRows);
        double colBorder = 1.0 / (2.0 * nBackgroundCols);

        minRowTexCoor = rowBorder;
        minColTexCoor = colBorder;
        maxRowTexCoor = (double) nRows / nBackgroundRows - rowBorder;
        maxColTexCoor = (double) nCols / nBackgroundCols - colBorder;
        dRowTexCoor = (double) (maxRowTexCoor - minRowTexCoor) / (nRows - 1);
        dColTexCoor = (double) (maxColTexCoor - minColTexCoor) / (nCols - 1);

        int dirNo = textureTiles.dir;
        if(dirNo < 0) return;

        double minRowF = (double) rowMin / (nTotalRows - 1);
        double minColF = (double) colMin / (nTotalCols - 1);
        double maxRowF = (double) rowMax / (nTotalRows - 1);
        double maxColF = (double) colMax / (nTotalCols - 1);

        double xTileMin = 0.0f, xTileMax = 0.0f;
        double yTileMin = 0.0f, yTileMax = 0.0f;
        double zTileMin = 0.0f, zTileMax = 0.0f;

        float zMin = axisRanges[0][0];
        float zMax = axisRanges[0][1];
        float xMin = axisRanges[1][0];
        float xMax = axisRanges[1][1];

        xTileMin = xMin + minRowF * (xMax - xMin);
        xTileMax = xMin + maxRowF * (xMax - xMin);
        zTileMin = zMin + minColF * (zMax - zMin);
        zTileMax = zMin + maxColF * (zMax - zMin);
        xyzPlane[0] = new double[]
                {xTileMin, 0.0f, zTileMin};
        xyzPlane[1] = new double[]
                {xTileMin, 0.0f, zTileMax};
        xyzPlane[2] = new double[]
                {xTileMax, 0.0f, zTileMax};
        xyzPlane[3] = new double[]
                {xTileMax, 0.0f, zTileMin};
        if(debug)
        {
            System.out.println("StsTextureTile.constructor() called.");
            System.out.println("    Geometry for tile in direction: " + dirNo);
            System.out.println("        xTileMin: " + xTileMin + " xTileMax: " + xTileMax);
            System.out.println("        yTileMin: " + yTileMin + " yTileMax: " + yTileMax);
            System.out.println("        zTileMin: " + zTileMin + " zTileMax: " + zTileMax);
            System.out.println("        maxRowTexCoor: " + maxRowTexCoor + " maxColTexCoor: " + maxColTexCoor);
        }
    }

    private void initializeCroppedRange()
    {
        croppedRowMin = rowMin;
        croppedColMin = colMin;
        croppedRowMax = rowMax;
        croppedColMax = colMax;
    }

    public void setAxesFlipped(boolean axesFlipped)
    {
        this.axesFlipped = axesFlipped;
    }

    public void display(StsTextureSurfaceFace textureSurface, float[] planeData, boolean isPixelMode, GL gl, boolean reversePolarity, int n, byte nullByte)
    {
        // convert back to unsigned bytes
        if(planeData == null)
            display(textureSurface, (byte[]) null, isPixelMode, gl, n, nullByte);
        else
        {
            byte[] bData = new byte[planeData.length];
            float val;
            if(reversePolarity)
                for(int i = 0; i < planeData.length; i++)
                {
                    val = planeData[i];
                    val = val > 1.0f ? 1.0f : val;
                    val = val < -1.0f ? -1.0f : val;
                    bData[i] = (byte) ((int) (127 * (-val)) + 127); // signed to unsigned byte, reverse polarity
                }

            else
                for(int i = 0; i < planeData.length; i++)
                {
                    val = planeData[i];
                    val = val > 1.0f ? 1.0f : val;
                    val = val < -1.0f ? -1.0f : val;
                    bData[i] = (byte) ((int) (127 * val) + 127); // signed to unsigned byte
                }
            display(textureSurface, bData, isPixelMode, gl, n, nullByte);
            bData = null;
        }
    }

    public void display(StsTextureSurfaceFace textureSurface, byte[] planeData, boolean isPixelMode, GL gl, int nTile, byte nullByte)
    {
        if(textureSurface == null) return;

        if(texture == 0)
        {
            if(planeData == null) return;
            createTexture(gl);
            if(!bindTexture(gl)) return;
            createTextureBackground(isPixelMode, gl);
        }
        else if(!bindTexture(gl)) return;

        if(planeData != null) addData(gl, planeData);

       // if(debug) printDebug("display", "", nullByte);

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        //gl.glDisable(GL.GL_LIGHTING); nog good -- lighted surfaces // jbw
        if(displayListNum != 0)
            gl.glCallList(displayListNum);
        else
            textureSurface.drawTextureTileSurface(this, gl, true);
        gl.glFlush();
    }

    private void printDebug(String method, String action, byte nullByte)
    {
        String tileDescription = " textureTiles[" + nTile + "] for: " + textureTiles.textureSurface.getName();
        String threadName = " " + Thread.currentThread().getName() + " ";
        String geometryDescription = " texture " + texture + " " + nRows + " rows(" + rowMin + "-" + rowMax + ") " + nCols + " cols(" + colMin + "-" + colMax + ")" + " maxRowTexCoor " + maxRowTexCoor + " maxColTexCoor " + maxColTexCoor;
        String message = tileDescription + threadName + geometryDescription;
        StsCursor3dTexture.textureDebug(tileData, nRows, nCols, this, method, message, nullByte);
    }

    public void display(StsTextureSurfaceFace textureSurface, ByteBuffer planeData, boolean isPixelMode, GL gl, byte nullByte)
    {
        if(textureSurface == null) return;

        if(texture == 0)
        {
            if(planeData == null) return;
            createTexture(gl);
            if(!bindTexture(gl)) return;
            createTextureBackground(isPixelMode, gl);
        }
        else if(!bindTexture(gl)) return;

        if(planeData != null) addData(gl, planeData, nullByte);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

        //gl.glDisable(GL.GL_LIGHTING); nog good -- lighted surfaces // jbw
        if(displayListNum != 0)
            gl.glCallList(displayListNum);
        else
            textureSurface.drawTextureTileSurface(this, gl, true);
        gl.glFlush();
    }

    public void display2d(StsTextureSurfaceFace textureSurface, boolean isPixelMode, byte[] planeData, GL gl, boolean axesFlipped, int nTile, byte nullByte)
    {
        if(textureSurface == null) return;

        if(texture == 0)
        {
            createTexture(gl);
            if(!bindTexture(gl)) return;
            createTextureBackground(isPixelMode, gl);
        }
        else if(!bindTexture(gl)) return;

        if(planeData != null) addData(gl, planeData);

       //  if(debug) printDebug("display2d", "", nullByte);

        this.axesFlipped=axesFlipped;

        gl.glColor4f(1.0f,1.0f,1.0f,1.0f);

    //		System.out.println("TextureTile lighting: " + gl.glIsEnabled(GL.GL_LIGHTING));
    //		gl.glDisable(GL.GL_LIGHTING); // jbw

        if(displayListNum2d!=0)
            gl.glCallList(displayListNum2d);
        else
            textureSurface.drawTextureTileSurface(this, gl, false);
        gl.glFlush();
    }

    public void display2d(StsTextureSurfaceFace textureSurface, boolean isPixelMode, ByteBuffer planeData, GL gl, boolean axesFlipped, int n, byte nullByte)
    {
        if(textureSurface == null) return;

        if(texture == 0)
        {
            if(planeData == null) return;
            createTexture(gl);
            if(!bindTexture(gl)) return;
            createTextureBackground(isPixelMode, gl);
        }
        else if(!bindTexture(gl)) return;

        if(planeData != null) addData(gl, planeData, nullByte);

        // if(debug) printDebug("display2d", "", nullByte);

        this.axesFlipped = axesFlipped;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

        gl.glDisable(GL.GL_LIGHTING); // jbw

        if(displayListNum2d != 0)
            gl.glCallList(displayListNum2d);
        else
            textureSurface.drawTextureTileSurface(this, gl, false);
        gl.glFlush();
    }

    /** A default method which draw this texture in 3d on a rectangle. */
    public void drawQuadSurface3d(GL gl, byte zDomainData)
    {
        int dir = textureTiles.dir;
        if(debug) textureTileDebug(this, "drawQuadSurface3d", dir, zDomainData);
        gl.glBegin(GL.GL_QUADS);
        gl.glTexCoord2d(minColTexCoor, minRowTexCoor);
        gl.glVertex3dv(xyzPlane[0], 0);
        gl.glTexCoord2d(maxColTexCoor, minRowTexCoor);
        gl.glVertex3dv(xyzPlane[1], 0);
        gl.glTexCoord2d(maxColTexCoor, maxRowTexCoor);
        gl.glVertex3dv(xyzPlane[2], 0);
        gl.glTexCoord2d(minColTexCoor, maxRowTexCoor);
        gl.glVertex3dv(xyzPlane[3], 0);
        gl.glEnd();
    }
    
    public String getTileDrawString()
    {
        String textureCoordinatesString = " rowTexture: " + minRowTexCoor + " - " + maxRowTexCoor + " colTexture: " + minRowTexCoor + " - " + maxRowTexCoor;
        String textureGeometryString = " dir " + textureTiles.dir + " quad min: " + StsMath.toString(xyzPlane[0]) + " max: " +  StsMath.toString(xyzPlane[0]);
        return textureCoordinatesString + textureGeometryString;
    }

    /** A default method which draw this texture in 2d on a rectangle. */
    public void drawQuadSurface2d(GL gl)
    {
        int dir = textureTiles.dir;
        if(debug) textureTileDebug(this, "drawQuadSurface2d", getTileDrawString());
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glDisable(GL.GL_LIGHTING); // jbw
 
        if(dir == StsCursor3d.XDIR)
        {
            gl.glBegin(GL.GL_QUADS);
            gl.glTexCoord2d(minColTexCoor, minRowTexCoor);
            gl.glVertex2d(xyzPlane[0][1], xyzPlane[0][2]);
            gl.glTexCoord2d(maxColTexCoor, minRowTexCoor);
            gl.glVertex2d(xyzPlane[1][1], xyzPlane[1][2]);
            gl.glTexCoord2d(maxColTexCoor, maxRowTexCoor);
            gl.glVertex2d(xyzPlane[2][1], xyzPlane[2][2]);
            gl.glTexCoord2d(minColTexCoor, maxRowTexCoor);
            gl.glVertex2d(xyzPlane[3][1], xyzPlane[3][2]);
            gl.glEnd();
        }
        else if(dir == StsCursor3d.YDIR)
        {
            gl.glBegin(GL.GL_QUADS);
            gl.glTexCoord2d(minColTexCoor, minRowTexCoor);
            gl.glVertex2d(xyzPlane[0][0], xyzPlane[0][2]);
            gl.glTexCoord2d(maxColTexCoor, minRowTexCoor);
            gl.glVertex2d(xyzPlane[1][0], xyzPlane[1][2]);
            gl.glTexCoord2d(maxColTexCoor, maxRowTexCoor);
            gl.glVertex2d(xyzPlane[2][0], xyzPlane[2][2]);
            gl.glTexCoord2d(minColTexCoor, maxRowTexCoor);
            gl.glVertex2d(xyzPlane[3][0], xyzPlane[3][2]);
            gl.glEnd();
        }
        else if(dir == StsCursor3d.ZDIR)
        {
            if(!axesFlipped)
            {
                gl.glBegin(GL.GL_QUADS);
                gl.glTexCoord2d(minColTexCoor, minRowTexCoor);
                gl.glVertex2dv(xyzPlane[0], 0);
                gl.glTexCoord2d(maxColTexCoor, minRowTexCoor);
                gl.glVertex2dv(xyzPlane[1], 0);
                gl.glTexCoord2d(maxColTexCoor, maxRowTexCoor);
                gl.glVertex2dv(xyzPlane[2], 0);
                gl.glTexCoord2d(minColTexCoor, maxRowTexCoor);
                gl.glVertex2dv(xyzPlane[3], 0);
                gl.glEnd();
            }
            else
            {
                gl.glBegin(GL.GL_QUADS);
                gl.glTexCoord2d(minColTexCoor, minRowTexCoor);
                gl.glVertex2d(xyzPlane[0][1], xyzPlane[0][0]);
                gl.glTexCoord2d(maxColTexCoor, minRowTexCoor);
                gl.glVertex2d(xyzPlane[1][1], xyzPlane[1][0]);
                gl.glTexCoord2d(maxColTexCoor, maxRowTexCoor);
                gl.glVertex2d(xyzPlane[2][1], xyzPlane[2][0]);
                gl.glTexCoord2d(minColTexCoor, maxRowTexCoor);
                gl.glVertex2d(xyzPlane[3][1], xyzPlane[3][0]);
                gl.glEnd();
            }
        }
    }

    /** A default method which draw this texture in 2d on a series of rectangular coordinates. */
    public void drawQuadStripCursorSurface2d(GL gl, double[] zCoordinates)
    {
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glDisable(GL.GL_LIGHTING); // jbw
        int dir = textureTiles.dir;
        if(dir == StsCursor3d.XDIR)
            drawQuadStripSurface2d(gl, xyzPlane[0][1], xyzPlane[3][1], zCoordinates);
        else if(dir == StsCursor3d.YDIR)
            drawQuadStripSurface2d(gl, xyzPlane[0][0], xyzPlane[3][0], zCoordinates);
    }

    public void drawQuadStripSurface2d(GL gl, double xMin, double xMax, double[] yCoordinates)
    {
        if(debug) StsException.systemDebug(this, "drawQuadStripSurface2d(gl, xMin, xMax, yCoordinates)", "draw called.");
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glBegin(GL.GL_QUAD_STRIP);
        double colTexCoor = minColTexCoor;
        for(int col = croppedColMin; col <= croppedColMax; col++, colTexCoor += dColTexCoor)
        {
            gl.glTexCoord2d(colTexCoor, minRowTexCoor);
            gl.glVertex2d(xMin, yCoordinates[col]);
            gl.glTexCoord2d(colTexCoor, maxRowTexCoor);
            gl.glVertex2d(xMax, yCoordinates[col]);
        }
        gl.glEnd();
    }

    public void drawQuadStripSurface2d(GL gl, double xMin, double xInc, double[][] yCoordinates)
    {
        if(yCoordinates == null) return;
        if(debug) StsException.systemDebug(this, "drawQuadStripSurface2d(gl, arcCoordinates, yCoordinates)", "draw called.");
         gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		gl.glDisable(GL.GL_LIGHTING);
        double rowTexCoor = minRowTexCoor;
        double x2 = xMin + rowMin * xInc;
        for(int row = croppedRowMin; row < croppedRowMax; row++, rowTexCoor += dRowTexCoor)
        {
            gl.glBegin(GL.GL_QUAD_STRIP);
            double colTexCoor = minColTexCoor;
            double x1 = x2;
            x2 += xInc;
            for(int col = croppedColMin; col <= croppedColMax; col++, colTexCoor += dColTexCoor)
            {
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex2d(x1, yCoordinates[row][col]);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex2d(x2, yCoordinates[row][col]);
            }
            gl.glEnd();
        }
    }

    public void drawQuadStripSurface2d(GL gl, double[] xCoordinates, double[] yCoordinates)
    {
        if(yCoordinates == null) return;
        if(debug) System.out.println("drawQuadStripSurface2d(gl, arcCoordinates, yCoordinates) called.");
         gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		gl.glDisable(GL.GL_LIGHTING);
        double rowTexCoor = minRowTexCoor;
        double x2 = xCoordinates[croppedRowMin];
        for(int row = croppedRowMin; row < croppedRowMax; row++, rowTexCoor += dRowTexCoor)
        {
            gl.glBegin(GL.GL_QUAD_STRIP);
            double colTexCoor = minColTexCoor;
            double x1 = x2;
            x2 = xCoordinates[row + 1];
            for(int col = croppedColMin; col <= croppedColMax; col++, colTexCoor += dColTexCoor)
            {
                double y = yCoordinates[col];
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex2d(x1, y);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex2d(x2, y);
            }
            gl.glEnd();
        }
    }

    private void createTexture(GL gl)
    {
        int[] textures = new int[1];
        gl.glGenTextures(1, textures, 0);
        texture = textures[0];
        if(StsTextureTiles.debug) textureTileDebug(this, "createTexture", "create", nTile);
    }

    public void textureTileDebug(Object object, String methodName, int dir, byte zDomainData)
    {
        textureTileDebug(object, methodName, getAxisRangeString(zDomainData, dir));
    }

    public void textureTileDebug(Object object, String methodName, String action, int nTile)
    {
        textureTileDebug(object, methodName, action);
    }

    public void textureTileDebug(Object object, String method, String action)
    {
        StsException.systemDebug(object, method, action + " textureTiles[" + nTile + "] texture " + texture);
    }

    public void textureTileDebug(Object object, String method, String action, byte zDomainData, byte zDomainProject, int dir, byte nullByte)
    {
        if(dir == ZDIR || zDomainData != StsProject.TD_DEPTH || zDomainProject != StsProject.TD_TIME) return;
        String rangeString = getAxisRangeString(zDomainData, dir);
        StsException.systemDebug(object, method, action + " textureTiles[" + nTile + "] texture " + texture + " range: " + rangeString);
        // debug check DEPTH to TIME domain cursor coordinates
        int[] rowColSample = StsCursor3dTexture.textureDebugGetFirstNonNullRowColSample(tileData, nRows, nCols, nullByte);
        if(rowColSample == null) return;
        float[] xyz = getTextureXYZCoordinates(rowColSample[0], rowColSample[1], dir);
        float time = StsCursor3dTexture.getTimeFromIntervalVelocities(xyz);
        StsException.systemDebug("    DEPTH->TIME debug. First nonNull sample " + rowColSample[2] + " is at x " + xyz[0] + " y " + xyz[1] + " depth " + xyz[2] + " time " + time);
    }

    private float[] getTextureXYZCoordinates(int row, int col, int dir)
    {
        StsRotatedGridBoundingBox boundingBox = textureTiles.boundingBox;
        float xInc, yInc;
        if(dir == XDIR)
        {
            xInc = 0;
            yInc = boundingBox.yInc;
        }
        else // dirNo == StsCursor3d.YDIR
        {
            xInc = boundingBox.xInc;
            yInc = 0;
        }
        float zInc = boundingBox.zInc;
        double[] xyzMin = xyzPlane[0];
        float x0 = (float) xyzMin[0];
        float y0 = (float) xyzMin[1];
        float z0 = (float) xyzMin[2];
        row += croppedRowMin;
        float x = x0 + row*xInc;
        float y = y0 + row*yInc;
        float z = z0 + col*zInc;
        return new float[] { x, y, z };
    }

    private String getAxisRangeString(byte zDomainData, int dir)
    {
        double x, y, z;
        double xMin, xMax, yMin, yMax, zMin, zMax;
        String timeOrDepthLabel = (zDomainData == 0 ? " depth " : " time " );
        if(dir == XDIR)
        {
            x = xyzPlane[0][0];
            yMin = xyzPlane[0][1];
            yMax = xyzPlane[3][1];
            zMin = xyzPlane[0][2];
            zMax = xyzPlane[1][2];
            return " XDIR x " +  x + " y " + yMin + " - " + yMax + timeOrDepthLabel + zMin + " - " + zMax;
        }
        else if(dir == YDIR)
        {
            y = xyzPlane[0][1];
            xMin = xyzPlane[0][0];
            xMax = xyzPlane[3][0];
            zMin = xyzPlane[0][2];
            zMax = xyzPlane[1][2];
            return " YDIR y " +  y + " x " + xMin + " - " + xMax + timeOrDepthLabel + zMin + " - " + zMax;
        }
        else
        {
            z = xyzPlane[0][2];
            xMin = xyzPlane[0][0];
            xMax = xyzPlane[1][0];
            yMin = xyzPlane[0][1];
            yMax = xyzPlane[3][1];
            return " ZDIR z " +  z + " x " + xMin + " - " + xMax + " y " + yMin + " - " + yMax;
        }
    }

    private boolean bindTexture(GL gl)
    {

        if(texture == 0)
        {
            StsException.systemError("Attempt to bind a 0 texture");
            return false;
        }

        gl.glActiveTexture(GL.GL_TEXTURE0);
        gl.glBindTexture(GL.GL_TEXTURE_2D, texture);
        return true;
    }

    private void createTextureBackground(boolean isPixelMode, GL gl)
    {
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        // ByteBuffer background = BufferUtil.newByteBuffer(nBackgroundRows*nBackgroundCols);
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE);
        if(!isPixelMode)
        {
            gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR);
            gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR);
        }
        else
        {
            gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
            gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);
        }
        /*
            int errorCode = glu.gluBuild2DMipmaps( GL.GL_TEXTURE_2D, GL.GL_RGBA8, nBackgroundCols, nBackgroundRows, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, background);
            if(errorCode != 0)
            {
             StsException.systemError("StsTextureTile.createTextureBackground() failed. GLU error code: " + glu.errorString(errorCode));
            }
           */
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1); // jbw essential to set june 27
        if((textureTiles.shader != StsJOGLShader.NONE) && StsJOGLShader.canDoARBShader(gl))
        {
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);
            gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_INTENSITY8 /*GL.GL_COMPRESSED_INTENSITY*/, nBackgroundCols, nBackgroundRows, 0, GL.GL_LUMINANCE, GL.GL_UNSIGNED_BYTE, null); //background
            //System.out.println("bg image2D "+gl.glGetError());
            //StsJOGLShader.reloadTLUT(gl);
            //StsJOGLShader.enableARBShader(gl, textureTiles.shader);
        }
        else
            gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, nBackgroundCols, nBackgroundRows, 0, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, null); //background

        //background = null;
    }

    protected void addData(GL gl, byte[] planeData)
    {
        int nData, nTile;
        if(planeData != null)
            try
            {
                if(planeData != null)
                {
                    if(runTimer) timer.start();

                    if(nTotalCols == nCols)
                    {
                        nData = rowMin * nTotalCols + colMin;
                        nTile = 0;
                        tileData = BufferUtil.newByteBuffer(nRows * nCols);
                        tileData.put(planeData, nData, nCols * nRows);
                        //                  System.arraycopy(planeData, nData, tileData, nTile, nTextureCols*nTextureRows);
                    }
                    else
                    {
                        nData = rowMin * nTotalCols + colMin;
                        nTile = 0;
                        tileData = BufferUtil.newByteBuffer(nRows * nCols);
                        int nDataRows = planeData.length / nTotalCols;
                        int nReadRows = Math.min(nDataRows, nRows);
                        for(int row = 0; row < nReadRows; row++)
                        {
                            tileData.put(planeData, nData, nCols);
                            nData += nTotalCols;
                            nTile += nCols;
                        }
                    }
                    if(runTimer) timer.stopPrint("        add data for tile " + texture);
                }
                tileData.rewind();
                if(runTimer) timer.start();
                gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1); // jbw essential to set june 27
                if((textureTiles.shader != StsJOGLShader.NONE) && StsJOGLShader.canDoARBShader(gl))
                {
                    gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);
                    gl.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, nCols, nRows, GL.GL_LUMINANCE, GL.GL_UNSIGNED_BYTE, tileData);
                }
                else
                    gl.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, nCols, nRows, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, tileData);
                if(!debug) tileData = null;
                if(runTimer) timer.stopPrint("        add subImage to texture");
            }
            catch(Exception e)
            {
                StsException.outputException("StsTextureTile.addData() failed.", e, StsException.WARNING);
            }
    }

    protected void addData(GL gl, ByteBuffer planeData, byte nullByte)
    {
        int nData;
        if(planeData != null)
            try
            {
                if(planeData != null)
                {
                    if(runTimer) timer.start();

                    if(nTotalCols == nCols)
                    {
                        nData = rowMin * nTotalCols + colMin;
                        tileData = planeData;
                        tileData.position(nData);
                    }
                    else
                    {
                        nData = rowMin * nTotalCols + colMin;
                        tileData = ByteBuffer.allocateDirect(nRows * nCols);
                        int nReadRows = Math.min(nTotalRows, nRows);
                        byte[] rowData = new byte[nCols];
                        for(int row = 0; row < nReadRows; row++)
                        {
                            planeData.position(nData);
                            planeData.get(rowData);
                            tileData.put(rowData);
                            nData += nTotalCols;
                        }
                        tileData.rewind();
                    }
                    if(debug) printDebug("addData", "adding new data", nullByte);

                    if(runTimer) timer.stopPrint("        add data for tile " + texture);
                }
                if(runTimer) timer.start();
                gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1); // jbw essential to set june 27
                if((textureTiles.shader == StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS) && StsJOGLShader.canUseShader)
                {
                    gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);
                    gl.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, nCols, nRows, GL.GL_LUMINANCE, GL.GL_UNSIGNED_BYTE, tileData);
                }
                else
                    gl.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, nCols, nRows, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, tileData);

                if(!debug) tileData = null;
                if(runTimer) timer.stopPrint("        add subImage to texture");
            }
            catch(Exception e)
            {
                StsException.outputException("StsTextureTile.addData() failed.", e, StsException.WARNING);
            }
    }

    public boolean deleteTexture(GL gl, int n)
    {
        if(texture == 0) return false;
        if(StsTextureTiles.debug) textureTileDebug(this, "deleteTexture", "deleteTileTexture");
        gl.glDeleteTextures(1, new int[]{texture}, 0);
        texture = 0;
        return true;
    }

    public boolean deleteDisplayList(GL gl, int n)
    {
        boolean deleted = false;
        if(displayListNum != 0)
        {
            if(StsTextureTiles.debug)
                textureTileDebug(this, "deleteDisplayList", "delete tile 3d displayList");
            gl.glDeleteLists(displayListNum, 1);
            displayListNum = 0;
            deleted = true;
        }
        if(displayListNum2d != 0)
        {
            if(StsTextureTiles.debug)
                textureTileDebug(this, "deleteDisplayList", "delete tile 2d displayList");
            gl.glDeleteLists(displayListNum2d, 1);
            displayListNum2d = 0;
            deleted = true;
        }
        return deleted;
    }

    public void setDirCoordinate(float dirCoordinate)
    {
        int dir = textureTiles.dir;
        for(int n = 0; n < 4; n++)
            xyzPlane[n][dir] = dirCoordinate;
    }

    public void adjust(int cropRowMin, int cropRowMax, int cropColMin, int cropColMax)
    {
        double rowBorder = 1.0 / (2.0 * nBackgroundRows);
        double colBorder = 1.0 / (2.0 * nBackgroundCols);

        croppedRowMin = StsMath.minMax(cropRowMin, rowMin, rowMax);
        croppedRowMax = StsMath.minMax(cropRowMax, rowMin, rowMax);
        croppedColMin = StsMath.minMax(cropColMin, colMin, colMax);
        croppedColMax = StsMath.minMax(cropColMax, colMin, colMax);

        minRowTexCoor = ((double) croppedRowMin - (double) rowMin) / nBackgroundRows + rowBorder;
        minColTexCoor = ((double) croppedColMin - (double) colMin) / nBackgroundCols + colBorder;
        maxRowTexCoor = ((double) croppedRowMax - (double) rowMin) / nBackgroundRows + rowBorder;
        maxColTexCoor = ((double) croppedColMax - (double) colMin) / nBackgroundCols + colBorder;

        double minRowF = (double) croppedRowMin / (nTotalRows - 1);
        double minColF = (double) croppedColMin / (nTotalCols - 1);
        double maxRowF = (double) croppedRowMax / (nTotalRows - 1);
        double maxColF = (double) croppedColMax / (nTotalCols - 1);

        double xTileMin = 0.0f, xTileMax = 0.0f;
        double yTileMin = 0.0f, yTileMax = 0.0f;
        double zTileMin = 0.0f, zTileMax = 0.0f;
        int dir = textureTiles.dir;
        StsRotatedGridBoundingBox boundingBox = textureTiles.boundingBox;
        float dirCoordinate = textureTiles.dirCoordinate;
        switch(dir)
        {
            case StsCursor3d.XDIR:
                yTileMin = boundingBox.yMin + minRowF * (boundingBox.yMax - boundingBox.yMin);
                yTileMax = boundingBox.yMin + maxRowF * (boundingBox.yMax - boundingBox.yMin);
                zTileMin = boundingBox.zMin + minColF * (boundingBox.zMax - boundingBox.zMin);
                zTileMax = boundingBox.zMin + maxColF * (boundingBox.zMax - boundingBox.zMin);
                xyzPlane[0] = new double[] {dirCoordinate, yTileMin, zTileMin};
                xyzPlane[1] = new double[] {dirCoordinate, yTileMin, zTileMax};
                xyzPlane[2] = new double[] {dirCoordinate, yTileMax, zTileMax};
                xyzPlane[3] = new double[] {dirCoordinate, yTileMax, zTileMin};
                break;
            case StsCursor3d.YDIR:
                xTileMin = boundingBox.xMin + minRowF * (boundingBox.xMax - boundingBox.xMin);
                xTileMax = boundingBox.xMin + maxRowF * (boundingBox.xMax - boundingBox.xMin);
                zTileMin = boundingBox.zMin + minColF * (boundingBox.zMax - boundingBox.zMin);
                zTileMax = boundingBox.zMin + maxColF * (boundingBox.zMax - boundingBox.zMin);
                xyzPlane[0] = new double[] {xTileMin, dirCoordinate, zTileMin};
                xyzPlane[1] = new double[] {xTileMin, dirCoordinate, zTileMax};
                xyzPlane[2] = new double[] {xTileMax, dirCoordinate, zTileMax};
                xyzPlane[3] = new double[] {xTileMax, dirCoordinate, zTileMin};
                break;
            case StsCursor3d.ZDIR:
                xTileMin = boundingBox.xMin + minColF * (boundingBox.xMax - boundingBox.xMin);
                xTileMax = boundingBox.xMin + maxColF * (boundingBox.xMax - boundingBox.xMin);
                yTileMin = boundingBox.yMin + minRowF * (boundingBox.yMax - boundingBox.yMin);
                yTileMax = boundingBox.yMin + maxRowF * (boundingBox.yMax - boundingBox.yMin);
                xyzPlane[0] = new double[] {xTileMin, yTileMin, dirCoordinate};
                xyzPlane[1] = new double[] {xTileMax, yTileMin, dirCoordinate};
                xyzPlane[2] = new double[] {xTileMax, yTileMax, dirCoordinate};
                xyzPlane[3] = new double[] {xTileMin, yTileMax, dirCoordinate};
                break;
        }

        if(debug)
        {
            System.out.println("StsTextureTile.adjust() called.");
            System.out.println("    Geometry for tile in direction: " + dir);
            System.out.println("        xTileMin: " + xTileMin + " xTileMax: " + xTileMax);
            System.out.println("        yTileMin: " + yTileMin + " yTileMax: " + yTileMax);
            System.out.println("        zTileMin: " + zTileMin + " zTileMax: " + zTileMax);
            System.out.println("        maxRowTexCoor: " + maxRowTexCoor + " maxColTexCoor: " + maxColTexCoor);
        }
    }

    public void constructSurface(StsTextureSurfaceFace surface, GL gl, boolean useDisplayLists, boolean is3d, int nTile)
    {
        deleteDisplayList(gl, nTile);
        if(!useDisplayLists) return;
        if(is3d)
        {
            displayListNum = gl.glGenLists(1);
            gl.glNewList(displayListNum, GL.GL_COMPILE);
            if(StsTextureTiles.debug)
                textureTileDebug(this, "constructSurface", "construct 3d displayList");
        }
        else
        {
            displayListNum2d = gl.glGenLists(1);
            gl.glNewList(displayListNum2d, GL.GL_COMPILE);
            if(StsTextureTiles.debug)
                textureTileDebug(this, "constructSurface", "construct 2d displayList");
        }
        surface.drawTextureTileSurface(this, gl, is3d);
        gl.glEndList();
    }
    /*
      public void draw(GL gl, boolean useDisplayLists)
      {
       if(displayListNum != 0) gl.deleteLists(displayListNum, 1);
       displayListNum = 0;
       if(!useDisplayLists) return;
       displayListNum = gl.genLists(1);
       gl.glNewList(displayListNum, GL.GL_COMPILE);
       surface.drawTextureTileSurface(this, gl);
       gl.glEndList();
      }
      */
}

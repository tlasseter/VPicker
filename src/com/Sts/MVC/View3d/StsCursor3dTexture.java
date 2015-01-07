package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.nio.ByteBuffer;

/** This class is used directly by all volume objects displayed as textures on cursor sections.
 *  It is sublclassed by StsCursor3dVolumeTexture for simple IJK-orthogonal volumes (seismic, stimulated volumes).
 */
abstract public class StsCursor3dTexture implements StsTextureSurfaceFace
{
    // These three members are initialized by initializeTransients call from viewPersistManager
    transient public StsCursor3d cursor3d;
    transient public int dirNo;
    transient public float dirCoordinate = -StsParameters.largeFloat;
    transient public StsRotatedGridBoundingBox cursorBoundingBox;

    /** convenience copy of the current model */
    transient protected StsModel model;
    transient public int nTextureRows, nTextureCols;
    transient boolean isVisible = true;
    /** Tiles on which texture is generated */
    transient public StsTextureTiles textureTiles;
    /** texture has changed because 1) cursor has moved; 2) texture on cursor has changed; 3) subVolume changed; pixelMode changed.  Run same displayList with new texture. */
    transient boolean textureChanged = true;
    /** geometry has changed because 1) domain changed (time or depth); 2) size of cursor has changed.  Rebuild displayLists. Use same texture unless cursor has moved. */
    transient boolean geometryChanged = true;
    transient public byte nullByte = StsParameters.nullByte;

    static public final int XDIR = StsCursor3d.XDIR;
    static public final int YDIR = StsCursor3d.YDIR;
    static public final int ZDIR = StsCursor3d.ZDIR;

    static public final byte TD_TIME = StsProject.TD_TIME;
    static public final byte TD_DEPTH = StsProject.TD_DEPTH;

    static public final boolean debug = false;

    abstract public boolean setDirCoordinate(float dirCoordinate);

    //    public void clearTextureDisplay();
    abstract public void cropChanged();

    abstract public boolean isDisplayableObject(Object object);

    abstract public boolean isDisplayingObject(Object object);

    abstract public void displayTexture(StsGLPanel3d glPanel3d, boolean is3d, StsCursorSection cursorSection);

    abstract public void display(StsGLPanel3d glPanel3d, boolean is3d);

    abstract public String propertyReadout(StsPoint point);

    abstract public boolean setObject(Object object);

    abstract public Object getObject();

    abstract public Class getDisplayableClass();

    abstract public boolean canDisplayClass(Class c);

    abstract public void subVolumeChanged();

    public boolean stsClassIsCursor3dTextureDisplayable(StsClass stsClass)
    {
        return canDisplayClass(stsClass.getClass());
    }

    public boolean initialize(StsModel model, StsCursor3d cursor3d, int dir)
    {
        this.model = model;
        this.dirNo = dir;
        this.cursor3d = cursor3d;
        StsCursorSection cursorSection = cursor3d.cursorSections[dir];
        this.dirCoordinate = cursorSection.dirCoordinate;
        this.cursorBoundingBox = cursor3d.rotatedBoundingBox;
        return true;
    }

    protected boolean isInRange(int row, int col)
    {
        if(row < 0 || row >= nTextureRows) return false;
        if(col < 0 || col >= nTextureCols) return false;
        return true;
    }

    protected boolean isInRange(float rowF, float colF)
    {
        if(rowF < 0 || rowF >= nTextureRows) return false;
        if(colF < 0 || colF >= nTextureCols) return false;
        return true;
    }

    public boolean isVisible()
    {
        return isVisible;
    }

    public void clearTextureDisplay()
    {
        textureChanged();
        geometryChanged();
    }

    public void geometryChanged()
    {
        geometryChanged = true;
        if(debug) StsException.systemDebug(this, "geometryChanged");
    }

    /**
     * This puts texture display on delete list.  Operation is performed
     * at beginning of next draw operation.
     */
/*
    public void addTextureToDeleteList()
    {
        if(textureTiles != null)
        {
            if(debug) System.out.println(debugMessage("clearTextureTileSurface() deleting textureTileSurface"));
            StsTextureList.addTextureToDeleteList(this);
            textureChanged = true;
        }
        if(lineSet3d == null) return;
        if((dirNo != StsCursor3d.ZDIR) || (lineSet3d.getDisplayAttribute().getName() == lineSet3d.ATTRIBUTE_NONE_STRING))
        {
            if(debug) System.out.println(debugMessage("cleareTextureTileSurface() display on X or Y section changed."));
//			planeData = null;
            textureChanged = true;
        }
    }
*/
    protected String debugMessage(String message)
    {
        return new String(StsToolkit.getSimpleClassname(this) + "[" + dirNo + "] " + message);
    }

    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        deleteTextures(gl);
        deleteDisplayLists(gl);
        if(debug) StsException.systemDebug(this, "deleteTexturesAndDisplayLists", debugMessage("deleted textures"));
    }

    public void deleteTextures(GL gl)
    {
        if(textureTiles == null) return;
        textureTiles.deleteTextures(gl);
        if(debug) StsException.systemDebug(this, "deleteTextures", debugMessage("deleted textures"));
        textureChanged = true;
    }

    protected void deleteDisplayLists(GL gl)
    {
        if(textureTiles != null) textureTiles.deleteDisplayLists(gl);
        if(debug) StsException.systemDebug(this, "deleteDisplayLists", debugMessage("deleted texture displayLists"));
    }

    public void disableShader(GL gl)
    {
        if(textureTiles == null) return;
        if(textureTiles.shader != StsJOGLShader.NONE)
            StsJOGLShader.disableARBShader(gl);
    }

    public void drawTextureTileSurface(StsRotatedGridBoundingBox volume, StsTextureTile tile, GL gl, boolean is3d, byte zDomainData)
    {
        byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
        int nTile = tile.nTile;
        if(projectZDomain == StsParameters.TD_TIME)
        {
            if(zDomainData == StsParameters.TD_TIME)
                drawTextureQuadTileSurface(tile, gl, is3d, zDomainData);
            else
            {
                StsSeismicVelocityModel velocityModel = model.getProject().velocityModel;
                drawTextureTileDepthSurfaceInTime(volume, velocityModel, tile, gl);
            }
        }
        else if(projectZDomain == StsParameters.TD_DEPTH)
        {
            StsModel model = StsObject.getCurrentModel();
            if(zDomainData == StsParameters.TD_DEPTH) // seismic already in depth, don't need to convert so draw as if in time
                drawTextureQuadTileSurface(tile, gl, is3d, zDomainData);
            else
            {
                StsSeismicVelocityModel velocityModel = model.getProject().velocityModel;
                if(velocityModel == null) return;
                drawTextureTileTimeSurfaceInDepth(volume, velocityModel, tile, gl);
            }
        }
    }

    public void drawTextureQuadTileSurface(StsTextureTile tile, GL gl, boolean is3d, byte zDomainData)
    {
        if(is3d)
            tile.drawQuadSurface3d(gl, zDomainData);
        else
            tile.drawQuadSurface2d(gl);
    }

    public void drawTextureTileTimeSurfaceInDepth(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        int dir = tile.textureTiles.dir;
        if(dir == StsCursor3d.ZDIR) return;
        if(velocityModel == null) return;

        if(velocityModel.inputVelocityVolume != null)
            drawTextureTileTimeSurfaceInDepthWithVolume(volume, velocityModel, tile, gl, dir);
        else
            drawTextureTileTimeSurfaceInDepthWithIntervalVelocities(volume, velocityModel, tile, gl, dir);
    }

    public void drawTextureTileTimeSurfaceInDepthWithVolume(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
    {
        int row = -1, col = -1;
        try
        {
            if(debug) tile.textureTileDebug(this, "drawTextureTileTimeSurfaceInDepthWithVolume", "", TD_TIME, TD_DEPTH, dir, nullByte);
            StsSeismicVolume velocityVolume = velocityModel.velocityVolume;
            gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            double rowTexCoor = tile.minRowTexCoor;
            double dRowTexCoor = tile.dRowTexCoor;
            double dColTexCoor = tile.dColTexCoor;
            double[] xyzMin = tile.xyzPlane[0];
            float displayXMin = (float) xyzMin[0];
            float displayYMin = (float) xyzMin[1];
            float displayTMin = (float) xyzMin[2];
            double[] xyzMax = tile.xyzPlane[2];

            float depthDatum = velocityModel.depthDatum;
            float timeDatum = velocityModel.timeDatum;
            int[] rowCol = velocityVolume.getCursorDataNRowCols(dir);
            //		int nDataRows = rowCol[0];
            int nDataCols = rowCol[1];
            //		int volumeRowInc = 0;
            //		int volumeColInc = 0;
            byte[] velocityBytes;
            int velocityRow, velocityCol;
            float cursorXInc, cursorYInc;
            float velocityTInc = velocityVolume.zInc;
            if(dir == XDIR)
            {
                velocityCol = velocityVolume.getNearestBoundedColCoor(displayXMin);
                velocityBytes = velocityVolume.readBytePlaneData(dir, velocityCol);
                cursorXInc = 0;
                cursorYInc = volume.getYInc();
                //			volumeRowInc = 1;
            }
            else // dirNo == StsCursor3d.YDIR
            {
                velocityRow = velocityVolume.getNearestBoundedRowCoor(displayYMin);
                velocityBytes = velocityVolume.readBytePlaneData(dir, velocityRow);
                cursorXInc = volume.getXInc();
                cursorYInc = 0;
                //			volumeColInc = 1;
            }
            double tInc = volume.getZInc();
            float x1 = displayXMin;
            float y1 = displayYMin;
            float[] velocityTrace1 = getVelocityTrace(velocityBytes, dir, x1, y1, velocityVolume);
            for(row = tile.croppedRowMin; row < tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
            {
                float x0 = x1;
                float y0 = y1;
                float[] velocityTrace0 = velocityTrace1;
                x1 += cursorXInc;
                y1 += cursorYInc;
                velocityTrace1 = getVelocityTrace(velocityBytes, dir, x1, y1, velocityVolume);
                int i = row * nDataCols + tile.colMin;
                gl.glBegin(GL.GL_QUAD_STRIP);
                double colTexCoor = tile.minColTexCoor;
                float t = displayTMin;

                for(col = tile.croppedColMin; col <= tile.croppedColMax; col++, t += tInc, colTexCoor += dColTexCoor, i++)
                {
                    float v0 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace0);
                    float z0 = (v0 * (t - timeDatum) + depthDatum);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor);
                    gl.glVertex3d(x0, y0, z0);

                    float v1 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace1);
                    float z1 = (v1 * (t - timeDatum) + depthDatum);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                    gl.glVertex3d(x1, y1, z1);
                }
                gl.glEnd();
                //			volumeRow += volumeRowInc;
                //			volumeCol += volumeColInc;
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawTextureTileDepthSurfaceInTime", " at row " + row + " col " + col, e);
        }
    }

    float[] getVelocityTrace(byte[] velocityData, int dir, float x, float y, StsSeismicVolume velocityVolume)
    {
        int nSlices = velocityVolume.nSlices;
        float[] velocities = new float[nSlices];
        if(dir == XDIR)
        {
            int row = velocityVolume.getNearestBoundedRowCoor(y);
            int i = row * nSlices;
            for(int k = 0; k < nSlices; k++, i++)
                velocities[k] = velocityVolume.computeFloatFromByte(velocityData[i]);
        }
        else // dirNo == XDIR
        {
            int col = velocityVolume.getNearestBoundedColCoor(x);
            int i = col * nSlices;
            for(int k = 0; k < nSlices; k++, i++)
                velocities[k] = velocityVolume.computeFloatFromByte(velocityData[i]);
        }
        // debug check
        /*
            float v1 = 0.0f;
            for(int k = 0; k < nCroppedSlices; k++)
            {
                float v0 = v1;
                v1 = velocities[k];
                if(v1 < v0)
                    System.out.println("Trace is not monotonically increasing at x: " + x + " y: " + y);
            }
         */
        return velocities;
    }

    private float getVelocityFromTrace(float[] velocities, double t, double timeDatum, float velocityTInc)
    {
        int index = (int) Math.round((t - timeDatum) / velocityTInc);
        if(index < 0)
            return velocities[0];
        int length = velocities.length;
        if(index >= length)
            return velocities[length - 1];
        else
            return velocities[index];
    }

    public void drawTextureTileTimeSurfaceInDepthWithIntervalVelocities(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
    {
        int n = 0, row = -1, col = -1;
        try
        {
            if(debug) tile.textureTileDebug(this, "drawTextureTileTimeSurfaceInDepthWithVolume", "", TD_TIME, TD_DEPTH, dir, nullByte);
            double rowTexCoor = tile.minRowTexCoor;
            double dRowTexCoor = tile.dRowTexCoor;
            double dColTexCoor = tile.dColTexCoor;
            double[] xyzMin = tile.xyzPlane[0];
            float displayXMin = (float) xyzMin[0];
            float displayYMin = (float) xyzMin[1];
            float displayTMin = (float) xyzMin[2];
            double[] xyzMax = tile.xyzPlane[2];

            float cursorXInc, cursorYInc;
            if(dir == XDIR)
            {
                cursorXInc = 0;
                cursorYInc = volume.getYInc();
            }
            else // dirNo == StsCursor3d.YDIR
            {
                cursorXInc = volume.getXInc();
                cursorYInc = 0;
            }
            float tInc = volume.getZInc();
            float x1 = displayXMin;
            float y1 = displayYMin;
            gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            for(row = tile.croppedRowMin; row < tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
            {
                float x0 = x1;
                float y0 = y1;
                x1 += cursorXInc;
                y1 += cursorYInc;
                gl.glBegin(GL.GL_QUAD_STRIP);
                double colTexCoor = tile.minColTexCoor;
                float t = displayTMin;
                int nSlices = tile.nCols;
                float[] depths0 = getDepthTraceFromIntervalVelocities(t, tInc, nSlices, x0, y0, velocityModel);
                float[] depths1 = getDepthTraceFromIntervalVelocities(t, tInc, nSlices, x1, y1, velocityModel);

                for(n = 0, col = tile.croppedColMin; col <= tile.croppedColMax; n++, col++, colTexCoor += dColTexCoor)
                {
                    gl.glTexCoord2d(colTexCoor, rowTexCoor);
                    gl.glVertex3f(x0, y0, depths0[n]);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                    gl.glVertex3f(x1, y1, depths1[n]);
                }
                gl.glEnd();
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawTextureTileDepthSurfaceInTime", " at row " + row + " col " + col, e);
        }
    }

    private float[] getDepthTraceFromIntervalVelocities(float timeMin, float timeInc, int nSlices, float x, float y, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            int row = velocityModel.getNearestBoundedRowCoor(y);
            int col = velocityModel.getNearestBoundedColCoor(x);
            int nSurfaces = velocityModel.getNSurfaces();
            float[] traceTimeValues = velocityModel.getSurfaceTimes(row, col);
            float[] traceIntervalVelocities = velocityModel.getIntervalVelocities(row, col);
            float[] traceDepthValues = velocityModel.getSurfaceDepths(row, col);
            float depthDatum = velocityModel.depthDatum;
            int nBotSurface = 0;
            float timeDatum = velocityModel.timeDatum;
            float topTime = timeDatum;
            float botTime = traceTimeValues[0];
            float topDepth = depthDatum;
            float botDepth = traceDepthValues[0];
            float intervalVelocity = traceIntervalVelocities[0];
            float[] depths = new float[nSlices];
            float t = timeMin;
            for(int n = 0; n < nSlices; n++, t += timeInc)
            {
                while (t > botTime && nBotSurface < nSurfaces - 1)
                {
                    topTime = botTime;
                    topDepth = botDepth;
                    nBotSurface++;
                    botTime = traceTimeValues[nBotSurface];
                    botDepth = traceDepthValues[nBotSurface];
                    intervalVelocity = traceIntervalVelocities[nBotSurface];
                }
                depths[n] = topDepth + (t - topTime) * intervalVelocity;
            }
            return depths;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getDepthTraceFromIntervalVelocities", e);
            return null;
        }
    }

    public void drawTextureTileDepthSurfaceInTime(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        int dir = tile.textureTiles.dir;
        if(dir == StsCursor3d.ZDIR) return;
        if(velocityModel == null) return;

        //TODO correctly implement drawTextureTileDepthSurfaceInTimeWithVolume: need smart interpolator
        // if(velocityModel.inputVelocityVolume != null)
        //     drawTextureTileDepthSurfaceInTimeWithVolume(volume, velocityModel, tile, gl, dir);
        // else
            drawTextureTileDepthSurfaceInTimeWithIntervalVelocities(volume, velocityModel, tile, gl, dir);
    }

    public void drawTextureTileDepthSurfaceInTimeWithVolume(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
    {
        int row = -1, col = -1;
        try
        {
            if(debug) tile.textureTileDebug(this, "drawTextureTileDepthSurfaceInTimeWithVolume", "", TD_DEPTH, TD_TIME, dir, nullByte);
            StsSeismicVolume velocityVolume = velocityModel.velocityVolume;
            gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            double rowTexCoor = tile.minRowTexCoor;
            double dRowTexCoor = tile.dRowTexCoor;
            double dColTexCoor = tile.dColTexCoor;
            double[] xyzMin = tile.xyzPlane[0];
            float displayXMin = (float) xyzMin[0];
            float displayYMin = (float) xyzMin[1];
            float displayZMin = (float) xyzMin[2];
            double[] xyzMax = tile.xyzPlane[2];

            float depthDatum = velocityModel.depthDatum;
            float timeDatum = velocityModel.timeDatum;
            int[] rowCol = velocityVolume.getCursorDataNRowCols(dir);
            //		int nDataRows = rowCol[0];
            int nDataCols = rowCol[1];
            //		int volumeRowInc = 0;
            //		int volumeColInc = 0;
            byte[] velocityBytes;
            int velocityRow, velocityCol;
            float cursorXInc, cursorYInc;
            float velocityTInc = velocityVolume.zInc;
            if(dir == XDIR)
            {
                velocityCol = velocityVolume.getNearestBoundedColCoor(displayXMin);
                velocityBytes = velocityVolume.readBytePlaneData(dir, velocityCol);
                cursorXInc = 0;
                cursorYInc = volume.getYInc();
                //			volumeRowInc = 1;
            }
            else // dirNo == StsCursor3d.YDIR
            {
                velocityRow = velocityVolume.getNearestBoundedRowCoor(displayYMin);
                velocityBytes = velocityVolume.readBytePlaneData(dir, velocityRow);
                cursorXInc = volume.getXInc();
                cursorYInc = 0;
                //			volumeColInc = 1;
            }
            double zInc = volume.getZInc();
            double tInc = velocityModel.getZInc();
            float x1 = displayXMin;
            float y1 = displayYMin;
            float[] velocityTrace1 = getVelocityTrace(velocityBytes, dir, x1, y1, velocityVolume);
            for(row = tile.croppedRowMin; row < tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
            {
                float x0 = x1;
                float y0 = y1;
                float[] velocityTrace0 = velocityTrace1;
                x1 += cursorXInc;
                y1 += cursorYInc;
                velocityTrace1 = getVelocityTrace(velocityBytes, dir, x1, y1, velocityVolume);
                int i = row * nDataCols + tile.colMin;
                gl.glBegin(GL.GL_QUAD_STRIP);
                double colTexCoor = tile.minColTexCoor;
                float z = displayZMin;

                for(col = tile.croppedColMin; col <= tile.croppedColMax; col++, z += zInc, colTexCoor += dColTexCoor, i++)
                {
                    velocityModel.getT(x0, y0, z, 0.0f);
                    //TODO implement
                    float z0 = 0.0f;
                    //float v0 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace0);
                    //float z0 = (v0 * (t - timeDatum) + depthDatum);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor);
                    gl.glVertex3d(x0, y0, z0);
                    //TODO implement
                    float z1 = 0.0f;
                    //float v1 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace1);
                    //float z1 = (v1 * (t - timeDatum) + depthDatum);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                    gl.glVertex3d(x1, y1, z1);
                }
                gl.glEnd();
                //			volumeRow += volumeRowInc;
                //			volumeCol += volumeColInc;
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawTextureTileDepthSurfaceInTime", " at row " + row + " col " + col, e);
        }
    }
    /*
    public void drawTextureTileDepthSurfaceInTime(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dirNo)
    {
        float cursorXInc, cursorYInc;
        if(dirNo == StsCursor3d.ZDIR) return;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        double[] xyz = tile.xyzPlane[0];
        double x1 = xyz[0];
        double y1 = xyz[1];
        double t1 = xyz[2];
        int volumeRow = volume.getNearestBoundedRowCoor((float) y1);
        int volumeCol = volume.getNearestBoundedColCoor((float) x1);
        StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
        if(velocityVolume == null) return;
        float depthMin = velocityModel.depthDatum;
        int volumeRowInc = 0;
        int volumeColInc = 0;
        if(dirNo == XDIR)
        {
            cursorXInc = 0;
            cursorYInc = volume.getYInc();
            volumeRowInc = 1;
        }
        else // dirNo == StsCursor3d.YDIR
        {
            cursorXInc = volume.getXInc();
            cursorYInc = 0;
            volumeColInc = 1;
        }
        double tInc = velocityModel.getZInc();
        for(int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            double x0 = x1;
            double y0 = y1;
            x1 += cursorXInc;
            y1 += cursorYInc;

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            double t = t1 + tile.colMin * tInc;

            for(int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
            {
                float v0 = velocityVolume.getValue(volumeRow, volumeCol, col);
                float z0 = (float) (v0 * t + depthMin);
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3d(x0, y0, z0);
                float v1 = velocityVolume.getValue(volumeRow + volumeRowInc, volumeCol + volumeColInc, col);
                float z1 = (float) (v1 * t + depthMin);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex3d(x1, y1, z1);
            }
            gl.glEnd();
            volumeRow += volumeRowInc;
            volumeCol += volumeColInc;
        }
    }
    */
    public void drawTextureTileDepthSurfaceInTimeWithIntervalVelocities(StsRotatedGridBoundingBox volume, StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
    {
        int n = 0, row = -1, col = -1;
        try
        {
            if(debug) tile.textureTileDebug(this, "drawTextureTileDepthSurfaceInTimeWithIntervalVelocities", "", TD_DEPTH, TD_TIME, dir, nullByte);
            double rowTexCoor = tile.minRowTexCoor;
            double dRowTexCoor = tile.dRowTexCoor;
            double dColTexCoor = tile.dColTexCoor;
            double[] xyzMin = tile.xyzPlane[0];
            float displayXMin = (float) xyzMin[0];
            float displayYMin = (float) xyzMin[1];
            float displayZMin = (float) xyzMin[2];
            double[] xyzMax = tile.xyzPlane[2];

            float cursorXInc, cursorYInc;
            if(dir == XDIR)
            {
                cursorXInc = 0;
                cursorYInc = volume.getYInc();
            }
            else // dirNo == StsCursor3d.YDIR
            {
                cursorXInc = volume.getXInc();
                cursorYInc = 0;
            }
            float zInc = volume.getZInc();
            float x1 = displayXMin;
            float y1 = displayYMin;
            gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            for(row = tile.croppedRowMin; row < tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
            {
                float x0 = x1;
                float y0 = y1;
                x1 += cursorXInc;
                y1 += cursorYInc;
                gl.glBegin(GL.GL_QUAD_STRIP);
                double colTexCoor = tile.minColTexCoor;
                float z = displayZMin;
                int nSlices = tile.nCols;
                float[] times0 = getTimeTraceFromIntervalVelocities(z, zInc, nSlices, x0, y0, velocityModel);
                float[] times1 = getTimeTraceFromIntervalVelocities(z, zInc, nSlices, x1, y1, velocityModel);

                for(n = 0, col = tile.croppedColMin; col <= tile.croppedColMax; n++, col++, colTexCoor += dColTexCoor)
                {
                    gl.glTexCoord2d(colTexCoor, rowTexCoor);
                    gl.glVertex3f(x0, y0, times0[n]);
                    gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                    gl.glVertex3f(x1, y1, times1[n]);
                }
                gl.glEnd();
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawTextureTileDepthSurfaceInTime", " at row " + row + " col " + col, e);
        }
    }

    static public float[] getTimeTraceFromIntervalVelocities(float depthMin, float depthInc, int nSlices, float x, float y, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            int row = velocityModel.getNearestBoundedRowCoor(y);
            int col = velocityModel.getNearestBoundedColCoor(x);
            int nSurfaces = velocityModel.getNSurfaces();
            float[] traceTimeValues = velocityModel.getSurfaceTimes(row, col);
            float[] traceIntervalVelocities = velocityModel.getIntervalVelocities(row, col);
            float[] traceDepthValues = velocityModel.getSurfaceDepths(row, col);
            float depthDatum = velocityModel.depthDatum;
            int nBotSurface = 0;
            float timeDatum = velocityModel.timeDatum;
            float topTime = timeDatum;
            float botTime = traceTimeValues[0];
            float topDepth = depthDatum;
            float botDepth = traceDepthValues[0];
            float intervalVelocity = traceIntervalVelocities[0];
            float[] times = new float[nSlices];
            float z = depthMin;
            for(int n = 0; n < nSlices; n++, z += depthInc)
            {
                while (z > botDepth && nBotSurface < nSurfaces - 1)
                {
                    topDepth = botDepth;
                    topTime = botTime;
                    nBotSurface++;
                    botDepth = traceDepthValues[nBotSurface];
                    botTime = traceTimeValues[nBotSurface];
                    intervalVelocity = traceIntervalVelocities[nBotSurface];
                }
                times[n] = topTime + (z - topDepth) / intervalVelocity;
            }
            return times;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsCursor3dTexture.class, "getDepthTraceFromIntervalVelocities", e);
            return null;
        }
    }
    static public float getTimeFromIntervalVelocities(float[] xyz)
    {
        try
        {
            StsSeismicVelocityModel velocityModel = StsModel.getCurrentProject().velocityModel;
            int row = velocityModel.getNearestBoundedRowCoor(xyz[1]);
            int col = velocityModel.getNearestBoundedColCoor(xyz[0]);
            int nSurfaces = velocityModel.getNSurfaces();
            float[] traceTimeValues = velocityModel.getSurfaceTimes(row, col);
            float[] traceIntervalVelocities = velocityModel.getIntervalVelocities(row, col);
            float[] traceDepthValues = velocityModel.getSurfaceDepths(row, col);
            float depthDatum = velocityModel.depthDatum;
            int nBotSurface = 0;
            float timeDatum = velocityModel.timeDatum;
            float topTime = timeDatum;
            float botTime = traceTimeValues[0];
            float topDepth = depthDatum;
            float botDepth = traceDepthValues[0];
            float intervalVelocity = traceIntervalVelocities[0];
            float z = xyz[2];
            while (z > botDepth && nBotSurface < nSurfaces - 1)
            {
                topDepth = botDepth;
                topTime = botTime;
                nBotSurface++;
                botDepth = traceDepthValues[nBotSurface];
                botTime = traceTimeValues[nBotSurface];
                intervalVelocity = traceIntervalVelocities[nBotSurface];
            }
            return topTime + (z - topDepth) / intervalVelocity;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsCursor3dTexture.class, "getDepthTraceFromIntervalVelocities", e);
            return StsParameters.nullValue;
        }
    }

    static public void textureDebug(ByteBuffer dataBuffer, int nTextureRows, int nTextureCols, Object object, String methodName, String message, byte nullByte)
    {
        int[] rowColSample = textureDebugGetFirstNonNullRowColSample(dataBuffer, nTextureRows, nTextureCols, nullByte);
        if(rowColSample != null)
            StsException.systemDebug(object, methodName, message + " First non-null byte: dataBuffer[" + rowColSample[0] + "][" + rowColSample[1] + "] = " + rowColSample[2]);
        else
            StsException.systemDebug(object, methodName, message + " NO non-null bytes");
    }

    static public int[] textureDebugGetFirstNonNullRowColSample(ByteBuffer dataBuffer, int nTextureRows, int nTextureCols, byte nullByte)
    {
        if(dataBuffer == null) return null;

        dataBuffer.rewind();
        byte sample = nullByte;
        int row = -1, col = -1;
        for(row = 0; row < nTextureRows; row++)
        {
            for(col = 0; col < nTextureCols; col++)
            {
                sample = dataBuffer.get();
                if(sample != nullByte) break;
            }
            if(sample != nullByte) break;
        }
        dataBuffer.rewind();
        if(sample != nullByte)
            return new int[] { row, col, sample };
        else
            return null;
    }

    static public void textureDebug(byte[] data, int nTextureRows, int nTextureCols, Object object, String methodName, String message, byte nullByte)
    {
        if(data == null)
        {
            StsException.systemDebug(object, methodName, message + " byte data is null!");
            return;
        }

        byte sample = nullByte;
        int row = -1, col = -1, i = 0;
        for(row = 0; row < nTextureRows; row++)
        {
            for(col = 0; col < nTextureCols; col++)
            {
                sample = data[i++];
                if(sample != nullByte) break;
            }
            if(sample != nullByte) break;
        }
        if(sample != nullByte)
            StsException.systemDebug(object, methodName, message + " First non-null byte: dataBuffer[" + row + "][" + col + "] = " + sample);
        else
            StsException.systemDebug(object, methodName, message + " NO non-null bytes");
    }
}

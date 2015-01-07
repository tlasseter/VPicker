package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.Crossplot.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.Seismic.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.nio.*;

public class StsSeismicCursorSection extends StsCursor3dVolumeTexture
{
    public StsSeismicVolume seismicVolume;
    transient StsTraceSet traces;
    transient public ByteBuffer dataBuffer = null;
    transient public ByteBuffer subVolumeDataBuffer = null;
    transient public boolean cropChanged = false;
    transient boolean isVisible = true;
    transient byte[] transparentData = null;
    transient private StsSeismicVolumeClass seismicVolumeClass;

    /** minimum size in pixels of a crossplot point */
    static private int minPointSize = 4;

    static private final boolean debug = false;

    public StsSeismicCursorSection()
    {
    }

    public StsSeismicCursorSection(StsModel model, StsSeismicVolume seismicVolume, StsCursor3d cursor3d, int dir)
    {
        this.seismicVolume = seismicVolume;
        initialize(model, cursor3d, dir);
    }

    public boolean initialize(StsModel model, StsCursor3d cursor3d, int dir)
    {
        super.initialize(model, cursor3d, dir);
        if(seismicVolume == null) return false;
        seismicVolumeClass = seismicVolume.getSeismicVolumeClass();
        isPixelMode = getIsPixelMode();
        return true;
    }

    protected boolean getIsPixelMode()
    {
        return seismicVolumeClass.getIsPixelMode();
    }

    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    public boolean getUseShader() { return seismicVolume.getSeismicVolumeClass().getContourColors(); }

    /**
     * This results in texture being deleted on next displayTexture call.
     * Delete the texture unconditionally unless plane is ZDIR.
     * If is ZDIR, delete texture if not displaying basemap or if basemap has changed.
     */
    public boolean textureChanged()
    {
        if (seismicVolume == null)
        {
            textureChanged = false;
            return false;
        }
        if(dirNo != StsCursor3d.ZDIR)
            textureChanged = true;
        else // dirNo == ZDIR
        {
            if(seismicVolume.displayingAttribute())
                textureChanged = seismicVolume.hasAttributeChanged();
            else
                textureChanged = true;
        }
        if(debug) StsException.systemDebug(this, "textureChanged");
        return textureChanged;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void cropChanged()
    {
        /*
                 if (textureTiles == null)
           {
            if(!checkTextureTiles()) return;
           }
                 else
           */
        if (textureTiles == null) return;
//		addTextureToDeleteList();
        textureTiles.cropChanged();
    }

    public void subVolumeChanged()
    {
    }

    public Class getDisplayableClass() { return StsSeismicVolume.class; }

    public boolean canDisplayClass(Class c) { return StsSeismicVolume.class.isAssignableFrom(c); }

    public boolean isDisplayableObject(Object object)
    {
        return (object instanceof StsSeismicVolume);
    }

    public boolean isDisplayingObject(Object object)
    {
        if (seismicVolume == object) return true;
        if (this.textureTiles != null && object == textureTiles.cropVolume) return true;
        return false;
    }

    public boolean setDirCoordinate(float dirCoordinate)
    {
        if (this.dirCoordinate == dirCoordinate) return false;
        this.dirCoordinate = dirCoordinate;
        textureChanged();
        StsXPolygonAction.clearSeismicPoints();
        return true;
    }

    public boolean setObject(Object object)
    {
        if (object == seismicVolume) return false;

        if (debug) StsException.systemDebug(this, "setObject", "Object changed from " +
            StsMainObject.getObjectName(seismicVolume) + " to: " +
            StsMainObject.getObjectName((StsMainObject) object));

        if (isDisplayableObject(object))
        {
            seismicVolume = (StsSeismicVolume) object;
            isVisible = true;
        }
        else if (object != null)
            return false;
        else // object == null
        {
            isVisible = false;
            seismicVolume = null;
        }
        textureChanged();
        return true;
    }

    public Object getObject()
    {
        return seismicVolume;
    }

    public boolean isVisible()
    {
        if (!isVisible) return false;
        if (seismicVolume == null) return false;
        if (!seismicVolume.canDisplayZDomain()) return false;
        return seismicVolume.getIsVisibleOnCursor();
    }

    protected void checkTextureAndGeometryChanges(StsGLPanel3d glPanel3d, boolean is3d)
    {
        if (isPixelMode != seismicVolume.getSeismicVolumeClass().getIsPixelMode())
        {
            textureChanged = true;
            isPixelMode = !isPixelMode;
        }

        if (seismicVolumeChanged())
        {
            textureChanged = true;
        }

        byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
        if (projectZDomain != zDomain)
        {
            geometryChanged = true;
            zDomain = projectZDomain;
        }
    }

    protected boolean initializeTextureTiles(StsGLPanel3d glPanel3d)
    {
        // if (!glPanel3d.initialized) return false;
        if (seismicVolume == null) return false;

        if (textureTiles == null)
        {
            StsCropVolume subVolume = glPanel3d.model.getProject().getCropVolume();
            textureTiles = StsTextureTiles.constructor(model, this, dirNo, seismicVolume, isPixelMode, subVolume);
            if (textureTiles == null) return false;
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
//            geometryChanged = true;
        }
        else if (!textureTiles.isSameSize(seismicVolume))
        {
            textureTiles.constructTiles(seismicVolume);
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
//            geometryChanged = true;
        }

        if (textureTiles.shaderChanged())
        {
            textureChanged = true;
        }
        textureTiles.setTilesDirCoordinate(dirCoordinate);

        return true;
    }

    protected void computeTextureData()
    {
        if (debug)
        {
            int nPlane = seismicVolume.getCursorPlaneIndex(dirNo, dirCoordinate);
            System.out.println("Texture changed, reading new seismic texture." + "dir: " + dirNo + " " + " coor: " +
                dirCoordinate + " plane: " + nPlane);
        }
        dataBuffer = seismicVolume.readByteBufferPlane(dirNo, dirCoordinate);
        if (dataBuffer == null)
        {
            textureChanged = false;
            return;
        }
        byte[] subVolumePlane = cursor3d.getSubVolumePlane(dirNo, dirCoordinate, seismicVolume, seismicVolume.getZDomain());
        subVolumeDataBuffer = applySubVolume(subVolumePlane);
        if(subVolumeDataBuffer == null)
        {
            StsException.systemError(this, "computeTextureData", "Failed to compute subVolumeDataBuffer");
            subVolumeDataBuffer = dataBuffer;
        }
        if (debug)
        {
            String message = "read dataBuffer for volume: " + seismicVolume.getName() + " dirNo: " + dirNo;
            String methodName = "displayTexture";
            textureDebug(subVolumeDataBuffer, nTextureRows, nTextureCols, this, methodName, message, nullByte);
            StsException.systemDebug(this, "displayTexture", "Displaying changed seismic texture 3d in dirNo " + dirNo + ".");
            {
                subVolumeDataBuffer.rewind();

                byte sample = -1;
                int row = -1, col = -1;
                for(row = 0; row < nTextureRows; row++)
                {
                    for(col = 0; col < nTextureCols; col++)
                    {
                        byte value = dataBuffer.get();
                        if(value > 0)
                        {
                            sample = value;
                            break;
                        }
                    }
                }
                subVolumeDataBuffer.rewind();
                if(sample != -1)
                    StsException.systemDebug(this, "displayTexture", "read dataBuffer for volume: " + seismicVolume.getName() + " dirNo: " + dirNo +
                        " First byte value > 0: dataBuffer[" + row + "][" + col + "] = " + sample);
                else
                    StsException.systemDebug(this, "displayTexture", "read dataBuffer for volume: " + seismicVolume.getName() + " dirNo: " + dirNo +
                        " NO positive bytes");
            }
        }
    }

    protected boolean enableGLState(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        super.enableGLState(glPanel3d, gl, is3d);
        if (debug) StsException.systemDebug(this, "displayTexture", "setGLColorList called with shader: " + textureTiles.shader);
        return seismicVolume.setGLColorList(gl, false, dirNo, textureTiles.shader);
    }

    protected void displayTiles3d(StsGLPanel3d glPanel3d, GL gl)
    {
        textureTiles.displayTiles(this, gl, isPixelMode, subVolumeDataBuffer, nullByte);
        subVolumeDataBuffer = null;
    }

    protected void displayTiles2d(StsGLPanel3d glPanel3d, GL gl)
    {
        StsViewCursor viewCursor = (StsViewCursor)glPanel3d.getView();
        textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, subVolumeDataBuffer, nullByte);
        subVolumeDataBuffer = null;
    }

    private boolean seismicVolumeChanged()
    {
        StsSeismicVolumeClass seismicVolumeClass = (StsSeismicVolumeClass) model.getStsClass(seismicVolume.getClass());
        StsSeismicVolume currentSeismicVolume = seismicVolumeClass.getCurrentSeismicVolume();
        if (seismicVolume == currentSeismicVolume)
            return seismicVolume.hasAttributeChanged();

        if (currentSeismicVolume == null) return true;
        return true;
    }

    private byte[] transparentData()
    {
        if (transparentData != null) return transparentData;
        transparentData = new byte[nTextureRows * nTextureCols];
        return transparentData;
    }

    /**
     * The subVolumePlane has the dimensions of the cursorBoundingBox it cuts through;
     * The dataBuffer has the dimensions of the corresponding seismic/virtual volume.
     * If they are the same size, apply the dataBuffer directly to the subVolumePlane;
     * if not get the grid coordinates of the dataBuffer relative to the subVolume plane
     * and apply.
     */
    protected ByteBuffer applySubVolume(byte[] subVolumePlane)
    {
        int nCursorCols;
        int cursorIndex, seismicIndex;

        if (dataBuffer == null) return null;
        if (subVolumePlane == null) return dataBuffer;
        if (!seismicVolume.getSeismicVolumeClass().getDisplayOnSubVolumes())
            return dataBuffer;

        int nPlanePoints = nTextureRows * nTextureCols;
        ByteBuffer subVolumeData = ByteBuffer.allocateDirect(nPlanePoints);

        if (nPlanePoints == subVolumePlane.length)
        {
            for (int n = 0; n < nPlanePoints; n++)
            {
                if (subVolumePlane[n] == 0)
                    subVolumeData.put((byte) -1);
                else
                    subVolumeData.put(dataBuffer.get(n));
            }
            return subVolumeData;
        }
        else
        {
            nCursorCols = cursorBoundingBox.getNCursorCols(dirNo);

            int seismicRowStart = cursorBoundingBox.getSubVolumeCursorRowMin(seismicVolume, dirNo);
            int seismicRowEnd = cursorBoundingBox.getSubVolumeCursorRowMax(seismicVolume, dirNo);
            int seismicColStart = cursorBoundingBox.getSubVolumeCursorColMin(seismicVolume, dirNo);
            int seismicColEnd = cursorBoundingBox.getSubVolumeCursorColMax(seismicVolume, dirNo);

            cursorIndex = 0;
            seismicIndex = 0;

            try
            {
                for (int row = seismicRowStart; row <= seismicRowEnd; row++)
                {
                    cursorIndex = row * nCursorCols + seismicColStart;
                    for (int col = seismicColStart; col <= seismicColEnd; col++, cursorIndex++, seismicIndex++)
                    {
                        byte planeByte = dataBuffer.get();
                        if (subVolumePlane[cursorIndex] != 0)
                            subVolumeData.put(planeByte);
                        else
                            subVolumeData.put(StsParameters.nullByte);
                    }
                }
                return subVolumeData;
            }
            catch (Exception e)
            {
                StsException.outputException("StsSeismicCursorSection.applySubVolume() failed.", e, StsException.WARNING);
                return null;
            }
        }
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        drawTextureTileSurface(seismicVolume, tile, gl, is3d, seismicVolume.getZDomain());
    }

    /*
      public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
      {
       line2d.drawTextureTileSurface(tile, gl, dirNo, false);
      }
      */
    public void display(StsGLPanel3d glPanel3d, boolean is3d)
    {
        if (seismicVolume == null)
            return;

        if (!seismicVolume.canDisplayZDomain()) return;

        GL gl = glPanel3d.getGL();
        if (is3d)
        {
            displayCrossplotPoints(gl, glPanel3d);
        }
        else
        {
            StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getView();
            displayCrossplotPoints2d(gl, glPanel3d, viewCursor);
            displaySeismicPolygons2d(gl, glPanel3d, viewCursor.axesFlipped);
            if (!displayWiggles()) return;
            StsWiggleDisplayProperties wiggleDisplayProperties = seismicVolume.getWiggleDisplayProperties();
                boolean smoothCurve = wiggleDisplayProperties.getWiggleSmoothCurve();
            if(dirNo == XDIR)
                displayCubicWiggleTracesCol(glPanel3d, gl, viewCursor);
            else if(dirNo == YDIR)
                displayCubicWiggleTracesRow(glPanel3d, gl, viewCursor);
        }
    }

    private void displayCubicWiggleTracesRow(StsGLPanel3d glPanel3d, GL gl, StsViewCursor viewCursor)
    {
        float[][] axisRanges = viewCursor.axisRanges;
        float xInc = cursorBoundingBox.xInc;
       // If density is less than 1:4 traces to pixels, display wiggles
        if (((axisRanges[0][1] - axisRanges[0][0]) / xInc) > glPanel3d.getWidth() / getWiggleToPixelRatio())
              return;

        StsWiggleDisplayProperties wiggleProperties = seismicVolume.getWiggleDisplayProperties();

        int nSlices = seismicVolume.nSlices;

        int row = seismicVolume.getNearestRowCoor(dirCoordinate);
        if(traces == null)
        {
            FloatBuffer floatBuffer = seismicVolume.getRowPlaneFloatBuffer(row);
            if(floatBuffer == null) return;
            int nCols = seismicVolume.nCols;
            float[][] traceData = new float[nCols][nSlices];
            for(int col = 0; col < nCols; col++)
                floatBuffer.get(traceData[col]);
            StsMath.computeNormalizedRMSAmplitude(traceData);
            traces = new StsTraceSet(traceData, seismicVolume, 5, wiggleProperties);
        }
        traces.initializeDraw(viewCursor, wiggleProperties);

        StsColor.BLACK.setGLColor(gl);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glLineWidth(0.5f);



        float displayXMin = axisRanges[0][0];
        float displayXMax = axisRanges[0][1];
        int colMin = seismicVolume.getFloorBoundedColCoor(displayXMin);
        int colMax = seismicVolume.getCeilingBoundedColCoor(displayXMax);
        float colXMin = seismicVolume.getXCoor(colMin);

        float x = colXMin;
        for (int col = colMin; col <= colMax; col++, x += xInc)
        {
            traces.displayInterpolatedPoints(gl, col, x);
        }
        gl.glEnable(GL.GL_LIGHTING);
    }

    private void displayCubicWiggleTracesCol(StsGLPanel3d glPanel3d, GL gl, StsViewCursor viewCursor)
    {
        float[][] axisRanges = viewCursor.axisRanges;
        float xInc = cursorBoundingBox.xInc;
       // If density is less than 1:4 traces to pixels, display wiggles
        if (((axisRanges[0][1] - axisRanges[0][0]) / xInc) > glPanel3d.getWidth() / getWiggleToPixelRatio())
              return;

        StsWiggleDisplayProperties wiggleProperties = seismicVolume.getWiggleDisplayProperties();

        int nSlices = seismicVolume.nSlices;

        int col = seismicVolume.getNearestColCoor(dirCoordinate);
        if(traces == null)
        {
            int nRows = seismicVolume.nRows;
            float[][] traceData = new float[nRows][nSlices];
            for(int row = 0; row < nRows; row++)
                traceData[row] = seismicVolume.getTraceValues(row, col);
            StsMath.computeNormalizedRMSAmplitude(traceData);
            traces = new StsTraceSet(traceData, seismicVolume, 5, wiggleProperties);
        }
        traces.initializeDraw(viewCursor, wiggleProperties);

        float displayYMin = axisRanges[0][0];
        float displayYMax = axisRanges[0][1];
        int rowMin = StsMath.floor(seismicVolume.getRowCoor(displayYMin));
        int rowMax = StsMath.ceiling(seismicVolume.getRowCoor(displayYMax));
        float rowYMin = seismicVolume.getYCoor(rowMin);

        StsColor.BLACK.setGLColor(gl);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glLineWidth(0.5f);
        float y = rowYMin;
        float yInc = cursorBoundingBox.yInc;
        for (int row = rowMin; row <= rowMax; row++, y += yInc)
            traces.displayInterpolatedPoints(gl, row, y);
        gl.glEnable(GL.GL_LIGHTING);
    }

    private boolean displayWiggles()
    {
        StsSeismicVolumeClass volumeClass = (StsSeismicVolumeClass) model.getStsClass(seismicVolume.getClass());
        return volumeClass.getDisplayWiggles();
    }

    private int getWiggleToPixelRatio()
    {
        StsSeismicVolumeClass volumeClass = (StsSeismicVolumeClass) model.getStsClass(seismicVolume.getClass());
        return volumeClass.getWiggleToPixelRatio();
    }

    public void displayCrossplotPoints(GL gl, StsGLPanel3d glPanel3d)
    {
        displayCrossplotPoints(gl, glPanel3d, StsXPolygonAction.getSeismicCrossplotPoints(dirNo, dirCoordinate));
    }

    public void displayCrossplotPoints(GL gl, StsGLPanel3d glPanel3d, StsCrossplotPoint[] crossplotPoints)
    {
        if (crossplotPoints == null)
        {
            return;
        }
        int nPoints = crossplotPoints.length;

        gl.glDisable(GL.GL_LIGHTING);

        glPanel3d.setViewShift(gl, 1.0);
        for (int n = 0; n < nPoints; n++)
        {
            StsGLDraw.drawPoint(gl, crossplotPoints[n].getVolumeXYZ(), StsColor.BLACK, minPointSize + 2);

        }
        glPanel3d.setViewShift(gl, 2.0);
        for (int n = 0; n < nPoints; n++)
        {
            StsGLDraw.drawPoint(gl, crossplotPoints[n].getVolumeXYZ(), crossplotPoints[n].stsColor, minPointSize);

        }
        glPanel3d.resetViewShift(gl);

        gl.glEnable(GL.GL_LIGHTING);

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPopMatrix();
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    public void displayCrossplotPoints2d(GL gl, StsGLPanel3d glPanel3d, StsViewCursor viewCursor)
    {
        displayCrossplotPoints2d(gl, glPanel3d, StsXPolygonAction.getSeismicCrossplotPoints(dirNo, dirCoordinate), viewCursor.axesFlipped);
    }

    public void displayCrossplotPoints2d(GL gl, StsGLPanel3d glPanel3d, StsCrossplotPoint[] crossplotPoints, boolean axesFlipped)
    {
        // Because of the data layout, rows are vertically down for XDIR and YDIR cursor planes
        if (crossplotPoints == null)
        {
            return;
        }

        float rowMinF, colMinF, rowMaxF, colMaxF;

        if (!(glPanel3d.getView() instanceof StsViewCursor))
        {
            return;
        }

        StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getView();
        float[][] axisRanges = viewCursor.getAxisRanges();
        boolean flip = false;
        if (dirNo == StsCursor3d.XDIR)
        {
            flip = true;
            colMinF = seismicVolume.getSliceCoor(axisRanges[1][0]);
            rowMinF = seismicVolume.getRowCoor(axisRanges[0][0]);
            colMaxF = seismicVolume.getSliceCoor(axisRanges[1][1]);
            rowMaxF = seismicVolume.getRowCoor(axisRanges[0][1]);
        }
        else if (dirNo == StsCursor3d.YDIR)
        {
            flip = true;
            colMinF = seismicVolume.getSliceCoor(axisRanges[1][0]);
            rowMinF = seismicVolume.getColCoor(axisRanges[0][0]);
            colMaxF = seismicVolume.getSliceCoor(axisRanges[1][1]);
            rowMaxF = seismicVolume.getColCoor(axisRanges[0][1]);
        }
        else
        {
            if (!axesFlipped)
            {
                flip = false;
                rowMinF = seismicVolume.getRowCoor(axisRanges[1][0]);
                colMinF = seismicVolume.getColCoor(axisRanges[0][0]);
                rowMaxF = seismicVolume.getRowCoor(axisRanges[1][1]);
                colMaxF = seismicVolume.getColCoor(axisRanges[0][1]);
            }
            else
            {
                flip = true;
                colMinF = seismicVolume.getRowCoor(axisRanges[1][0]);
                rowMinF = seismicVolume.getColCoor(axisRanges[0][0]);
                colMaxF = seismicVolume.getRowCoor(axisRanges[1][1]);
                rowMaxF = seismicVolume.getColCoor(axisRanges[0][1]);
            }

        }

        GLU glu = glPanel3d.getGLU();

        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPushMatrix();
        gl.glLoadIdentity();

        float viewPortHeight, viewPortWidth;
        float pixelsPerRow, pixelsPerCol;
        float halfBoxWidth = 0.5f;
        float halfBoxHeight = 0.5f;
        float halfOutlineBoxWidth, halfOutlineBoxHeight;
        if (flip)
        {
            viewPortWidth = glPanel3d.getWidth();
            pixelsPerRow = Math.abs(viewPortWidth / (rowMaxF - rowMinF + 1));
            if (pixelsPerRow < minPointSize)
            {
                halfBoxWidth = minPointSize / pixelsPerRow / 2;
            }
            halfOutlineBoxWidth = halfBoxWidth + 1.0f / pixelsPerRow;

            viewPortHeight = glPanel3d.getHeight();
            pixelsPerCol = Math.abs(viewPortHeight / (colMaxF - colMinF + 1));
            if (pixelsPerCol < minPointSize)
            {
                halfBoxHeight = minPointSize / pixelsPerCol / 2;
            }
            halfOutlineBoxHeight = halfBoxHeight + 1.0f / pixelsPerCol;

            glu.gluOrtho2D(rowMinF, rowMaxF, colMinF, colMaxF);
        }
        else
        {
            viewPortHeight = glPanel3d.getHeight();
            pixelsPerRow = Math.abs(viewPortHeight / (rowMaxF - rowMinF + 1));
            if (pixelsPerRow < minPointSize)
            {
                halfBoxHeight = minPointSize / pixelsPerRow / 2;
            }
            halfOutlineBoxHeight = halfBoxHeight + 1.0f / pixelsPerRow;

            viewPortWidth = glPanel3d.getWidth();
            pixelsPerCol = Math.abs(viewPortWidth / (colMaxF - colMinF + 1));
            if (pixelsPerCol < minPointSize)
            {
                halfBoxWidth = minPointSize / pixelsPerCol / 2;
            }
            halfOutlineBoxWidth = halfBoxWidth + 1.0f / pixelsPerCol;

            glu.gluOrtho2D(colMinF, colMaxF, rowMinF, rowMaxF);
        }
//        System.out.println("halfBoxWidth: " + halfBoxWidth + " halfBoxHeight: " + halfBoxHeight);

        // save the projection matrix: needs to be set correctly for StsGLDraw.drawPoint
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);

        float rowCenter, colCenter;
        int nPoints = crossplotPoints.length;
//        float[] xyz = new float[3];
        gl.glDisable(GL.GL_LIGHTING);
        for (int n = 0; n < nPoints; n++)
        {
            //TODO fix!  confusion in StsCursorPoint use of rowNum/colNum cursorRow/Col
            //TODO should eliminate from StsCursorPoint and compute locally
            // data row and col are bin numbers so are the lower left corner
            if (!flip)
            {
                rowCenter = crossplotPoints[n].getCursorRow();
                colCenter = crossplotPoints[n].getCursorCol();
            }
            else
            {
                rowCenter = crossplotPoints[n].getCursorCol();
                colCenter = crossplotPoints[n].getCursorRow();
            }

            gl.glBegin(GL.GL_QUADS);
            StsColor.BLACK.setGLColor(gl);
            gl.glVertex2f(colCenter - halfOutlineBoxWidth, rowCenter - halfOutlineBoxHeight);
            gl.glVertex2f(colCenter + halfOutlineBoxWidth, rowCenter - halfOutlineBoxHeight);
            gl.glVertex2f(colCenter + halfOutlineBoxWidth, rowCenter + halfOutlineBoxHeight);
            gl.glVertex2f(colCenter - halfOutlineBoxWidth, rowCenter + halfOutlineBoxHeight);

            crossplotPoints[n].stsColor.setGLColor(gl);
            gl.glVertex2f(colCenter - halfBoxWidth, rowCenter - halfBoxHeight);
            gl.glVertex2f(colCenter + halfBoxWidth, rowCenter - halfBoxHeight);
            gl.glVertex2f(colCenter + halfBoxWidth, rowCenter + halfBoxHeight);
            gl.glVertex2f(colCenter - halfBoxWidth, rowCenter + halfBoxHeight);
            gl.glEnd();
        }
        gl.glEnable(GL.GL_LIGHTING);

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPopMatrix();
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    public void displaySeismicPolygons2d(GL gl, StsGLPanel3d glPanel3d, boolean axesFlipped)
    {
        StsXPolygon[] seismicPolygons = StsXPolygonAction.getSeismicPolygons(dirNo, dirCoordinate);
        if (seismicPolygons == null)
        {
            return;
        }

        if (!(glPanel3d.getView() instanceof StsViewCursor))
        {
            return;
        }

        int nPolygons = seismicPolygons.length;
        for (int n = 0; n < nPolygons; n++)
        {
            StsXPolygon polygon = seismicPolygons[n];
            polygon.draw(gl, null, false);
        }
    }

    public String propertyReadout(StsPoint point)
    {
        StringBuffer stringBuffer = null;
        if (seismicVolume == null)
            return new String("Nothing is currently displayed on the cursor plane.");

        StsObject[] volumes = model.getObjectList(seismicVolume.getClass());
        for (int i = 0; i < volumes.length; i++)
        {
            StsSeismicVolume volume = (StsSeismicVolume) volumes[i];
            if (!volume.getReadoutEnabled() && volume != seismicVolume)
            {
                continue;
            }
            if (stringBuffer == null)
            {
                stringBuffer = new StringBuffer();
            }
            byte byteValue = (byte) volume.getPlaneValue(dirNo, point.v);
            if (byteValue == StsParameters.nullByte)
                stringBuffer.append(" " + volume.getName() + ": null");
            else
                stringBuffer.append(" " + volume.getName() + ": " + volume.getScaledValue(byteValue)+ " "+volume.getUnits());
        }
        if (stringBuffer == null)
        {
            return null;
        }
        return stringBuffer.toString();
    }

    public String logReadout2d(StsPoint2D point)
    {
        int[] dataRowCol = getDataRowCol(point);
        return logReadout2d(dataRowCol);
    }

    public String logReadout2d(float[] rowColF)
    {
        int row = Math.round(rowColF[0]);
        int col = Math.round(rowColF[1]);
        int[] rowCol = new int[]
            {row, col};
        return logReadout2d(rowCol);
    }

    public String logReadout2d(int[] rowCol)
    {
        String valueString;

        if (rowCol != null && isInRange(rowCol[0], rowCol[1]))
        {
            if (dataBuffer == null)
            {
                return "no data";
            }
            byte byteValue = dataBuffer.get(rowCol[0] * nTextureCols + rowCol[1]);
            float value = seismicVolume.getScaledValue(byteValue);
            return Float.toString(value);
        }
        else
        {
            return "not in range";
        }
        /*
                  int row = rowCol[0];
                  int col = rowCol[1];
                  if(dirNo == StsCursor3d.XDIR)
                  {
                      float yLabel = line2d.getNumFromIndex(StsCursor3d.YDIR, (float)col);
                      float zLabel = line2d.getNumFromIndex(StsCursor3d.ZDIR, (float)row);
                      return new String("Line: " + yLabel + " Slice: " + zLabel + " value: " + valueString);
                  }
                  else if(dirNo == StsCursor3d.YDIR)
                  {
                      float xLabel = line2d.getNumFromIndex(StsCursor3d.XDIR, (float)col);
                      float zLabel = line2d.getNumFromIndex(StsCursor3d.ZDIR, (float)row);
                      return new String("Crossline: " + xLabel + " Slice: " + zLabel + " Value: " + valueString);
                  }
                  else
                  {
                      float xLabel = line2d.getNumFromIndex(StsCursor3d.XDIR, (float)col);
                      float yLabel = line2d.getNumFromIndex(StsCursor3d.YDIR, (float)row);
                      return new String("Line: " + yLabel + " Crossline: " + xLabel + " Value: " + valueString);
                  }
           */
    }

    public int[] getDataRowCol(StsPoint2D point)
    {
        int row, col;

        if (dirNo == StsCursor3d.XDIR) // first row is vertically down
        {
            row = (int) seismicVolume.getRowCoor(point.x);
            col = (int) seismicVolume.getSliceCoor(point.y);
        }
        else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
        {
            row = (int) seismicVolume.getColCoor(point.x);
            col = (int) seismicVolume.getSliceCoor(point.y);
        }
        else // dirNo == ZDIR  first row is horizontally across
        {
            row = (int) seismicVolume.getRowCoor(point.y);
            col = (int) seismicVolume.getColCoor(point.x);
        }
        if (isInRange(row, col))
        {
            return new int[]
                {row, col};
        }
        else
        {
            return null;
        }
    }

    /** given 2D coordinates on this 2D cursorSection, return the row col coordinates */
    public float[] getDataRowColF(float x, float y)
    {
        float rowF, colF;

        if (dirNo == StsCursor3d.XDIR) // first row is vertically down
        {
            rowF = seismicVolume.getRowCoor(x);
            colF = seismicVolume.getSliceCoor(y);
        }
        else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
        {
            rowF = seismicVolume.getColCoor(x);
            colF = seismicVolume.getSliceCoor(y);
        }
        else // dirNo == ZDIR  first row is horizontally across
        {
            rowF = seismicVolume.getRowCoor(y);
            colF = seismicVolume.getColCoor(x);
        }
        if (isInRange(rowF, colF))
        {
            return new float[]
                {rowF, colF};
        }
        else
        {
            return null;
        }
    }

    /** given 3D coordinates on this 3d cursorSection, return the row col coordinates */
    public float[] getDataRowColF(StsPoint point)
    {
        float rowF, colF;

        if (dirNo == StsCursor3d.XDIR) // first row is vertically down
        {
            rowF = seismicVolume.getRowCoor(point.getY());
            colF = seismicVolume.getSliceCoor(point.getZ());
        }
        else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
        {
            rowF = seismicVolume.getColCoor(point.getX());
            colF = seismicVolume.getSliceCoor(point.getZ());
        }
        else // dirNo == ZDIR  first row is horizontally across
        {
            rowF = seismicVolume.getRowCoor(point.getY());
            colF = seismicVolume.getColCoor(point.getX());
        }
        if (isInRange(rowF, colF))
        {
            return new float[]
                {rowF, colF};
        }
        else
        {
            return null;
        }
    }

    public String getName()
    {
        return "Cursor view[" + dirNo + "] of: " + seismicVolume.getName();
    }
}

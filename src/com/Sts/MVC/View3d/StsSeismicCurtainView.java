package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining two-dimensional cursor view. The cursor is a reference to the
 * three planes that are isVisible in the 3D view. An object based on this class wuold result in a 2D
 * view of any one of the three cursor planes.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.StsAction;
import com.Sts.Actions.StsDefaultAction3d;
import com.Sts.DB.StsSerializable;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsClassCursorDisplayable;
import com.Sts.Interfaces.StsTextureSurfaceFace;
import com.Sts.MVC.Main;
import com.Sts.MVC.StsActionManager;
import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.StsModel;
import com.Sts.Types.*;
import com.Sts.UI.StsMenuItem;
import com.Sts.UI.StsViewItem;
import com.Sts.Utilities.*;

import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

// TODO should really be subclassed into ViewCursorSeismic and viewCursorPreStack to avoid display confusion
public class StsSeismicCurtainView extends StsView2d implements StsSerializable, StsTextureSurfaceFace, ActionListener
{
    /** well defining the seismic curtain */
    private StsWell well;
    StsSeismicVolume seismicVolume = null;
    /** number of traces in Seismic */
    public int nTraces = -1;
    /** number of samples in each trace */
    public int nSamples = -1;
    /** z coordinates (top down) for texture */
    public float[] zCoordinates;
    /** Display lists should be used (controlled by View:Display Options) */
    boolean useDisplayLists;
    /** Display lists currently being used for surface geometry */
    boolean usingDisplayLists = true;
    transient private StsSeismicCurtain seismicCurtain;
    transient StsSeismicVolumeClass seismicClass = null;
    transient StsTextureTiles textureTiles = null;
    transient boolean textureChanged = true;
    transient boolean isPixelMode = false;
    public byte zDomain = StsParameters.TD_NONE;

    private transient JPopupMenu tp = null;
    private transient StsMenuItem aspectBtn = new StsMenuItem();
    /** vertical plane at a constant x going y+ and z+ directions */
    static final int XDIR = StsCursor3d.XDIR;
    /** vertical plane at a constant y going x+ and z+ directions */
    static final int YDIR = StsCursor3d.YDIR;
    /** horizontal plane at a constant z going x+ and y+ directions */
    static final int ZDIR = StsCursor3d.ZDIR;

    static final byte CORE_LABS_LEGEND = 1;
    static final byte S2S_LEGEND = 2;

    static final String[] compassDirections = new String[]{"E", "NE", "N", "NW", "W", "SW", "S", "SE"};
    /**
     * reverse and flip XY axes so they are closest to looking in a north direction with east to the right.
     * See setViewOrientation() for details of reverse and flip depending on view quadrant (N, S, E, W)
     */
    static final boolean reorient = false;


    static public final String viewNameCurtain = "Seismic Curtain View";
    static public final String shortViewNameCurtain = "Curtain";
    static final long serialVersionUID = 1l;

    public StsSeismicCurtainView()
    {
    }

    public StsSeismicCurtainView(StsModel model, StsWell well)
    {
        this.well = well;
        initialize();
    }

    public boolean initializeView(StsModel model, StsActionManager actionManager)
    {
        try
        {
            this.model = model;

            this.seismicVolume = seismicCurtain.getSeismicVolume();
            this.seismicClass = seismicCurtain.getSeismicVolumeClass();
            // seismicVolume.addActionListener(this);

            int width = seismicCurtain.getNCols() * 2;
            int height = 500;
            glPanel = new StsGLPanel(model, actionManager, width, height, this);
            gl = glPanel.getGL();
            glu = glPanel.getGLU();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initializeView", e);
            return false;
        }
    }


    protected void initializeRange()
    {
        nSamples = seismicCurtain.getNSlices();
        nTraces = seismicCurtain.getNCols();

        float tMin = seismicCurtain.getZMin();
        float tInc = seismicCurtain.getZInc();
        double[] zCoordinates = new double[nSamples];
        float t = tMin;
        for(int n = 0; n < nSamples; n++, t += tInc)
            zCoordinates[n] = t;
        this.zCoordinates = StsMath.convertDoubleToFloatArray(zCoordinates);
        StsGridPoint[] arcGridPoints = seismicCurtain.cellGridCrossingPoints;

        totalAxisRanges = new float[2][2];
        totalAxisRanges[0][0] = 0.f;
        totalAxisRanges[0][1] = arcGridPoints[nTraces - 1].point.getF();
        totalAxisRanges[1][1] = (float) zCoordinates[0];
        totalAxisRanges[1][0] = (float) zCoordinates[nSamples - 1];

        axisRanges = new float[2][2];
        axisRanges[0][0] = 0.f;
        axisRanges[0][1] = totalAxisRanges[0][1];
        axisRanges[1][1] = totalAxisRanges[1][1];
        axisRanges[1][0] = totalAxisRanges[1][0];
    }

    public void initializeTransients(StsGLPanel3d glPanel)
    {
        super.initializeTransients(glPanel);
        initialize();
    }

    static public String getStaticViewName()
    {
        return viewNameCurtain;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameCurtain;
    }

    public void initialize()
    {
        aspectBtn.setMenuActionListener("1:1 Aspect Ratio", this, "aspectRatio", null);
        seismicCurtain = well.getCreateSeismicCurtain();
        seismicClass = seismicCurtain.getSeismicVolumeClass();
        seismicVolume = seismicCurtain.getSeismicVolume();
        seismicVolume.addActionListener(this);
        initializeRange();
    }

    public void computeProjectionMatrix()
    {
        if(axisRanges == null) return;
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    public boolean isViewable()
    {
        return getCursor3d().selectedDirection != StsParameters.NO_MATCH;
    }

    public void viewChanged()
    {
        if(glPanel3d == null) return;
        glPanel3d.viewChanged = true;
    }

    public void display(GLAutoDrawable component)
    {

        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        if(Main.isGLDebug) System.out.println("StsWellSeismicCurtainView.display() called.");
        display2d(glPanel3d);
        displayTexture();
    }


    public void display2d(StsGLPanel3d glPanel3d)
    {
        StsArrayList displayableClasses = model.getCursorDisplayableClasses();
        int nDisplayableClasses = displayableClasses.size();
        for(int n = 0; n < nDisplayableClasses; n++)
        {
            StsClassCursorDisplayable displayableClass = (StsClassCursorDisplayable) displayableClasses.get(n);
            displayableClass.drawOn3dCurtain(glPanel3d, seismicCurtain);
        }
    }

    protected void displayTexture()
    {
        if(seismicVolumeChanged())
        {
            if(seismicVolume == null) return;
            deleteTexturesAndDisplayLists(gl);
        }

        if(textureTiles == null)
        {
            if(!initializeTextureTiles()) return;
            textureChanged = true;

        }
        if(isPixelMode != seismicClass.getIsPixelMode())
        {
            deleteTexturesAndDisplayLists(gl);
            isPixelMode = !isPixelMode;
            textureChanged = true;
        }
        byte projectZDomain = model.getProject().getZDomain();
        if(projectZDomain != seismicCurtain.zDomain)
        {
            deleteDisplayLists(gl);
            zDomain = projectZDomain;
            usingDisplayLists = false;
        }
        useDisplayLists = model.useDisplayLists;
        if(!useDisplayLists && usingDisplayLists)
        {
            deleteDisplayLists(gl);
            usingDisplayLists = false;
        }

        //computeProjectionMatrix(gl, glu); // should call only if view changed
        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);

        if(isPixelMode != getIsPixelMode())
        {
            textureChanged = true;
            isPixelMode = !isPixelMode;
        }

        setGLColorList(gl);

        byte[] data = null;
        if(textureChanged)
        {
            data = getData();
        }

        if(textureChanged || useDisplayLists && !usingDisplayLists)
        {
            if(textureChanged) textureChanged = false;
            else if(useDisplayLists) usingDisplayLists = true;
            if(textureTiles == null)
                StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
            textureTiles.constructSurface(this, gl, useDisplayLists, false);
        }
        textureTiles.displayTiles2d(this, gl, false, isPixelMode, data, StsParameters.nullByte);

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

        if(textureTiles.shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
    }

    private boolean seismicVolumeChanged()
    {
        StsSeismicVolume currentSeismicVolume = seismicClass.getCurrentSeismicVolume();
        if(seismicVolume == currentSeismicVolume) return false;
        if(currentSeismicVolume == null) return true;
        seismicVolume.removeActionListener(this);
        seismicVolume = currentSeismicVolume;
        seismicVolume.addActionListener(this);
        return true;
    }

    private boolean initializeTextureTiles()
    {

        float[][] textureRanges;
        textureRanges = totalAxisRanges;
        // if(textureTiles != null) deleteTextureAndSurface(gl);
        textureTiles = StsTextureTiles.constructor(model, this, nTraces, nSamples, true, textureRanges);
        if(textureTiles == null) return false;
        textureChanged = true;
        return true;
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        if(seismicCurtain.zDomain == seismicVolume.getZDomain())
            drawTextureTileTimeSurface2d(tile, gl);
        else
        {
            StsSeismicVelocityModel velocityVolume = model.getProject().velocityModel;
            if(velocityVolume == null) return;
            drawTextureTileDepthSurface2d(velocityVolume, tile, gl);
        }
    }

    private void drawTextureTileTimeSurface2d(StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        StsGridPoint[] arcGridPoints = seismicCurtain.cellGridCrossingPoints;
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
        double dRowTexCoor = tile.dRowTexCoor;
        float tileZMin = seismicCurtain.tMin + tile.colMin * seismicCurtain.tInc;
        float tileZMax = seismicCurtain.tMin + tile.colMax * seismicCurtain.tInc;
        gl.glBegin(GL.GL_QUAD_STRIP);
        for(int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float arcLength = arcGridPoints[row].point.getF();;
            gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
            gl.glVertex2f(arcLength, tileZMin);
            gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
            gl.glVertex2f(arcLength, tileZMax);
        }
        gl.glEnd();
    }

    public void drawTextureTileDepthSurface2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        if(velocityModel == null) return;

        if(velocityModel.getInputVelocityVolume() != null)
            drawTextureTileDepthSurfaceFromVolume2d(velocityModel, tile, gl);
        else
            drawTextureTileDepthSurfaceFromIntervalVelocities2d(velocityModel, tile, gl);
    }

    public void drawTextureTileDepthSurfaceFromVolume2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        StsGridPoint[] arcGridPoints = seismicCurtain.cellGridCrossingPoints;
        float tt = zCoordinates[tile.rowMin];
        float[] xyz = seismicCurtain.cellGridCrossingPoints[tile.rowMin].getXYZorT();
        float x1 = xyz[0];
        float y1 = xyz[1];
        StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
        int r1 = velocityVolume.getNearestBoundedRowCoor(y1);
        int c1 = velocityVolume.getNearestBoundedColCoor(x1);
        float velocityTInc = velocityVolume.zInc;
        float[] velocityTrace1 = velocityVolume.getTraceValues(r1, c1);
        float depthDatum = velocityModel.depthDatum;
        float timeDatum = velocityModel.timeDatum;
        float tMin = seismicCurtain.tMin;
        float tInc = seismicCurtain.tInc;
        float arc1 = arcGridPoints[tile.rowMin].point.getF();
        for(int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float arc0 = arc1;
            arc1 = arcGridPoints[row].point.getF();
            float[] velocityTrace0 = velocityTrace1;

            xyz = seismicCurtain.cellGridCrossingPoints[row].getXYZorT();
            x1 = xyz[0];
            y1 = xyz[1];
            r1 = velocityVolume.getNearestBoundedRowCoor(y1);
            c1 = velocityVolume.getNearestBoundedColCoor(x1);
            velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            float t = tMin + tile.colMin * tInc;
            for(int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
            {
                float v0 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace0);
                float z0 = (v0 * (t - timeDatum) + depthDatum);
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex2f(arc0, z0);
                float v1 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace1);
                float z1 = (v1 * (t - timeDatum) + depthDatum);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex2f(arc1, z1);
            }
            gl.glEnd();
        }
    }

    public void drawTextureTileDepthSurfaceFromIntervalVelocities2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        StsPoint point = seismicCurtain.cellGridCrossingPoints[tile.rowMin].getPoint();
        float x1 = point.getX();
        float y1 = point.getY();
        float tMin = seismicCurtain.tMin;
        float tInc = seismicCurtain.tInc;
        float displayTMin = tMin + tile.colMin * tInc;
        StsGridPoint[] arcGridPoints = seismicCurtain.cellGridCrossingPoints;
        float arc1 = arcGridPoints[tile.rowMin].point.getF();
        for(int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float arc0 = arc1;
            arc1 = arcGridPoints[row].point.getF();

            float x0 = x1;
            float y0 = y1;
            point = seismicCurtain.cellGridCrossingPoints[row].getPoint();
            x1 = point.getX();
            y1 = point.getY();

            int nSlices = tile.colMax - tile.colMin + 1;
            float[] depths0 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x0, y0);
            float[] depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x1, y1);

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            for(int n = 0, col = tile.colMin; col <= tile.colMax; n++, col++, colTexCoor += dColTexCoor)
            {
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex2f(arc0, depths0[n]);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex2f(arc1, depths1[n]);
            }
            gl.glEnd();
        }
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void geometryChanged()
    {
    }

    private void deleteDisplayLists(GL gl)
    {
        if(textureTiles != null)
            textureTiles.deleteDisplayLists(gl);
    }

    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if(textureTiles == null) return;
        if(debug) StsException.systemDebug(this, "deleteTextureTileSurface");
        textureTiles.deleteTextures(gl);
        textureTiles.deleteDisplayLists(gl);
        textureChanged = true;
    }

    protected void initializeView()
    {
        gl.glShadeModel(GL.GL_FLAT);
        gl.glEnable(GL.GL_ALPHA_TEST);
        gl.glAlphaFunc(GL.GL_GREATER, 0.1f);
        initializeFont(gl);
        initializeInsets();
        computePixelScaling();
        computeProjectionMatrix();
        computeModelViewMatrix();
    }

    /** Set the aspect ratio to 1:1. */
    public void aspectRatio()
    {
        if(getCursor3d().selectedDirection != getCursor3d().ZDIR) return;

        float halfXRange = Math.abs(axisRanges[0][1] - axisRanges[0][0]) / 2.0f;
        float yPos = axisRanges[1][0] + Math.abs(axisRanges[1][1] - axisRanges[1][0]) / 2.0f;
        axisRanges[1][0] = yPos - halfXRange;
        axisRanges[1][1] = yPos + halfXRange;
        computePixelScaling();
        viewChangedRepaint();
    }

    /**
     * Display the popup menu)
     *
     * @param mouse
     */
    public void showPopupMenu(StsMouse mouse)
    {
        // chec if well was picked
        JPopupMenu tp = new JPopupMenu("2D View Popup");
        glPanel3d.add(tp);

        tp.add(aspectBtn);
        tp.show(glPanel3d, mouse.getX(), mouse.getY());
        clearCurrentKeyPressed();
    }

    public void cancelPopupMenu(StsMouse mouse)
    {
        if(tp != null) tp.setVisible(false);
    }

    private float[][] resetDisplayRanges(float[][] axisRanges, float[][] totalAxisRanges, int selectedDirection)
    {
        float[][] axisDisplayRanges = new float[2][2];
        float xMin = totalAxisRanges[0][0];
        axisDisplayRanges[0][0] = axisRanges[0][0] - xMin;
        axisDisplayRanges[0][1] = axisRanges[0][1] - xMin;
        if(selectedDirection != ZDIR)
        {
            axisDisplayRanges[1][0] = axisRanges[1][0];
            axisDisplayRanges[1][1] = axisRanges[1][1];
        }
        else
        {
            float yMin = totalAxisRanges[1][0];
            axisDisplayRanges[1][0] = axisRanges[1][0] - yMin;
            axisDisplayRanges[1][1] = axisRanges[1][1] - yMin;
        }
        return axisDisplayRanges;
    }
/*
    public void displayForeground()
    {
       	if(cursorButtonState == StsMouse.PRESSED)
        {
//            System.out.println("Pressed");
            previousCursorPoint = null;
            currentCursorPoint = cursorPoint;
            drawForeground();
        }
        else if(leftButtonState == StsMouse.DRAGGED)
        {
//            System.out.println("Dragged");
            currentCursorPoint = cursorPoint;
            drawForeground();
        }
        else // StsMouse.RELEASED
        {
//            System.out.println("Released");
            drawForeground();
            cursorButtonState = StsMouse.CLEARED;
        }
    }
*/

    /**
     * Draws in foreground XORed against current view.  Must be called
     * again and draw same objects to erase.
     */
    public void drawForeground()
    {
        gl.glDisable(GL.GL_BLEND);
        gl.glDrawBuffer(GL.GL_FRONT);
        gl.glEnable(GL.GL_COLOR_LOGIC_OP);
        gl.glLogicOp(GL.GL_XOR);
        gl.glDepthMask(false);
        gl.glDepthFunc(GL.GL_ALWAYS);
        drawForegroundCursor();
        gl.glFlush();
        if(cursorButtonState == StsMouse.RELEASED)
        {
            gl.glLogicOp(GL.GL_COPY);
            gl.glDisable(GL.GL_COLOR_LOGIC_OP);
            gl.glDepthMask(true);
            gl.glDepthFunc(GL.GL_LESS);
            gl.glDrawBuffer(GL.GL_BACK);
            cursorButtonState = StsMouse.CLEARED;
        }
    }

    /**
     * Draw the cursor in the front buffer.  If previously drawn,
     * draw previous one again to erase it; then drawn new one.
     */
    private void drawForegroundCursor()
    {
        if(previousXYZ != null)
            drawCursorPoint(previousXYZ, this.glPanel3d);
        if(currentXYZ != null)
            drawCursorPoint(currentXYZ, glPanel3d);
        previousXYZ = currentXYZ;
        currentXYZ = null;
    }

    /**
     * Draw vertical and horizontal lines thru the cursor point.
     * because the crossed area is drawn twice it will be clear.
     */
    private void drawCursorPoint(float[] xyz, StsGLPanel3d glPanel3d)
    {
        switch(getCursor3d().selectedDirection)
        {
            case XDIR:
                StsGLDraw.drawLine2d(xyz[1], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                StsGLDraw.drawLine2d(xyz[1], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                break;
            case YDIR:
                StsGLDraw.drawLine2d(xyz[0], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                StsGLDraw.drawLine2d(xyz[0], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                break;
            case ZDIR:
                if(!axesFlipped)
                {
                    StsGLDraw.drawLine2d(xyz[0], xyz[1], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                    StsGLDraw.drawLine2d(xyz[0], xyz[1], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                }
                else
                {
                    StsGLDraw.drawLine2d(xyz[1], xyz[0], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                    StsGLDraw.drawLine2d(xyz[1], xyz[0], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                }
        }
    }

    private void setLabels()
    {
        StsCursor3d cursor3d = getCursor3d();

        int selectedDirection = cursor3d.selectedDirection;
        if(!axesFlipped)
        {
            axisLabels[0] = cursor3d.getHorizontalAxisLabel(selectedDirection);
            axisLabels[1] = cursor3d.getVerticalAxisLabel(selectedDirection);
        }
        else
        {
            axisLabels[1] = cursor3d.getHorizontalAxisLabel(selectedDirection);
            axisLabels[0] = cursor3d.getVerticalAxisLabel(selectedDirection);
        }
    }

    public boolean initializeDefaultAction()
    {
        return setDefaultAction3d();
    }

    public boolean setDefaultAction3d()
    {
        StsActionManager actionManager = glPanel3d.window.actionManager;
        if(actionManager == null) return true;
        StsAction defaultAction = getDefaultAction();
        if(defaultAction instanceof StsDefaultAction3d) return false;
        setDefaultAction(new StsDefaultAction3d(glPanel3d, actionManager));
        return true;
    }

    /** StsView2d */
    public boolean viewObjectChanged(Object source, Object object)
    {
        if(!(object instanceof StsSeismicVolume)) return false;
        if(object == seismicVolume) return false;
        seismicVolume = (StsSeismicVolume) object;
        textureChanged();
        return true;
    }

    /**
     * object being viewed is changed. Repaint this view if affected.
     * Implement as needed in concrete subclasses.
     */
    public boolean viewObjectRepaint(Object source, Object object)
    {
        if(!(object instanceof StsSeismicVolume)) return false;
        if(object == seismicVolume) return false;
        viewChangedRepaint();
        return true;
    }


    public void actionPerformed(ActionEvent actionEvent)
    {
        textureChanged = true;
        repaint();
    }

    /** view2d */
    public void doInitialize()
    {

    }

    public void setAxisRanges()
    {

    }

    public void setDefaultView()
    {

    }

    /*
        public void init()
        {
            if(isGLInitialized)return;
            initGL();
            gl.glShadeModel(GL.GL_FLAT);
            gl.glEnable(GL.GL_ALPHA_TEST);
            gl.glAlphaFunc(GL.GL_GREATER, 0.1f);
            initializeFont(gl);
            isGLInitialized = true;
        }
    */
    public void setInsets(boolean axisOn)
    {
        if(axisOn)
        {
            int leftInset = halfWidth + majorTickLength + 2 * verticalFontDimension.width + 2 * fontLineSpace;
            ;
            int bottomInset = 0;
            int topInset = halfWidth + majorTickLength + 4 * horizontalFontDimension.height + 3 * fontLineSpace;
            int rightInset = 0;
            insets = new Insets(topInset, leftInset, bottomInset, rightInset);
        }
        else
        {
            insets = new Insets(0, 0, 0, 0);
        }
    }

    public void reshape(GLAutoDrawable component, int x, int y, int width, int height)
    {

    }

    protected StsColor getGridColor()
    {
        return StsColor.BLACK;
    }

    public void setDisplayPanelSize(Dimension size)
    {
        glPanel.setPreferredSize(size);
    }

    public void viewChangedRepaint()
    {
        glPanel.repaint();
    }

    public StsSeismicCurtain getSeismicCurtain()
    {
        return seismicCurtain;
    }

    public void resetToOrigin()
    {
    }

    public byte getHorizontalAxisType()
    {
        return AXIS_TYPE_NONE;
    }

    public byte getVerticalAxisType()
    {
        return AXIS_TYPE_NONE;
    }

    public boolean getUseShader()
    {
        return seismicClass.getContourColors();
    }

    public int getDefaultShader()
    {
        return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
    }

    protected boolean getIsPixelMode()
    {
        return seismicVolume.getIsPixelMode();
    }

    protected void setGLColorList(GL gl)
    {
        seismicVolume.seismicColorList.setGLColorList(gl, false, StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS);
    }

    protected void clearShader()
    {
        seismicVolume.seismicColorList.setGLColorList(gl, false, StsJOGLShader.NONE);
    }

    protected byte[] getData()
    {
        return seismicCurtain.getData2D(seismicVolume);
    }

    protected float compute2dValue(GL gl, double x, double y)
    {
        return StsParameters.nullValue;
    }

    protected String getValueLabel(double xCoordinate)
    {
        return "";
    }

    public boolean textureChanged()
    {
        if(seismicVolume == null)
        {
            textureChanged = false;
            return false;
        }
        return textureChanged;
    }

    public Class getDisplayableClass()
    {
        return StsSeismicVolume.class;
    }
}



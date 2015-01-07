package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Processes the crossplot data solely for mainDebug purposes</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.nio.*;

public class StsViewXPData extends StsView2d
{
    StsCrossplot crossplot = null;
    byte[] xData, yData;
    StsSpectrum spectrum;
    int colorListNum;

    FloatBuffer data = null;

    static final int width = 256;
    static final int height = 256;

    static StsViewXPData viewXPData = null;
    /** Name of this object XPData */
//    public final String name = "XPData";

	static final String viewNameXPData = "XP View Data";
    static final String shortViewNameXPData = "XPD";
    /**
     * StsViewXPData constuctor
     * @param glPanel3d the graphics context
     * @param crossplot the crossplot object and associated data
     */
    StsViewXPData(StsGLPanel3d glPanel3d, StsCrossplot crossplot)
    {
        super(glPanel3d);
        this.crossplot = crossplot;
        initialize();
    }

    /**
     * Create an instance of this object, set the current view to the object
     * @param model the current model
     * @param crossplot the crossplot object to present
     * @return this object
     */
    static public StsViewXPData setView(StsModel model, StsCrossplot crossplot)
    {
        if(viewXPData == null) viewXPData = new StsViewXPData(model.getGlPanel3d(), crossplot);
        model.getGlPanel3d().checkAddView(StsViewXPData.class);
        return viewXPData;
    }

    static public String getStaticViewName()
    {
        return viewNameXPData;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameXPData;
    }

    public void initialize()  {}

	public void setAxisRanges()
	{
		totalAxisRanges = crossplot.getTotalAxisRanges();
	}
    /**
     * Load the local data buffer
     */
    public void copyToBuffer( GL gl )
    {
        data = BufferUtil.newFloatBuffer(4*width*height);
        gl.glReadBuffer(GL.GL_FRONT);
        gl.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_FLOAT, data);
    }
    /**
     * Get the local data buffer
     * @return data
     */
    public float[] getData() { return data.array(); }
    /**
     * Reset the view back to the crossplot view, disabling the crossplot data view
     */
    public void unsetView()
    {
        glPanel3d.checkAddView(StsViewXP.class);
    }

    public void computeProjectionMatrix()
    {
        if(totalAxisRanges == null) return;
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(totalAxisRanges[0][0], totalAxisRanges[0][1], totalAxisRanges[1][0], totalAxisRanges[1][1]);
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

	public boolean viewObjectChanged(Object source, Object object) { return false; }
    public boolean viewObjectRepaint(Object source, Object object) { return false; }

    /**
     * Is the crossplot data viewable
     * @return true is viewable
     */
    public boolean isViewable() { return crossplot != null; }
    public void reshape( GLAutoDrawable drawable, int x, int y, int w, int h ) {}
    /**
     * Display method for crossplot data view
     * @param component drawable component
     */
    public void display(GLAutoDrawable component)
    {
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }  
        if(crossplot == null || totalAxisRanges == null) return;

//        glPanel3d.setViewPort(0, 0, width, height);

        if(glPanel3d.viewChanged)
        {
            computeProjectionMatrix();
            computeModelViewMatrix();
            glPanel3d.viewChanged = false;
        }

        gl.glDrawBuffer(GL.GL_BACK);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT|GL.GL_COLOR_BUFFER_BIT);
        crossplot.drawPolygons(gl, glu, true);
//        glPanel3d.swapBuffers();
        copyToBuffer( gl );
        return;
    }
    /**
     * Key Release event handling
     *
     *
     * @param mouse mouse object
     * @param e key event
     * @return true if successful
     */
    public void keyReleased(StsMouse mouse, KeyEvent e)
    {
        int releasedKeyCode = e.getKeyCode();
        super.keyReleased(mouse, e);

        if(releasedKeyCode == KeyEvent.VK_UP)
        {
            model.getGlPanel3d().checkAddView(StsView3d.class);
        }
        else if (releasedKeyCode == KeyEvent.VK_DOWN)
        {
            model.getGlPanel3d().checkAddView(StsViewCursor.class);
        }
        else if(releasedKeyCode == KeyEvent.VK_LEFT)
        {
            StsCrossplotClass crossplotClass = (StsCrossplotClass)StsCrossplot.getCrossplotClass();
            crossplotClass.previousCrossplot();
 //           model.win3d.getViewSelectToolbar().setButtonEnabled(model.glPanel3d.getCurrentView(), true);
        }
        else if(releasedKeyCode == KeyEvent.VK_RIGHT)
        {
            StsCrossplotClass crossplotClass = (StsCrossplotClass)StsCrossplot.getCrossplotClass();
            crossplotClass.nextCrossplot();
 //           model.win3d.getViewSelectToolbar().setButtonEnabled(model.glPanel3d.getCurrentView(), true);
        }
    }
    /** for abstract class compatibility */
    public void setDefaultView() { }


	/*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;

	public void resurrect(StsView o)
	{
	}

    public void resetToOrigin() { }

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_NONE; }
}

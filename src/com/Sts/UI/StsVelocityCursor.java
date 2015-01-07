package com.Sts.UI;

/**
 * StsVelocityCursor
 * 
 * Class for storing and drawing a cursor that displays current time/velocity in all velocity-picking related windows.
 * 
 * Singleton patter ensures that all displays will show the same time/velocity pair at all times
 * MouseMotionListener interface is used to get mouseMoved events from display panels.
 * X'es are draw in Stack and Semblance displays
 * Red NMO curve for current time/velocity is displayed on Gather display
 * Time/Vel displayed in lower-left corner of Semblance display as ascii text
 * 
 * To use:
 * Views wanting to add this cursor need to call:
 * createVelocityCursorMouseMotionListener()
 * upon initialization
 * 
 * then add:
 * displayCursor(GL gl, StsViewPreStack view)
 * to the views display() method.
 * 
 * @author Scott Cook - August 20, 2009
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.util.*;

public class StsVelocityCursor
{
    public static final int CursorSize = 10;
    private static StsVelocityCursor stsVelocityCursor = null;
    protected double time;
    protected double vel;
    protected ArrayList<StsView> views;
    
    private StsVelocityCursor()
    {
        views = new ArrayList<StsView>();
    }

    
    public void addView(StsView view)
    {
        views.add(view);
    }

    // using singleton pattern. There can only be one stsVelocityCursor in existence at a time per JVM
    // This way, all views will show the same time/velocity pairs to aid in picking
    public static synchronized StsVelocityCursor getStsVelocityCursor()
    {
        if (stsVelocityCursor == null) stsVelocityCursor = new StsVelocityCursor();
        return stsVelocityCursor;
    }
    
    public static MouseMotionListener createVelocityCursorMouseMotionListener(StsViewPreStack view)
    {
        MouseMotionListener listener = getStsVelocityCursor().newListener(view);
        getStsVelocityCursor().addView(view);
        return listener;
    }
    
    private MouseMotionListener newListener(StsViewPreStack view)
    {
        return new VelocityCursorMouseMotionListener(view);
    }

    private class VelocityCursorMouseMotionListener implements MouseMotionListener
    {
        StsViewPreStack view;
        
        public VelocityCursorMouseMotionListener(StsViewPreStack view)
        {
            this.view = view;
        }

        public void mouseDragged(MouseEvent e)
        {
            mouseMoved(e);
        }

        public void mouseMoved(MouseEvent e)
        {
            StsMouse mouse = StsGLPanel.mouse;
            mouse.setMouseCoordinates(e, e.getID());
            StsPoint point = view.computePickPoint(mouse);
            if (point != null)
            {
                time = point.getY();
                time = StsMath.minMax(time, view.lineSet.zMin, view.lineSet.zMax);
                vel = point.getX();
            }
            //System.out.println("time: "+time+" vel: "+vel);
            for (int i=0; i<views.size(); i++)
            {
                views.get(i).repaint();
            }
        }
        

    }

    public void displayCursor(GL gl, StsViewPreStack view)
    {
        if (view == null) return;
        view.computeHorizontalPixelScaling();
        view.computeVerticalPixelScaling();
        if (view instanceof StsViewSemblance)
        {
            drawX(view, vel, time, CursorSize);
            float[][] ranges = view.axisRanges;
            float xMin = ranges[0][0];
            float yMin = ranges[1][0];
            float shiftX = 5/view.pixelsPerXunit;
            float shiftY = 5/view.pixelsPerYunit;
            drawLabel(gl, view, StsColor.WHITE, xMin+shiftX, yMin+shiftY);
            drawLabel(gl, view, StsColor.BLACK, (float)vel+0.3f, (float)time);
        }
        else if (view instanceof StsViewCVStacks)
        {
            drawX(view, vel, time, CursorSize);
            float shiftX = 12/view.pixelsPerXunit;
            drawLabel(gl, view, StsColor.BLACK, (float)vel+shiftX, (float)time);
        }
        else if (view instanceof StsViewVVStacks)
        {
            StsVelocityProfile currentVelocityProfile = view.superGather.velocityProfile;
            if (currentVelocityProfile == null) return;
            StsPoint[] vvsInitialProfilePoints = currentVelocityProfile.getSetVvsInitialProfilePoints();
            if (vvsInitialProfilePoints == null || vvsInitialProfilePoints.length == 0) return;
            float velocity = 0;
            if (time < vvsInitialProfilePoints[0].v[1]) 
                velocity = vvsInitialProfilePoints[0].v[0];
            else
                velocity = StsMath.interpolateValue(vvsInitialProfilePoints, (float)time, 1, 0);
            double percent = (vel/velocity - 1)/100;
            drawX(view, percent*10000, time, CursorSize);
            float shiftX = 12/view.pixelsPerXunit;
            drawLabel(gl, view, StsColor.BLACK, (float)percent*10000+shiftX, (float)time);
        }
        else if (view instanceof StsViewGather)
        {
            drawNMOCurve(gl,view);
        }
    }

    private void drawNMOCurve(GL gl, StsViewPreStack view)
    {
        StsSuperGather gather = view.superGather;
        if (gather == null) return;
        if (!gather.checkInitializeGather()) return;
        if (!gather.checkComputeVelocitiesAndOffsetTimes()) return;
        StsPoint point = new StsPoint((float)vel, (float)time);
        StsVelocityProfile velocityProfile = view.superGather.velocityProfile;
        if (velocityProfile == null) velocityProfile = view.getVelocityProfile();
        if (velocityProfile == null) return;
        StsWiggleDisplayProperties wiggleProperties = view.lineSet.getWiggleDisplayProperties();
        if (wiggleProperties == null) return;
        int offsetAxisType = wiggleProperties.getOffsetAxisType();
        gl.glLineWidth(2);
        gather.displayPickNMOPoint(point, velocityProfile, StsGather.NONE_INDEX, offsetAxisType, false, false, gl, StsColor.FORESTGREEN);
    }

    /**
     * draw label showing velocity and time of current mouse position.
     * label drawn at position x,y
     * @param gl
     * @param view
     * @param stsColor 
     * @param x
     * @param y
     */
    private void drawLabel(GL gl, StsViewPreStack view, StsColor stsColor, float x, float y)
    {
        String text = (int)Math.round(time) + ", " + (int)Math.round(vel*1000);
        stsColor.setGLColor(gl);
        StsGLDraw.fontOutput(gl, x, y, text, com.magician.fonts.GLHelvetica18BitmapFont.getInstance(gl)); //larger text - easier to read
        gl.glDrawBuffer(GL.GL_FRONT_AND_BACK);
    }               


    private void drawX(StsViewPreStack view, double x, double y, int size)
    {
        double width = size/view.pixelsPerXunit;
        double height = size/view.pixelsPerYunit;
        GL gl = view.glPanel3d.getGL();
        gl.glLineWidth(3);
        StsColor.BLACK.setGLColor(gl);
        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_LINE_STRIP); 
            gl.glVertex2d(x - width, y - height);
            gl.glVertex2d(x + width, y + height);
            gl.glEnd();
            gl.glBegin(GL.GL_LINE_STRIP); 
            gl.glVertex2d(x - width, y + height);
            gl.glVertex2d(x + width, y - height);
        }
        catch(Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glDisable(GL.GL_LINE_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }
    
    public double getTime()
    {
        return time;
    }


    public double getVel()
    {
        return vel;
    }

}
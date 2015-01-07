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
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.util.*;

public class StsGatherCursor
{
    public static final int CursorSize = 10;
    private static StsGatherCursor stsGatherCursor = null;
    protected double time;
    protected double x;
    protected ArrayList<StsView> views;
    
    private StsGatherCursor()
    {
        views = new ArrayList<StsView>();
    }

    
    public void addView(StsView view)
    {
        views.add(view);
    }

    // using singleton pattern. There can only be one stsGatherCursor in existence at a time per JVM
    // This way, all views will show the same time/velocity pairs to aid in picking
    public static synchronized StsGatherCursor getStsGatherCursor()
    {
        if (stsGatherCursor == null) stsGatherCursor = new StsGatherCursor();
        return stsGatherCursor;
    }
    
    public static MouseMotionListener createGatherCursorMouseMotionListener(StsViewPreStack view)
    {
        MouseMotionListener listener = getStsGatherCursor().newListener(view);
        getStsGatherCursor().addView(view);
        return listener;
    }
    
    private MouseMotionListener newListener(StsViewPreStack view)
    {
        return new GatherCursorMouseMotionListener(view);
    }

    private class GatherCursorMouseMotionListener implements MouseMotionListener
    {
        StsViewPreStack view;
        
        public GatherCursorMouseMotionListener(StsViewPreStack view)
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
                x = point.getX();
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
        if (view instanceof StsViewGather)
        {
            drawX(view, x, time, CursorSize);
            float shift = CursorSize/view.pixelsPerXunit;
            drawLabel(gl, view, StsColor.RED, (float)x+shift, (float)time);
        }
        
    }

    private void drawLabel(GL gl, StsViewPreStack view, StsColor stsColor, float x, float y)
    {
        double val = ((StsViewGather)view).getCursorValue(x);
        String text = (int)Math.round(time) + ", " + (int)Math.round(val);
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


    public double getX()
    {
        return x;
    }

}
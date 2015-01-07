package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.GL;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 8, 2008
 * Time: 3:29:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsLineMeasure extends StsLineSegment
{
    private String label = null;
    private byte placement = START;
    private boolean showLabel = true;
    private boolean addEnds = false;
	public float radius = 50.0f;
	public static final byte START = 0;
    public static final byte END = 1;
    public static final byte MIDDLE = 2;

    public StsLineMeasure()
    {
        super();
    }

    public StsLineMeasure(StsColor color)
    {
        super(color);
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        if((getLastPoint() == null) || (getFirstPoint() == null))
            return;

        super.display(glPanel3d);

        // Add graphical output
        if(!showLabel)
            return;
        GL gl = glPanel3d.getGL();
        if(getLabel() != null)
            StsGLDraw.fontHelvetica12(gl, getTextPosition(), getLabel());

        if(addEnds)
            drawEndPoints(glPanel3d);
    }

    public void drawEndPoints(StsGLPanel glPanel)
    {
        StsGLDraw.drawSphere((StsGLPanel3d)glPanel, getFirstPoint().getXYZ(),StsColor.ALICEBLUE, radius);
        StsGLDraw.drawSphere((StsGLPanel3d)glPanel, getLastPoint().getXYZ(),StsColor.ALICEBLUE, radius);
    }

    public void pickPoint(StsGLPanel3d glPanel)
    {
        GL gl = glPanel.getGL();
        gl.glInitNames();

            gl.glPushName(START);
            StsGLDraw.drawSphere(glPanel, getFirstPoint().getXYZ(),StsColor.ALICEBLUE, radius);
            gl.glPopName();

            gl.glPushName(END);
            StsGLDraw.drawSphere(glPanel, getLastPoint().getXYZ(),StsColor.ALICEBLUE, radius);
            gl.glPopName();
    }

    public float[] getTextPosition()
    {
       switch(getPlacement())
       {
           case START:
              return getFirstPoint().getXYZorT();
           case END:
              return getLastPoint().getXYZorT();
           case MIDDLE:
              float[] xyz = new float[3];
              xyz[0] = getFirstPoint().getX() + Math.abs((getFirstPoint().getX() - getLastPoint().getX()) / 2.0f);
              xyz[1] = getFirstPoint().getY() + Math.abs((getFirstPoint().getY() - getLastPoint().getY()) / 2.0f);
              xyz[2] = getFirstPoint().getZorT() + Math.abs((getFirstPoint().getZorT() - getLastPoint().getZorT()) / 2.0f);
              return xyz;
           default:
               return getLastPoint().getXYZorT();
       }
    }
    public StsPoint getPoint(int index)
    {
        if(index == START)
            return getFirstPoint();
        else if(index == END)
            return getLastPoint();
        else
            return null;
    }
    public String getLabel() { return label; }
    public void setLabel(String lbl) { label = lbl; }
    public byte getPlacement() { return placement; }
    public void setPlacement(byte place) { placement = place; }
    public boolean getAddEnds() { return addEnds; }
    public void setAddEnds(boolean ends) { addEnds = ends; }
    public boolean getShowLabel() { return showLabel; }
    public void setShowLabel(boolean show) { showLabel = show; }
}
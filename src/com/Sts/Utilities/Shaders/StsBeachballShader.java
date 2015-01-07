package com.Sts.Utilities.Shaders;

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 14, 2009
 * Time: 10:54:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsBeachballShader extends StsShader
{
    boolean drawBeachball = true;

    StsUniformVariableF darkMaterialColor = new StsUniformVariableF("DarkMaterialColor", new float[] { 0.0f, 0.6f, 0.0f, 1.0f } );
    StsUniformVariableF lightMaterialColor = new StsUniformVariableF("LightMaterialColor", new float[] { 0.5f, 0.5f, 0.5f, 1.0f } );
    StsUniformVariableInt drawType = new StsUniformVariableInt("DrawType", new int[] { 0 } );

    public StsBeachballShader(GL gl) throws GLException
    {
        super(gl, "Beachball");
    }

    static public StsBeachballShader getShader(GL gl)
    {
        return (StsBeachballShader)getEnableShader(gl, StsBeachballShader.class);
    }

    public void display(StsGLPanel3d glPanel3d, float[] xyz, float[] dsr)
    {
        drawBeachball(glPanel3d, xyz, dsr, StsColor.RED, 10.0f);
    }

    public void setBeachballColors(GL gl, StsColor color)
    {
        darkMaterialColor.setValues(gl, shaderProgram, color.getRGBA(), 4);
        lightMaterialColor.setValues(gl, shaderProgram);
    }

    public void setdrawType(GL gl, int type)
    {
        drawType.setValues(gl, shaderProgram, new int[] { type }, 1);
    }

    public void drawBeachball(StsGLPanel3d glPanel3d, float[] xyz, float[] dsr, StsColor color, float size)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        try
        {
            int displayListSphere = StsGLDraw.getDisplayListSphere(glPanel3d);
            if(displayListSphere == 0) return;
            
            float zscale = glPanel3d.getZScale();
            setPosition(gl, xyz, dsr, size, zscale);

            StsBeachballShader.getShader(gl);
            setBeachballColors(gl, color);
            gl.glCallList(displayListSphere);
            disableCurrentShader(gl);
            drawBeachballAxes(gl);
        }
        finally
        {
            gl.glPopMatrix();
        }
    }

     private void setPosition(GL gl, float[] xyz, float[] dsr, float size, float zscale)
     {
        gl.glEnable(GL.GL_LIGHTING);
        gl.glShadeModel(GL.GL_SMOOTH);
        gl.glPushMatrix();
        gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
        float sphereSize = StsGLDraw.sphereSize;
        float scaleFactor = size/sphereSize;
        gl.glScalef(scaleFactor, scaleFactor, scaleFactor / zscale);
        gl.glRotatef(90.0f-dsr[1], 0.0f, 0.0f, 1.0f);
        gl.glRotatef(90.0f-dsr[0], 1.0f, 0.0f, 0.0f);
        gl.glRotatef(-dsr[2], 0.0f, 1.0f, 0.0f);
     }

    static float[] beachballCenter = new float[] { 0.0f, 0.0f, 0.0f };
    static float[] beachballXPoint = new float[] { StsGLDraw.sphereSize*1.25f, 0.0f, 0.0f };
    static float[] beachballYPoint = new float[] { 0.0f, StsGLDraw.sphereSize*1.25f, 0.0f };
    static float[] beachballZPoint = new float[] { 0.0f, 0.0f, -StsGLDraw.sphereSize*1.25f };

    public void drawBeachballAxes(GL gl)
    {
        try
        {
            // color.setGLColor(gl);
            gl.glDisable(GL.GL_LIGHTING);
            StsGLDraw.drawLineSegment(gl, StsColor.RED, beachballCenter, beachballXPoint, 3.0f);
            StsGLDraw.drawLineSegment(gl, StsColor.GREEN, beachballCenter, beachballYPoint, 3.0f);
            StsGLDraw.drawLineSegment(gl, StsColor.BLUE, beachballCenter, beachballZPoint, 3.0f);
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);
        }
    }
}

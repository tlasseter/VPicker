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
public class StsSphereShader extends StsShader
{
    StsUniformVariableF materialColor = new StsUniformVariableF("materialColor", new float[] { 0.0f, 0.6f, 0.0f, 1.0f } );


    public StsSphereShader(GL gl) throws GLException
    {
        super(gl, "Sphere");
    }

    static public StsSphereShader getShader(GL gl)
    {
        return (StsSphereShader)getEnableShader(gl, StsSphereShader.class);
    }

    public void display(StsGLPanel3d glPanel3d, float[] xyz, float[] dsr)
    {
        StsGLDraw.drawSphere(glPanel3d, xyz, StsColor.RED, 10.0f);
    }

    public void setColor(GL gl, StsColor color)
    {
        materialColor.setValues(gl, shaderProgram, color.getRGBA(), 4);
    }

    public void drawSphere(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        try
        {
            getEnableShader(gl, StsSphereShader.class);
            int displayListSphere = StsGLDraw.getDisplayListSphere(glPanel3d);
            if(displayListSphere == 0) return;
            //gl.glEnable(GL.GL_LIGHTING);
            //gl.glShadeModel(GL.GL_SMOOTH);
            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);

            float zscale = glPanel3d.getZScale();
            float sphereSize = StsGLDraw.sphereSize;
            float scaleFactor = size/sphereSize;
            gl.glScalef(scaleFactor, scaleFactor, scaleFactor / zscale);

            setColor(gl, color);
            gl.glCallList(displayListSphere);
        }
        finally
        {
            gl.glPopMatrix();
        }
    }
}

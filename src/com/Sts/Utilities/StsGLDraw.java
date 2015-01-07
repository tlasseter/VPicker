//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.magician.fonts.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.nio.*;
import java.util.*;

public class StsGLDraw
{
    static int displayListSphere = 0;
    static public float sphereSize = 100.0f;

    static final float nullValue = StsParameters.nullValue;
    static float DEG2RAD = 3.14159f / 180;

	static public void initialize()
	{
		displayListSphere = 0;
    }

    static public void drawFilledEllipse(float[] xyz, StsColor color, StsGLPanel glPanel, int xRadius, int yRadius, float orientation, double viewShift)
    {
        GL gl = glPanel.getGL();

        if (gl == null) return;

        try
        {
            if (xyz == null) return;

            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            color.setGLColor(gl);
            gl.glLineWidth(1.0f);
            gl.glEnable(GL.GL_POLYGON_STIPPLE);
            gl.glPolygonStipple(StsGraphicParameters.halftone, 1);

            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            gl.glRotatef(orientation, 0.0f, 0.0f, 1.0f);
            gl.glBegin(GL.GL_POLYGON);
            for (int i = 0; i < 360; i = i + 10)
            {
                float degInRad = i * DEG2RAD;
                gl.glVertex3f((float) (Math.cos(degInRad) * xRadius), (float) (Math.sin(degInRad) * yRadius), 0.0f);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glPopMatrix();
            gl.glDisable(GL.GL_POLYGON_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
            gl.glFlush();
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    static public void drawEllipse(float[] xyz, StsColor color, StsGLPanel glPanel, int xRadius, int yRadius, float orientation, double viewShift)
    {
        GL gl = glPanel.getGL();

        if (gl == null) return;

        try
        {
            if (xyz == null) return;

            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            color.setGLColor(gl);
            gl.glLineWidth(2.0f);

            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            gl.glRotatef(orientation, 0.0f, 0.0f, 1.0f);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int i = 0; i < 360; i = i + 10)
            {
                float degInRad = i * DEG2RAD;
                gl.glVertex3f((float) (Math.cos(degInRad) * xRadius), (float) (Math.sin(degInRad) * yRadius), 0.0f);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glPopMatrix();
            gl.glEnable(GL.GL_LIGHTING);
            gl.glFlush();
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    private static void drawLineStrip(GL gl, StsColor color, StsPoint[] points)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            gl.glLineWidth(1.0f);
            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < points.length; n++)
                gl.glVertex3fv(points[n].v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawLineStrip(GL gl, StsColor color, float[][] points)
    {
        drawLineStrip(gl, color, points, -1);
    }

    public static void drawLineStrip(GL gl, StsColor color, double[][] points)
    {
        drawLineStrip(gl, color, points, -1);
    }

    public static void drawLineStrip(GL gl, StsColor color, float[][] points, int width)
    {
        if (gl == null || color == null || points == null) return;

        if (width != -1)
        {
            gl.glLineWidth(width);
        }
        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < points.length; n++)
                gl.glVertex3fv(points[n], 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawLineStrip(GL gl, StsColor color, double[][] points, int width)
    {
        if (gl == null || color == null || points == null) return;

        if (width != -1)
        {
            gl.glLineWidth(width);
        }
        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < points.length; n++)
                gl.glVertex3dv(points[n], 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private static void drawLineStrip(GL gl, StsColor color, StsPoint[] points, int zIndex)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < points.length; n++)
            {
                float[] v = points[n].v;
                gl.glVertex3f(v[0], v[1], v[zIndex]);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private static void drawLineStrip(GL gl, StsColor color, StsPoint[] points, int min, int max, StsPoint topPoint,
                                      StsPoint botPoint, int zIndex)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            if (topPoint != null)
            {
                float[] v = topPoint.v;
                gl.glVertex3f(v[0], v[1], v[zIndex]);
            }
            for (int n = min; n <= max; n++)
            {
                float[] v = points[n].v;
                gl.glVertex3f(v[0], v[1], v[zIndex]);
            }
            if (botPoint != null)
            {
                float[] v = botPoint.v;
                gl.glVertex3f(v[0], v[1], v[zIndex]);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }
    private static void drawLineStrip(GL gl, StsColor color, StsPoint[] points, int min, int max, int zIndex)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);

            for (int n = min; n <= max; n++)
            {
                float[] v = points[n].v;
                gl.glVertex3f(v[0], v[1], v[zIndex]);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawLineStrip2d(GL gl, StsColor color, StsPoint[] points)
    {
        drawLineStrip2d(gl, color, points, 2);
    }

    /**
     * draw in 2d using 2 of 3 xyz coordinates
     *
     * @param gl
     * @param color  color of line
     * @param points list of points
     * @param nCoor  coordinate direction not used (orthogonal to 2d display)
     */
    static public void drawLineStrip2d(GL gl, StsColor color, StsPoint[] points, int nCoor, int width)
    {
        gl.glLineWidth(width);
        drawLineStrip2d(gl, color, points, nCoor);
    }

    static public void drawLineStrip2d(GL gl, StsColor color, StsPoint[] points, int nCoor)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (points.length < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);

            switch (nCoor)
            {
                case 0:
                    for (int n = 0; n < points.length; n++)
                        gl.glVertex2f(points[n].v[1], points[n].v[2]);
                    break;
                case 1:
                    for (int n = 0; n < points.length; n++)
                        gl.glVertex2f(points[n].v[0], points[n].v[2]);
                    break;
                case 2:
                    for (int n = 0; n < points.length; n++)
                        gl.glVertex2f(points[n].v[0], points[n].v[1]);
            }

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawLineStrip2d(GL gl, StsColor color, StsPoint[] points, int min, int max, int nXAxis2d, int nYAxis2d)
    {
       if (gl == null || color == null || points == null) return;

        try
        {
            if (max - min < 1) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);

            for (int n = min; n <= max; n++)
                gl.glVertex2f(points[n].v[nXAxis2d], points[n].v[nYAxis2d]);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawLineStrip2d(GL gl, StsColor color, float[] x, float[] adjustments, float multiplier, float zMin, float zInc)
    {
        if (gl == null || color == null || x == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            int nValues = x.length;
            if (nValues < 2) return;
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            float z = zMin;
            for (int n = 0; n < nValues; n++, z += zInc)
                gl.glVertex2f(x[n] + multiplier * adjustments[n], z);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawDottedLineStrip(GL gl, StsColor color, float[][] points, int width)
    {
        if (gl == null || color == null || points == null) return;
        if (points.length < 2) return;
        if (width != -1) gl.glLineWidth(width);
        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_LINE_STIPPLE);

            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            drawLineStrip(gl, color, points, width);

            color = StsColor.BLACK;
            gl.glLineStipple(1, StsGraphicParameters.altDottedLine);

            drawLineStrip(gl, color, points, width);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glDisable(GL.GL_LINE_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawDottedLineStrip2d(GL gl, StsColor color, StsPoint[] points, int lineWidth)
    {
        drawDottedLineStrip2d(gl, color, points, lineWidth, 2);
    }

    static public void drawDottedLineStrip2d(GL gl, StsColor color, float[] x, float[] adjustments, float multiplier, float zMin, float zInc, int lineWidth)
    {
        if (gl == null || color == null || x == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
            if (x.length < 2) return;

            gl.glLineWidth(lineWidth);

            gl.glEnable(GL.GL_LINE_STIPPLE);

            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            drawLineStrip2d(gl, color, x, adjustments, multiplier, zMin, zInc);

            color = StsColor.BLACK;
            gl.glLineStipple(1, StsGraphicParameters.altDottedLine);

            drawLineStrip2d(gl, color, x, adjustments, multiplier, zMin, zInc);

            gl.glDisable(GL.GL_LINE_STIPPLE);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
    }

    static public void drawDottedLineStrip2d(GL gl, StsColor color, StsPoint[] points, int lineWidth, int nCoor)
    {
        if (gl == null || color == null || points == null) return;

        try
        {
            if (points.length < 2) return;
            gl.glLineWidth(lineWidth);
            gl.glEnable(GL.GL_LINE_STIPPLE);

            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            drawLineStrip2d(gl, color, points, nCoor);

            gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
            drawLineStrip2d(gl, StsColor.BLACK, points, nCoor);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glDisable(GL.GL_LINE_STIPPLE);
        }
    }

    public static void pickLineStrip(GL gl, StsColor color, StsPoint[] points)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */
        //        gl.glDisable(GL.GL_LIGHTING);

        /** Note that PushName and PopName must be OUTSIDE the gl Begin/End loop
         * This is more expensive than drawing method previously which is why
         * we have separate methods (in some cases) for displaying and picking.
         */
        try
        {
            for (int n = 0; n < points.length - 1; n++)
            {
                gl.glPushName(n);

                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex3fv(points[n].v, 0);
                gl.glVertex3fv(points[n + 1].v, 0);
                gl.glEnd();

                gl.glPopName();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.pickLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            //            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void pickLineStrip(GL gl, StsColor color, StsPoint[] points, int zIndex)
    {
        if (gl == null || color == null || points == null) return;

        /** IMPORTANT: turn LIGHTING off and then back on for lines */
        //        gl.glDisable(GL.GL_LIGHTING);

        /** Note that PushName and PopName must be OUTSIDE the gl Begin/End loop
         * This is more expensive than drawing method previously which is why
         * we have separate methods (in some cases) for displaying and picking.
         */
        try
        {
            float[] v1 = points[0].v;
            if (v1.length > (zIndex))
            { // jbw
                for (int n = 0; n < points.length - 1; n++)
                {
                    gl.glPushName(n);

                    gl.glBegin(GL.GL_LINE_STRIP);
                    float[] v0 = v1;
                    v1 = points[n + 1].v;

                    if ((v0.length > (zIndex)) && (v1.length > (zIndex)))
                    {
                        gl.glVertex3f(v0[0], v0[1], v0[zIndex]);
                        gl.glVertex3f(v1[0], v1[1], v1[zIndex]);
                        gl.glEnd();
                    }
                    else
                    {
                        gl.glEnd();
                        System.out.println("Pickline error: vertex index " + n + " is only " + v0.length);
                        // jbw no T entry, for example

                    }
                    gl.glPopName();
                }
            }
            else
            {
                System.out.println("Pickline error: vertex index out of range " + zIndex);
            }

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.pickLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            //            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip(gl, color, points);
    }

    public static void drawLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int zIndex)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip(gl, color, points, zIndex);
    }

    public static void drawLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max,
                                StsPoint topPoint, StsPoint botPoint, int zIndex)
    {
        if (highlighted)
        {
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
            //StsGLDraw.drawHighlightedLine(gl, color, highlighted, points, min, max, topPoint, botPoint, zIndex);
        }
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);
    }

    public static void drawLine(GL gl, StsColor color, boolean highlighted, double[][] points)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip(gl, color, points);
    }

    public static void drawLine2d(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int nCoor)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip2d(gl, color, points, nCoor);
    }

    public static void drawLine2d(GL gl, StsColor color, float width, StsPoint[] points, int nCoor)
    {
        gl.glLineWidth(width);
        drawLineStrip2d(gl, color, points, nCoor);
    }

    public static void drawLine2d(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max, int nXAxis2d, int nYAxis2d)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawLineStrip2d(gl, color, points, min, max, nXAxis2d, nYAxis2d);
    }

    public static void pickLine(GL gl, StsColor color,
                                boolean highlighted, StsPoint[] points)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        pickLineStrip(gl, color, points);
    }

    public static void pickLine(GL gl, StsColor color,
                                boolean highlighted, StsPoint[] points, int zIndex)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        pickLineStrip(gl, color, points, zIndex);
    }

    public static void drawDottedLine2d(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max, int nXAxis2d, int nYAxis2d)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawDottedLineStrip2d(gl, color, StsColor.BLACK, points, min, max, nXAxis2d, nYAxis2d);
    }

    public static void drawHighlightedLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max, StsPoint topPoint, StsPoint botPoint, int zIndex)
    {
        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);

        color = StsColor.BLACK;
        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);

        gl.glDisable(GL.GL_LINE_STIPPLE);

        color = StsColor.GRAY;
        gl.glLineWidth(StsGraphicParameters.edgeLineWidthHighlighted);
        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);
    }

    public static void drawDottedLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max, StsPoint topPoint, StsPoint botPoint, int zIndex)
    {
        float edgeWidth;
        if (highlighted)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        gl.glLineWidth(edgeWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);

        color = StsColor.BLACK;
        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, topPoint, botPoint, zIndex);

        gl.glDisable(GL.GL_LINE_STIPPLE);
    }

    public static void drawDottedLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int min, int max, int zIndex)
    {
        float edgeWidth;
        if (highlighted)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        gl.glLineWidth(edgeWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, zIndex);

        color = StsColor.BLACK;
        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, min, max, zIndex);

        gl.glDisable(GL.GL_LINE_STIPPLE);
    }

    public static void drawDottedLine(GL gl, StsColor color, boolean highlighted, StsPoint[] points, int zIndex)
    {
        float edgeWidth;
        if (highlighted)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;
        drawDottedLine(gl, color, edgeWidth, points, zIndex);
    }

    public static void drawDottedLine(GL gl, StsColor color, float edgeWidth, StsPoint[] points, int zIndex)
    {
        gl.glLineWidth(edgeWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, zIndex);

        color = StsColor.BLACK;
        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points, zIndex);

        gl.glDisable(GL.GL_LINE_STIPPLE);
    }

    public static void drawDottedLine(GL gl, StsColor color, float lineWidth, StsPoint[] points)
    {
        gl.glLineWidth(lineWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points);

        color = StsColor.BLACK;
        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        //        drawLineStrip(gl, color, points);
    }

    public static void drawDottedLine(GL gl, StsColor color, float edgeWidth, StsPoint[] points, StsColor backgroundColor)
    {
        gl.glLineWidth(edgeWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, color, points);

        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLineStrip(gl, backgroundColor, points);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        //        drawLineStrip(gl, color, points);
    }

    static public void drawDottedLineStrip2d(GL gl, StsColor color, StsColor backgroundColor, StsPoint[] points, int min, int max, int nXAxis2d, int nYAxis2d)
    {
       if (gl == null || color == null || points == null) return;

        try
        {
            if (max - min < 1) return;
            gl.glDisable(GL.GL_LIGHTING);

            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
            drawLineStrip2d(gl, color, points, min, max, nXAxis2d, nYAxis2d);

            gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
            drawLineStrip2d(gl, backgroundColor, points, min, max, nXAxis2d, nYAxis2d);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        gl.glEnable(GL.GL_LIGHTING);

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawDottedLine2d(GL gl, StsColor color, StsColor backgroundColor, float edgeWidth,
                                        StsPoint point0, StsPoint point1, int[] coorIndexes)
    {
        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(edgeWidth);

        float[] xy0 = new float[]
            {
                point0.v[coorIndexes[0]], point0.v[coorIndexes[1]]};
        float[] xy1 = new float[]
            {
                point1.v[coorIndexes[0]], point1.v[coorIndexes[1]]};

        //        StsGLDraw.drawPoint2d(xy0, color, gl, 10);
        //        StsGLDraw.drawPoint2d(xy1, color, gl, 10);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, color, xy0, xy1);

        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, backgroundColor, xy0, xy1);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void drawDottedLine2d(GL gl, StsColor color, StsColor backgroundColor, float edgeWidth,
                                        StsPoint point0, StsPoint point1)
    {
        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(edgeWidth);

        float[] xy0 = new float[]{point0.v[0], point0.v[1]};
        float[] xy1 = new float[]{point1.v[0], point1.v[1]};

        //        StsGLDraw.drawPoint2d(xy0, color, gl, 10);
        //        StsGLDraw.drawPoint2d(xy1, color, gl, 10);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, color, xy0, xy1);

        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, backgroundColor, xy0, xy1);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void drawDottedLine2d(GL gl, StsColor color, StsColor backgroundColor, float edgeWidth,
                                        float x0, float y0, float x1, float y1)
    {
        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(edgeWidth);

        float[] xy0 = new float[]{x0, y0};
        float[] xy1 = new float[]{x1, y1};

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, color, xy0, xy1);

        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        drawLine2d(gl, backgroundColor, xy0, xy1);

        gl.glDisable(GL.GL_LINE_STIPPLE);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void drawLine2d(GL gl, StsColor color, float[] xy0, float[] xy1)
    {
        try
        {
            color.setGLColor(gl);
            
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex2fv(xy0, 0);
            gl.glVertex2fv(xy1, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLine2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
        }
    }

    public static void drawSolidEdge(GL gl, StsColor color, float edgeWidth,
                                     StsPoint[] points)
    {
        gl.glLineWidth(edgeWidth);
        drawLineStrip(gl, color, points);

        gl.glLineWidth(2.0f * edgeWidth);
        color = StsColor.BLACK;
        drawLineStrip(gl, color, points);
    }

    public static void pickEdge(GL gl, StsColor color, float edgeWidth,
                                StsPoint[] points)
    {
        gl.glLineWidth(edgeWidth);
        pickLineStrip(gl, color, points);
    }

    public static void pickEdge(GL gl, StsColor color, float edgeWidth,
                                StsPoint[] points, int zIndex)
    {
        gl.glLineWidth(edgeWidth);
        pickLineStrip(gl, color, points, zIndex);
    }

    public static void drawGridLines(GL gl, StsColor color, float lineWidth,
                                     StsPoint[][] points, int nRows, int nCols)
    {
        try
        {
            if (points == null || color == null) return;

            /** IMPORTANT: turn LIGHTING off and then back on for lines */
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glLineWidth(lineWidth);

            for (int i = 0; i < nRows; i++)
            {
                gl.glBegin(GL.GL_LINE_STRIP);
                for (int j = 0; j < nCols; j++)
                    gl.glVertex3fv(points[i][j].v, 0);
                gl.glEnd();
            }

            for (int j = 0; j < nCols; j++)
            {
                gl.glBegin(GL.GL_LINE_STRIP);
                for (int i = 0; i < nRows; i++)
                    gl.glVertex3fv(points[i][j].v, 0);
                gl.glEnd();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawGridLines() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    /** draw first and last rows and cols, but skip some of the inner rows and cols */
    public static void drawGridLines(GL gl, StsColor color, float lineWidth,
                                     StsPoint[][] points, int nRows, int nCols, int rowStride, int colStride)
    {
        int i, j;

        try
        {
            if (points == null) return;

            /** IMPORTANT: turn LIGHTING off and then back on for lines */
            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glLineWidth(lineWidth);

            gl.glBegin(GL.GL_LINE_STRIP);
            for (j = 0; j < nCols; j++)
                gl.glVertex3fv(points[0][j].v, 0);
            gl.glEnd();

            for (i = 1; i < nRows - 1; i += rowStride)
            {
                gl.glBegin(GL.GL_LINE_STRIP);
                for (j = 0; j < nCols; j++)
                    gl.glVertex3fv(points[i][j].v, 0);
                gl.glEnd();
            }

            gl.glBegin(GL.GL_LINE_STRIP);
            for (j = 0; j < nCols; j++)
                gl.glVertex3fv(points[nRows - 1][j].v, 0);
            gl.glEnd();

            gl.glBegin(GL.GL_LINE_STRIP);
            for (i = 0; i < nRows; i++)
                gl.glVertex3fv(points[i][0].v, 0);
            gl.glEnd();

            for (j = 1; j < nCols - 1; j += colStride)
            {
                gl.glBegin(GL.GL_LINE_STRIP);
                for (i = 0; i < nRows; i++)
                    gl.glVertex3fv(points[i][j].v, 0);
                gl.glEnd();
            }

            gl.glBegin(GL.GL_LINE_STRIP);
            for (i = 0; i < nRows; i++)
                gl.glVertex3fv(points[i][nCols - 1].v, 0);
        }
        catch (Exception e)
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    /**
     * Draw a square point "width" pixels wide and "height" pixels high centered
     * at point on the screen, shifted "shift" screen units away from the viewer.
     * In general, don't do viewShifting on a point by point basis, i.e, call with
     * shift == 0; for lots of points, call it once before starting, and then
     * resetViewShift afterwards;
     */
    static public void drawPoint(float[] xyz, StsColor color, StsGLPanel glPanel,
                                 int width, int height, double viewShift)
    {
        drawPoint(xyz, color, glPanel, width, height, viewShift, 0.0, 0.0);
    }

    static public void drawPoint(float[] xyz, StsColor color, StsGLPanel glPanel,
                                 int width, int height, double viewShift,
                                 double leftRightShift, double upDownShift)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {
            if (xyz == null) return;
            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            double halfHeight = height / 2;
            double[] screenPoint = glPanel.getScreenCoordinates(xyz);
            screenPoint[0] += leftRightShift;
            screenPoint[1] += halfHeight + upDownShift;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint topPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            screenPoint[1] -= height;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint botPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            gl.glDisable(GL.GL_LIGHTING);

            if (color != null) color.setGLColor(gl);
            gl.glLineWidth((float) width);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex3fv(topPoint.v, 0);
            gl.glVertex3fv(botPoint.v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    static public void drawBitMap(float[] xyz, StsColor color, StsGLPanel glPanel)
    {
        StsBitmap.fontBitmapF.glBitmap(glPanel.getGL(), xyz, color);
    }

    static public void drawTriangle(float[] xyz, StsColor color, StsGLPanel glPanel, byte size)
    {
        StsBitmap.drawTriangle(glPanel.getGL(), xyz, color, size);
    }

    static public void drawTriangle2d(float[] xyz, StsColor color, GL gl, byte size)
    {
        StsBitmap.drawTriangle2d(gl, xyz, color, size);
    }

    static public void drawDiamond2d(float[] xyz, StsColor color, GL gl, byte size)
    {
        StsBitmap.drawDiamond2d(gl, xyz, color, size);
    }

    static public void drawDiamond(float[] xyz, StsColor color, StsGLPanel glPanel, byte size)
    {
        StsBitmap.drawDiamond(glPanel.getGL(), xyz, color, size);
    }

    static public void drawCircle2d(float[] xyz, StsColor color, GL gl, byte size)
    {
        StsBitmap.drawCircle2d(gl, xyz, color, size);
    }

    static public void drawCircle(float[] xyz, StsColor color, StsGLPanel glPanel, byte size)
    {
        StsBitmap.drawCircle(glPanel.getGL(), xyz, color, size);
    }

    static public void drawStar(float[] xyz, StsColor color, StsGLPanel glPanel, byte size)
    {
        StsBitmap.drawStar(glPanel.getGL(), xyz, color, size);
    }

    static public void drawStar2d(float[] xyz, StsColor color, GL gl, byte size)
    {
        StsBitmap.drawStar2d(gl, xyz, color, size);
    }

    static public void drawEquipment2d(float[] xyz, GL gl, StsColor color, byte type)
    {
        StsBitmap.drawEquipmentBitmap2d(gl, xyz, color, type);
    }

    static public void drawEquipment2d(int x, int y, GL gl, StsColor color, byte type)
    {
        StsBitmap.drawEquipmentBitmap2d(gl, x, y, color, type);
    }

    static public void drawEquipment(float[] xyz, GL gl, StsColor color, byte type)
    {
        StsBitmap.drawEquipmentBitmap(gl, xyz, color, type);
    }

    static public void drawWellhead2d(float[] xyz, GL gl, StsColor color, byte type)
    {
        StsBitmap.drawWellBitmap2d(gl, xyz, color, type);
    }

    static public void drawWellhead(float[] xyz, GL gl, StsColor color, byte type)
    {
        StsBitmap.drawWellBitmap(gl, xyz, color, type);
    }

    static public void drawEquipment2d(float[] xyz, StsGLPanel glPanel, StsColor color, byte type)
    {
        StsBitmap.drawEquipmentBitmap2d(glPanel.getGL(), xyz, color, type);
    }

    static public void drawEquipment(float[] xyz, StsGLPanel glPanel, StsColor color, byte type)
    {
        StsBitmap.drawEquipmentBitmap(glPanel.getGL(), xyz, color, type);
    }

    static public void drawWellhead2d(float[] xyz, StsGLPanel glPanel, StsColor color, byte type)
    {
        StsBitmap.drawWellBitmap2d(glPanel.getGL(), xyz, color, type);
    }

    static public void drawWellhead(float[] xyz, StsGLPanel glPanel, StsColor color, byte type)
    {
        StsBitmap.drawWellBitmap(glPanel.getGL(), xyz, color, type);
    }

    /**
     * Draw a line in 3d starting at this xyz with width(pixels) and length(pixels) in
     * 1 of 4 cardinal directions.  The line can be view shifted and/or offset in pixels
     * by leftRightShift and updownShift from the original starting point.
     */
    static public void drawLine(float[] xyz, StsColor color, StsGLPanel glPanel,
                                int width, int length, int direction, double viewShift,
                                double leftRightShift, double upDownShift)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {
            if (xyz == null) return;
            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            double[] screenPoint = glPanel.getScreenCoordinates(xyz);
            screenPoint[0] += leftRightShift;
            screenPoint[1] += upDownShift;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint topPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            if (direction == StsParameters.NORTH)
                screenPoint[1] += length;
            else if (direction == StsParameters.SOUTH)
                screenPoint[1] -= length;
            else if (direction == StsParameters.EAST)
                screenPoint[0] += length;
            else if (direction == StsParameters.WEST)
                screenPoint[0] -= length;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint botPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glLineWidth((float) width);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex3fv(topPoint.v, 0);
            gl.glVertex3fv(botPoint.v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    /**
     * Draw a line in 3d starting at this xyz with width(pixels) and length(pixels) in
     * 1 of 4 cardinal directions.  The line can be view shifted and/or offset in pixels
     * by leftRightShift and updownShift from the original starting point.
     */
    static public void drawLine2d(float x, float y, StsColor color, StsGLPanel glPanel,
                                  int width, int length, int direction,
                                  double leftRightShift, double upDownShift)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {

            double[] screenPoint = glPanel.getScreenCoordinates(x, y);
            screenPoint[0] += leftRightShift;
            screenPoint[1] += upDownShift;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint topPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            if (direction == StsParameters.NORTH)
                screenPoint[1] += length;
            else if (direction == StsParameters.SOUTH)
                screenPoint[1] -= length;
            else if (direction == StsParameters.EAST)
                screenPoint[0] += length;
            else if (direction == StsParameters.WEST)
                screenPoint[0] -= length;
            //TODO instead, use getWorldCoordinates which returns double[]
            StsPoint botPoint = glPanel.getWorldCoordinatesPoint(screenPoint);

            gl.glDisable(GL.GL_LIGHTING);

            color.setGLColor(gl);
            gl.glLineWidth((float) width);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex3fv(topPoint.v, 0);
            gl.glVertex3fv(botPoint.v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawCylinder(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size, float height)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        GLU glu = glPanel3d.getGLU();
        if ((gl == null) || (glu == null)) return;

        try
        {
            gl.glEnable(gl.GL_LIGHTING);
            gl.glEnable(GL.GL_BLEND);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);
            // gl.glShadeModel(gl.GL_SMOOTH);
            color.setGLColor(gl);
            GLUquadric qobj = glu.gluNewQuadric();

            gl.glPushMatrix();

            glu.gluQuadricDrawStyle(qobj, GLU.GLU_FILL);
            glu.gluQuadricNormals(qobj, GLU.GLU_SMOOTH);
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            float zscale = glPanel3d.getZScale();
            gl.glScalef(1.0f, 1.0f, 1.0f / zscale);

            glu.gluDisk(qobj, 0.0, size, 40, 1);
            gl.glTranslatef(0.0f, 0.0f, -height);
            glu.gluCylinder(qobj, size, size, height, 40, 1);
            glu.gluDisk(qobj, 0.0, size, 40, 1);

            gl.glPopMatrix();

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            gl.glFlush();
        }
    }

    static public void drawRotatedCylinder(StsGLPanel3d glPanel3d, float[][] xyz, StsColor color, float size, float height)
    {
        float dx = xyz[0][0] - xyz[1][0];
        float dy = xyz[0][1] - xyz[1][1];
        float dz = xyz[0][2] - xyz[1][2];

        // Compute rotation
        double xydist = Math.sqrt(dx * dx + dy * dy);
        int xRotation = 0;
        int yRotation = (int) (StsMath.atan2(dx, dy));
        int zRotation = (int) (Math.tan(xydist / dz) - 90);
        //System.out.println("Rotation X=" + xRotation + " Y=" + yRotation + " Z=" + zRotation);
        drawRotatedCylinder(glPanel3d, xyz[0], color, size, height, xRotation, yRotation, zRotation);
    }

    static public void drawRotatedCylinder(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size, float height,
                                           int xRotation, int yRotation, int zRotation)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        GLU glu = glPanel3d.getGLU();
        if ((gl == null) || (glu == null)) return;

        try
        {
            gl.glEnable(gl.GL_LIGHTING);
            gl.glEnable(GL.GL_BLEND);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);
            color.setGLColor(gl);
            GLUquadric qobj = glu.gluNewQuadric();

            gl.glPushMatrix();

            glu.gluQuadricDrawStyle(qobj, GLU.GLU_FILL);
            glu.gluQuadricNormals(qobj, GLU.GLU_SMOOTH);
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            float zscale = glPanel3d.getZScale();
            gl.glScalef(1.0f, 1.0f, 1.0f / zscale);

            gl.glRotatef(xRotation, 1, 0, 0);
            gl.glRotatef(yRotation, 0, 1, 0);
            gl.glRotatef(zRotation, 0, 0, 1);

            glu.gluDisk(qobj, 0.0, size, 10, 1);
            gl.glTranslatef(0.0f, 0.0f, -height);
            glu.gluCylinder(qobj, size, size, height, 10, 1);
            glu.gluDisk(qobj, 0.0, size, 10, 1);

            gl.glRotatef(0, 1, 0, 0);
            gl.glRotatef(0, 0, 1, 0);
            gl.glRotatef(0, 0, 0, 1);

            gl.glPopMatrix();
            //gl.glPopMatrix();

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            gl.glFlush();
        }
    }
    static public void drawSphere(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        try
        {
            // gl.glEnable(GL.GL_LIGHTING);
            // gl.glShadeModel(GL.GL_SMOOTH);
            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            float zscale = glPanel3d.getZScale();
            float scaleFactor = size/sphereSize;
            gl.glScalef(scaleFactor, scaleFactor, scaleFactor / zscale);
            color.setGLColor(gl);
            if (displayListSphere == 0 && !constructSphereDisplayList(glPanel3d, sphereSize)) return;
            gl.glCallList(displayListSphere);
        }
        finally
        {
            gl.glFlush();
            gl.glPopMatrix();
        }
    }

    static public void drawEllipsoid(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size, float azimuth)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        try
        {
            gl.glEnable(GL.GL_LIGHTING);
            gl.glShadeModel(GL.GL_SMOOTH);
            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            gl.glRotatef(azimuth, 1.0f, 0.0f, 1.0f);
            float zscale = glPanel3d.getZScale();
            float scaleFactor = size/sphereSize;
            gl.glScalef(0.8f*scaleFactor, 0.1f*scaleFactor, 0.1f*scaleFactor / zscale);
            color.setGLColor(gl);
            if (displayListSphere == 0 && !constructSphereDisplayList(glPanel3d, sphereSize)) return;
            gl.glCallList(displayListSphere);
        }
        finally
        {
            gl.glFlush();
            gl.glPopMatrix();
        }
    }
    static public void drawDisk2d(StsGLPanel3d glPanel3d, float[] xy, StsColor color, float size)
    {
        GL gl = glPanel3d.getGL();
        gl.glDisable(GL.GL_LIGHTING);
        float[] xyz = new float[] { xy[0], xy[1], 0.0f };
        drawDisk3d(glPanel3d, xyz, StsColor.BLACK, size+2, 0.0f, 0.0f, false);
        drawDisk3d(glPanel3d, xyz, color, size, 0.0f, 0.0f, true);
        gl.glEnable(GL.GL_LIGHTING);
    }

    static public void drawDisk(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size, float azimuth, float elevation)
    {
        drawDisk3d(glPanel3d, xyz, color, size, 0.0f, 0.0f, true);
    }
    static public void drawDisk3d(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size, float azimuth, float elevation, boolean filled)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        GLU glu = glPanel3d.getGLU();
        if ((gl == null) || (glu == null)) return;
        try
        {
            color.setGLColor(gl);
            gl.glShadeModel(gl.GL_SMOOTH);

            GLUquadric qobj = glu.gluNewQuadric();

            gl.glPushMatrix();
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            gl.glRotatef(elevation, 1.0f, 0.0f, 0.0f);
            gl.glRotatef(azimuth, 1.0f, 0.0f, 1.0f);
            gl.glScalef(1.0f, 1.0f, 1.0f);

            gl.glLineWidth(1.0f);
            if(filled)
                glu.gluQuadricDrawStyle(qobj, GLU.GLU_FILL);
            else
                glu.gluQuadricDrawStyle(qobj, GLU.GLU_LINE);
            glu.gluQuadricNormals(qobj, GLU.GLU_SMOOTH);
            glu.gluDisk(qobj, 0, size, 32, 8);
            gl.glPopMatrix();
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawSphere() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glFlush();
        }
    }

    static public void displayTestBlobs(StsGLPanel3d glPanel3d)
    {
        float[] xyz;
        xyz = new float[]{10.0f, 10.0f, 0.0f};
        StsGLDraw.drawSphere(glPanel3d, xyz, StsColor.GREEN, 10);
        xyz = new float[]{25.0f, 25.0f, 50.0f};
        StsGLDraw.drawCylinder(glPanel3d, xyz, StsColor.PURPLE, 10, 10);
        xyz = new float[]{50.0f, 50.0f, 100.0f};
        StsGLDraw.drawCube(glPanel3d, xyz, StsColor.CYAN, 10);
    }

    static public void drawCube(StsGLPanel3d glPanel3d, float[] xyz, StsColor color, float size)
    {
        if (xyz == null) return;
        GL gl = glPanel3d.getGL();
        GLU glu = glPanel3d.getGLU();
        if ((gl == null) || (glu == null)) return;
        try
        {
            color.setGLColor(gl);
            gl.glEnable(gl.GL_LIGHTING);
            gl.glShadeModel(gl.GL_SMOOTH);

            GLUquadric qobj = glu.gluNewQuadric();
            gl.glPushMatrix();
            glu.gluQuadricDrawStyle(qobj, GLU.GLU_FILL);
            glu.gluQuadricNormals(qobj, GLU.GLU_SMOOTH);
            gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
            float zscale = glPanel3d.getZScale();
            gl.glScalef(1.0f, 1.0f, 1.0f / zscale);
            glu.gluDisk(qobj, 0.0, size, 4, 1);
            gl.glTranslatef(0.0f, 0.0f, -size);
            glu.gluCylinder(qobj, size, size, size, 4, 4);
            glu.gluDisk(qobj, 0.0, size, 4, 1);
            gl.glPopMatrix();
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawSphere() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glFlush();
        }
    }

    static public void drawPoint(GL gl, float[] xyz, StsColor color, int size)
    {
        try
        {
            if (xyz == null) return;

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3fv(xyz, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawTriangle(float[] xyz, StsColor color, StsGLPanel glPanel, int size)
    {
        drawTriangle(xyz, color, glPanel, size, 0);
    }

    static public void drawTriangle(float[] xyz, StsColor color, StsGLPanel glPanel, int size, double viewShift)
    {
        GL gl = glPanel.getGL();

        if (gl == null) return;

        try
        {
            if (xyz == null) return;

            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3fv(xyz, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    static public int getDisplayListSphere(StsGLPanel3d glPanel3d)
    {
        if(displayListSphere == 0) constructSphereDisplayList(glPanel3d, sphereSize);
        return displayListSphere;
    }

    static private boolean constructSphereDisplayList(StsGLPanel3d glPanel3d, float size)
    {
        if (displayListSphere > 0) return true;
        GL gl = glPanel3d.getGL();
        if (gl == null) return false;
        GLU glu = glPanel3d.getGLU();
        if (glu == null) return false;
        displayListSphere = gl.glGenLists(1);
        if (displayListSphere == 0)
        {
            StsMessageFiles.logMessage("System Error in StsGrid.displaySurface: " + "Failed to allocate a display list");
            return false;
        }

        gl.glNewList(displayListSphere, GL.GL_COMPILE);
        if(!drawQuadricSphere(glu, size)) return false;
        gl.glEndList();
        return true;
    }

    static boolean drawQuadricSphere(GLU glu, float size)
    {
         try
        {
            GLUquadric qobj = glu.gluNewQuadric();
            glu.gluQuadricDrawStyle(qobj, GLU.GLU_FILL);
            glu.gluQuadricNormals(qobj, GLU.GLU_SMOOTH);
            glu.gluSphere(qobj, size, 8, 8);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsGLDraw.class, "drawSphere", e);
            return false;
        }
    }

    static public void drawPoint(float[] xyz, StsColor color, StsGLPanel glPanel, int size)
    {
        drawPoint(xyz, color, glPanel, size, 0);
    }

    static public void drawPoint(float[] xyz, StsColor color, StsGLPanel glPanel, int size, double viewShift)
    {
        GL gl = glPanel.getGL();

        if (gl == null) return;

        try
        {
            if (xyz == null) return;

            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3fv(xyz, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    public static void drawPoint(float[] xyz, GL gl, int size)
    {

        drawPoint(gl, xyz, StsColor.YELLOW, size);
    }

    static public void drawPoint(double[] xyz, StsColor color, GL gl, int size)
    {
        try
        {
            if (xyz == null) return;

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3dv(xyz, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint(float[] v, int zIndex, StsColor color, StsGLPanel glPanel, int size, double viewShift)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {
            if (v == null) return;

            color.setGLColor(gl);
            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3f(v[0], v[1], v[zIndex]);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint(float x, float y, float z, StsColor color, GL gl, int size)
    {
        if (gl == null) return;

        try
        {
            color.setGLColor(gl);
            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3f(x, y, z);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint(float[] v, int zIndex, StsColor color, GL gl, int size)
    {
        if (gl == null) return;

        try
        {
            if (v == null) return;

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex3f(v[0], v[1], v[zIndex]);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }


    static public void drawPoint2d(StsPoint point, GL gl, int size)
    {
        drawPoint2d(point.v, null, gl, size);
    }

    static public void drawPoint2d(StsPoint point, StsColor color, GL gl, int size)
    {
        drawPoint2d(point.v, color, gl, size);
    }

    static public void drawPoint2d(float[] xy, GL gl, int size)
    {
        drawPoint2d(xy, null, gl, size);
    }

    static public void drawRectangle2d(float[] xy, StsColor color, GL gl, int width, int height)
    {
        try
        {
            if (xy == null) return;

            if (color != null) color.setGLColor(gl);

            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POLYGON);
            gl.glVertex2f(xy[0], xy[1]);
            gl.glVertex2f(xy[0] + width, xy[1]);
            gl.glVertex2f(xy[0] + width, xy[1] + height);
            gl.glVertex2f(xy[0], xy[1] + height);
            gl.glVertex2f(xy[0], xy[1]);

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawRectangle2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint2d(float[] xy, StsColor color, GL gl, int size)
    {
        try
        {
            if (xy == null) return;

            if (color != null) color.setGLColor(gl);
            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex2fv(xy, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint2d(float x, float y, GL gl, int size)
    {
        try
        {
            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex2f(x, y);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint2d(double[] xy, StsColor color, GL gl, int size)
    {
        try
        {
            if (xy == null) return;

            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex2dv(xy, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoint2d(float x, float y, StsColor color, StsGLPanel glPanel, int size, double viewShift)
    {
        GL gl = glPanel.getGL();
        if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);
        drawPoint2d(x, y, color, gl, size);
        if (viewShift != 0.0) glPanel.resetViewShift(gl);
    }

    static public void drawPoint2d(float x, float y, StsColor color, GL gl, int size)
    {
        if (gl == null) return;

        try
        {
            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            gl.glVertex2f(x, y);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoints2d(double[][] xy, StsColor color, GL gl, int size)
    {
        try
        {
            if (xy == null) return;
            int nPoints = xy.length;
            if (nPoints == 0) return;
            color.setGLColor(gl);
            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            for (int n = 0; n < nPoints; n++)
                gl.glVertex2dv(xy[n], 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoint2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    static public void drawPoints(StsPoint[] points, StsColor color, StsGLPanel glPanel, int size, double viewShift)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {
            if (points == null) return;
            int nPnts = points.length;
            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);
            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            for (int n = 0; n < nPnts; n++)
                gl.glVertex3fv(points[n].v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoints() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    static public void drawPoints(StsPoint[] points, StsColor color, StsGLPanel glPanel,
                                  int size, double viewShift, int zIndex)
    {
        GL gl = glPanel.getGL();
        if (gl == null) return;

        try
        {
            if (points == null) return;
            int nPnts = points.length;
            if (viewShift != 0.0) glPanel.setViewShift(gl, viewShift);
            color.setGLColor(gl);

            gl.glPointSize((float) size);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_POINTS);
            for (int n = 0; n < nPnts; n++)
                gl.glVertex3f(points[n].v[0], points[n].v[1], points[n].v[zIndex]);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawPoints() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            if (viewShift != 0.0) glPanel.resetViewShift(gl);
        }
    }

    /** draw a well line segment */
    public static void drawWellLineSegment(GL gl, StsColor color,
                                           boolean highlighted, StsPoint[] points, int start, int end,
                                           StsPoint topPoint, StsPoint basePoint)
    {
        if (highlighted)
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        else
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);

        drawSegmentedLineWithEnds(gl, color, points, start, end, topPoint,
            basePoint);
    }

    /** draw a line segment */
    public static void drawSegmentedLineWithEnds(GL gl, StsColor color,
                                                 StsPoint[] points, int start, int end,
                                                 StsPoint firstPoint, StsPoint lastPoint)
    {
        try
        {
            if (gl == null || color == null || points == null) return;
            if (start < 0) start = 0;
            if (end > points.length - 1) end = points.length - 1;

            gl.glLineWidth(1.0f);
            /** IMPORTANT: turn LIGHTING off and then back on for lines */
            gl.glDisable(GL.GL_LIGHTING);
            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);

            if (firstPoint != null) gl.glVertex3fv(firstPoint.v, 0);
            for (int n = start; n <= end; n++) gl.glVertex3fv(points[n].v, 0);
            if (lastPoint != null) gl.glVertex3fv(lastPoint.v, 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawSegmentedLineWithEnds() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawSegmentedLine(GL gl, float[][] xyzPoints, int zIndex)
    {
        try
        {
            gl.glLineWidth(1.0f);
            if (zIndex == 2) drawSegmentedLine(gl, xyzPoints);

            if (xyzPoints == null) return;
            int nPoints = xyzPoints.length;

            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_LINE_STRIP);

            for (int n = 0; n < nPoints; n++)
                if (xyzPoints[n] != null) gl.glVertex3f(xyzPoints[n][0], xyzPoints[n][1], xyzPoints[n][zIndex]);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawSegmentedLine() failed. ",
                e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawSegmentedLine(GL gl, float[][] xyzPoints)
    {
        try
        {
            gl.glLineWidth(1.0f);
            if (xyzPoints == null) return;
            int nPoints = xyzPoints.length;

            gl.glDisable(GL.GL_LIGHTING);
            gl.glBegin(GL.GL_LINE_STRIP);

            for (int n = 0; n < nPoints; n++)
                if (xyzPoints[n] != null) gl.glVertex3fv(xyzPoints[n], 0);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawSegmentedLine() failed. ",
                e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void drawLineSegment(GL gl, StsColor color, StsPoint firstPoint, StsPoint lastPoint)
    {
        drawLineSegment(gl, color, firstPoint.v,  lastPoint.v, 1.0f);
    }

    public static void drawLineSegment(GL gl, StsColor color, StsPoint firstPoint, StsPoint lastPoint, float lineWidth)
    {
        drawLineSegment(gl, color, firstPoint.v,  lastPoint.v, lineWidth);
    }

    public static void drawLineSegment(GL gl, StsColor color, float[] p1, float[] p2, float lineWidth)
    {
        try
        {
            gl.glLineWidth(lineWidth);
            gl.glDisable(GL.GL_LIGHTING);
            color.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex3fv(p1, 0);
            gl.glVertex3fv(p2, 0);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsGLDraw.class, "drawLineSegment", e);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
            gl.glLineWidth(1.0f);
        }
    }

    public static void drawQuadStrips(GL gl, StsGridBoundingBox boundingBox, StsList quadStrips,
                                      float[][] zValues, float[][][] normals, boolean displayGaps)
    {
        QuadStrip q = null;
        float[] point = new float[3];
        float[] normal;
        int i = -1, j = -1; // rowcol coordinates relative to boundingBox
        int ii = -1, jj = -1; // global rowcol coordinates relative to boundingBox

        try
        {
            float xMin = boundingBox.xMin;
            float xInc = boundingBox.xInc;
            float yMin = boundingBox.yMin;
            float yInc = boundingBox.yInc;
            int rowMin = boundingBox.rowMin;
            int colMin = boundingBox.colMin;

            // shouldn't have to do this:
            // someone is turning off the lights and not turning them back on!!
            gl.glEnable(GL.GL_LIGHTING);

            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                // if (!displayGaps && q.cellType == StsParameters.CELL_GAP) continue;

                i = q.rowNumber;
                ii = i - rowMin;
                point[1] = yMin + i * yInc;

                int jMin = q.firstCol;
                int jMax = q.lastCol;

                point[0] = xMin + jMin * xInc;

                gl.glBegin(GL.GL_QUAD_STRIP);
                for (j = jMin, jj = jMin - colMin; j <= jMax; j++, jj++)
                {
                    normal = normals[ii][jj];
                    if (normal != null) gl.glNormal3fv(normal, 0);
                    point[2] = zValues[ii][jj];
                    if (point[2] != nullValue) gl.glVertex3fv(point, 0);

                    normal = normals[ii + 1][jj];
                    if (normal != null) gl.glNormal3fv(normal, 0);
                    point[1] += yInc;
                    point[2] = zValues[ii + 1][jj];
                    if (point[2] != nullValue) gl.glVertex3fv(point, 0);

                    point[1] -= yInc;
                    point[0] += xInc;
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawQuadStrips() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
        }
    }

    public static void drawQuadStripProps(GL gl, StsGridBoundingBox boundingBox, StsList quadStrips,
                                          float[][] zValues, float[][][] normals, boolean displayGaps,
                                          int nLayer, StsColor[][] colors)
    {
        QuadStrip q = null;

        int i = -1, j = -1; // rowcol coordinates relative to boundingBox
        int ii = -1, jj = -1; // global rowcol coordinates relative to boundingBox

        float[] n0, n1;
        float[] p0 = new float[3];
        float[] p1 = new float[3];

        try
        {
            float xMin = boundingBox.xMin;
            float xInc = boundingBox.xInc;
            float yMin = boundingBox.yMin;
            float yInc = boundingBox.yInc;
            int rowMin = boundingBox.rowMin;
            int colMin = boundingBox.colMin;

            // shouldn't have to do this:
            // someone is turning off the lights and not turning them back on!!
            gl.glEnable(GL.GL_LIGHTING);

            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                // if (!displayGaps && q.cellType == StsParameters.CELL_GAP) continue;

                i = q.rowNumber;
                ii = i - rowMin;

                j = q.firstCol;
                jj = j - colMin;

                float x = xMin + xInc * j;
                float y = yMin + yInc * i;

                p0[0] = x;
                p0[1] = y;
                p0[2] = zValues[ii][jj];

                p1[0] = x;
                p1[1] = y + yInc;
                p1[2] = zValues[ii + 1][jj];

                n0 = normals[ii][jj];
                n1 = normals[ii + 1][jj];

                gl.glBegin(GL.GL_QUADS);
                for (; j < q.lastCol;)
                {
                    colors[ii][jj].setGLColor(gl);

                    if (n0 != null) gl.glNormal3fv(n0, 0);
                    if (p0[2] != nullValue) gl.glVertex3fv(p0, 0);
                    if (n1 != null) gl.glNormal3fv(n1, 0);
                    if (p1[2] != nullValue) gl.glVertex3fv(p1, 0);

                    j++;
                    jj++;

                    p0[0] += xInc;
                    p1[0] += xInc;
                    p0[2] = zValues[ii][jj];
                    p1[2] = zValues[ii + 1][jj];
                    /*
                                       if(p0[2] > 10000.0f)
                                       {
                                           System.out.println("bad point - i: " + i + " j: " + j + " z: " + p0[2]);
                                       }
                    */
                    n0 = normals[ii][jj];
                    n1 = normals[ii + 1][jj];

                    if (n1 != null) gl.glNormal3fv(n1, 0);
                    if (p1[2] != nullValue) gl.glVertex3fv(p1, 0);
                    if (n0 != null) gl.glNormal3fv(n0, 0);
                    if (p0[2] != nullValue) gl.glVertex3fv(p0, 0);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawQuadStripProps() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
        }
    }

    public static void drawPolygonZContours(StsPoint[] points, float contourInterval, StsColor color)
    {
        ArrayList<StsPoint> contourPoints = new ArrayList<StsPoint>();
        int nPoints = points.length;
        StsPoint point = points[nPoints-1];
        for(int n = 0; n < nPoints; n++)
        {
            StsPoint prevPoint = point;
            point = points[n];
            StsPoint[] segmentContourPoints = getInterpolatedZPoints(prevPoint, point, 0.0f, contourInterval);
        }
    }

    private static StsPoint[] getInterpolatedZPoints(StsPoint point0, StsPoint point1, float zMin, float zInc)
    {
        int fi0 = 0, fi1 = 0, finc = 1;
        int nPoints = 0;
        float f0 = (point0.getZ() - zMin)/zInc;
        float f1 = (point1.getZ() - zMin)/zInc;
        f0 -= (int)f0;
        f1 -= (int)f0;
        if(f0 == f1)
        {
            if(f0 == 0.0f)
                return new StsPoint[] { point0, point1 };
        }
        else if(f0 < f1)
        {
            fi0 = StsMath.ceiling(f0);
            fi1 = StsMath.floor(f1);
            finc = 1;
            nPoints = fi1 - fi0 + 1;
        }
        else
        {
            fi0 = StsMath.floor(f0);
            fi1 = StsMath.ceiling(f1);
            finc = 1;
            nPoints = fi0 - fi1 + 1;
        }
        StsPoint[] points = new StsPoint[nPoints];
        int nPoint = 0;
        for(float fi = fi0; fi != fi1; fi += finc)
        {
            points[nPoint++] = StsPoint.staticInterpolatePoints(point0, point1, fi);
        }
        return points;
    }

    public static void fontTimesRoman24(GL gl, float[] xyz, String string)
    {

        fontOutput(gl, xyz, string, GLTimesRoman24BitmapFont.getInstance(gl));
    }

    public static void fontGL9x15(GL gl, float[] xyz, String string)
    {
        fontOutput(gl, xyz, string, GL9x15BitmapFont.getInstance(gl));
    }

    public static void fontHelvetica12(GL gl, float[] xyz, String string)
    {
        fontOutput(gl, xyz, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void fontHelvetica12WithBackground(GL gl, float[] xyz, String string)
    {
        fontOutput(gl, xyz, string, GLHelvetica12BitmapFont.getInstance(gl), true);
    }

    public static void fontHelvetica12(GL gl, float x, float y, float z, String string)
    {
        fontOutput(gl, x, y, z, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void fontHelvetica12(GL gl, float x, float y, String string)
    {
        fontOutput(gl, x, y, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void highlightedFontHelvetica12(GL gl, float x, float y, String string)
    {
        fontOutput(gl, x, y, string, GLHelvetica12BitmapFont.getInstance(gl));
        StsColor.WHITE.setGLColor(gl);
        fontOutput(gl, x + 1, y + 1, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void highlightedFontHelvetica12(GL gl, int x, int y, String string)
    {
        fontOutput(gl, x, y, string, GLHelvetica12BitmapFont.getInstance(gl));
        StsColor.WHITE.setGLColor(gl);
        fontOutput(gl, x + 1, y + 1, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void fontHelvetica12(GL gl, int x, int y, String string)
    {
        fontOutput(gl, x, y, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static void fontHelvetica12(GL gl, double x, double y, String string)
    {
        fontOutput(gl, x, y, string, GLHelvetica12BitmapFont.getInstance(gl));
    }

    public static double getFontHelvetica12Height(GL gl)
    {
        return GLHelvetica12BitmapFont.getInstance(gl).getCharacterHeight('A');
    }

    public static void fontHelvetica18(GL gl, float[] xyz, String string)
    {
        fontOutput(gl, xyz, string, GLHelvetica18BitmapFont.getInstance(gl));
    }

    /** Draws a bitmapped string */
    public static void fontOutput(GL gl, float[] xyz, String str, GLBitmapFont font)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glRasterPos3f(xyz[0], xyz[1], xyz[2]);
        font.drawString(str);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void fontOutput(GL gl, double[] xyz, String str, GLBitmapFont font)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glRasterPos3d(xyz[0], xyz[1], xyz[2]);
        font.drawString(str);
        gl.glEnable(GL.GL_LIGHTING);
    }
    
    /** Draws a bitmapped string */
    public static void fontOutput(GL gl, float[] xyz, String str, GLBitmapFont font, boolean drawBackground)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glRasterPos3f(xyz[0], xyz[1], xyz[2]);
        if (drawBackground)
        {
            int fontHeight = font.getCharacterHeight('A');
            int fontWidth = font.getCharacterWidth('Z');
            int width = (int) (fontWidth * str.length() * 1.5);
            int height = (int) (fontHeight * 1.5);
            ByteBuffer image = BufferUtil.newByteBuffer(height * width * 3);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    byte c = StsMath.unsignedIntToUnsignedByte(0);
                    image.put(c);
                    image.put(c);
                    image.put(c);
                }
            }
            image.rewind();
            gl.glDrawPixels(width, height, gl.GL_RGB, gl.GL_UNSIGNED_BYTE, image);
        }
        font.drawString(str);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void fontOutput(GL gl, float x, float y, float z, String str, GLBitmapFont font)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glRasterPos3f(x, y, z);
        font.drawString(str);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public static void fontOutput(GL gl, float x, float y, String str, GLBitmapFont font)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glRasterPos2f(x, y);
        font.drawString(str);
        gl.glEnable(GL.GL_LIGHTING);
    }

    /*
        public static GLBitmapFont getFontTimesRoman24() { return new GLTimesRoman24BitmapFont(); }
        public static GLBitmapFont getFontGL9x15() { return new GL9x15BitmapFont(); }
        public static GLBitmapFont getFontHelvetica12() { return new GLHelvetica12BitmapFont(); }
        public static GLBitmapFont getFontHelvetica18() { return new GLHelvetica18BitmapFont(); }
    */
    public static void fontOutput(GL gl, int x, int y, String str, GLBitmapFont font)
    {
        gl.glRasterPos2i(x, y);
        font.drawString(str);
    }

    public static void verticalFontOutput(GL gl, int x, int y, String str, StsVerticalFont font)
    {
        gl.glRasterPos2i(x, y);
        font.drawString(gl, str);
    }

    public static void fontOutput(GL gl, double x, double y, String str, GLBitmapFont font)
    {
        gl.glRasterPos2d(x, y);
        font.drawString(str);
    }

    public static int getFontStringLength(GLBitmapFont font, String str)
    {
        if (str == null) return 0;
        int length = 0;
        char[] chars = str.toCharArray();
        for (int n = 0; n < chars.length; n++)
            length += font.getCharData(chars[n]).advance;
        return length;
    }

    public static void enableTransparentOverlay(GL gl)
    {
        gl.glEnable(GL.GL_ALPHA_TEST);
        gl.glAlphaFunc(GL.GL_GREATER, 0.1f);
        gl.glDepthFunc(GL.GL_LEQUAL);
    }

    public static void disableTransparentOverlay(GL gl)
    {
        gl.glDisable(GL.GL_ALPHA_TEST);
        gl.glAlphaFunc(GL.GL_ALWAYS, 0.0f);
        gl.glDepthFunc(GL.GL_LESS);
    }

    public static void enableLineAntiAliasing(GL gl)
    {
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_DONT_CARE);
    }

    public static void disableLineAntiAliasing(GL gl)
    {
        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
    }
}

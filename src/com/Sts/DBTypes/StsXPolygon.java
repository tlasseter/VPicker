package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Actions.Crossplot.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.*;

public class StsXPolygon extends StsMainObject implements StsSelectable, IPolygon
{
    protected double[][] points; /** allocated array of points: [nAllocated][] */
    protected int nPoints = 0;  /** number of points in polygon */
    protected StsType libraryType;
    protected int ID = -1;

    transient boolean closed = true;
    transient boolean drawingConcave = false;
    transient StsBoundingBox boundingBox = null;

    static final int inc = 10; /** points array increment when resizing */
    static StsConcavePolygon concavePolygon = null;
    static final boolean debug = false;
    static final float[] normal = new float[] { 0, 0, 1 };
    public StsXPolygon()
    {
    }

    public StsXPolygon(StsType type)
    {
        this(type, false);
		this.addToModel();
    }

    public StsXPolygon(StsType type, boolean persistent)
    {
        super(persistent);
        closed = false;
        this.libraryType = type;
        points = new double[inc][];
        if(persistent) ID = getIndex();
    }

    public boolean initialize(StsModel model)
    {
        closed = true;
        return true;
    }

    public void setID(int ID) { this.ID = ID; }

    public void addPoint(StsPoint point, boolean debug)
    {
        if(nPoints == points.length)
        {
            double[][] newPoints = new double[points.length+inc][];
            System.arraycopy(points, 0, newPoints, 0, nPoints);
            points = newPoints;
        }
        points[nPoints++] = new double[] { point.v[0], point.v[1], point.v[2] };
        if(debug) System.out.println("Point added: " + points[nPoints-1][0] + " " + points[nPoints-1][1] + " " + points[nPoints-1][2]);

    }

    public void moveLastPoint(StsPoint point, boolean debug)
    {
        points[nPoints-1][0] = point.v[0];
        points[nPoints-1][1] = point.v[1];
        points[nPoints-1][2] = point.v[2];
        if(debug) System.out.println("Last point moved.");
    }

    public void movePoint(int index, StsPoint point, boolean debug)
    {
        points[index][0] = point.v[0];
        points[index][1] = point.v[1];
        points[index][2] = point.v[2];
        if(debug) System.out.println("Selected point moved.");

    }

	public void movePoint(int index, double[] point, boolean debug)
	{
		points[index][0] = point[0];
		points[index][1] = point[1];
		points[index][2] = point[2];
		if(debug) System.out.println("Selected point moved.");
    }

    public void deleteLastPoint()
    {
        if(nPoints > 0) nPoints--;
    }
/*
    public void removePoint(StsXPolygonVertex vertex)
    {
        StsMessageFiles.infoMessage("Remove point " + vertex + " from Polygon.");
        return;
    }
*/

    public StsType getStsType() { return libraryType; }
    public double[][] getPoints() { return points; }
    public int getNPoints() { return nPoints; }
    public StsColor getStsColor()
    {
        if(libraryType == null) return null;
        return libraryType.getStsColor();
    }
    public Color getColor()
    {
        if(libraryType == null) return null;
        return libraryType.getStsColor().getColor();
    }

    public double[] getPoint(int index) { return points[index]; }

    public void setLastPoint(double[] point) { points[nPoints-1] = point; }

    public int insertPoint(int index, StsPoint point, boolean debug)
    {
        if(nPoints+1 > points.length)
        {
            double[][] newPoints = new double[nPoints+1][];
            System.arraycopy(points, 0, newPoints, 0, nPoints);
            points = newPoints;
        }
        nPoints++;
        for(int n = nPoints-1; n > index+1; n--)
            points[n] = points[n-1];
        points[index+1] = new double[] { point.v[0], point.v[1], point.v[2] };

        if(debug) System.out.println("Point inserted.");

        boundingBox = null;

        return index+1;
    }

    public void deletePoint(int index)
    {
        for(int n = index; n < nPoints-1; n++)
            points[n] = points[n+1];
        nPoints--;
        boundingBox = null;
    }
/*
    public void addPoint(double[] point)
    {
        if(nPoints == points.length)
        {
            double[][] newPoints = new double[points.length+inc][];
            System.arraycopy(points, 0, newPoints, 0, nPoints);
            points = newPoints;
        }
        points[nPoints++] = point;
    }
*/
    public void close()
    {
        closed = true;
        points = (double[][])StsMath.trimArray(points, nPoints);
    }
/*
    public void draw(GL gl)
    {
        StsColor.setGLColor(gl, color);

        if(closed)
            gl.glBegin(GL.GL_LINE_LOOP);
        else
            gl.glBegin(GL.GL_LINE_STRIP);

            for(int n = 0; n < nPoints; n++)
                gl.glVertex3fv(points[n]);
        gl.glEnd();
    }
*/
    public void pick(GL gl, StsGLPanel glPanel)
    {
        gl.glPushName(StsXPolygonAction.TYPE_VERTEX);
        for (int n = 0; n < nPoints; n++)
        {
            gl.glPushName(n);
            StsGLDraw.drawPoint2d(points[n], StsColor.WHITE, gl, 4);
            gl.glPopName();
        }
        gl.glPopName();

        gl.glPushName(StsXPolygonAction.TYPE_EDGE);
        for (int n = 0; n < nPoints; n++)
        {
            gl.glPushName(n);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex2dv(points[n], 0);
            gl.glVertex2dv(points[(n+1)%nPoints], 0);
            gl.glEnd();
            gl.glPopName();
        }
        gl.glPopName();
    }

    public void draw(GL gl, GLU glu, boolean filled)
    {
        float edgeWidth;

        int mode = StsXPolygonAction.getPolygonActionMode(this);
        if(mode != StsXPolygonAction.NONE)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        StsColor color = libraryType.getStsColor();
        if(color == null)
            return;
        if(!filled || !closed || !drawConcave(gl, glu))
            drawDottedEdge(gl, color, edgeWidth, points, nPoints, closed);

        if(mode == StsXPolygonAction.NONE) return;
        {
            int type = StsXPolygonAction.getPolygonActionType(this);

            for(int n = 0; n < nPoints-1; n++)
                StsGLDraw.drawPoint2d(points[n], StsColor.WHITE, gl, 4);

            if(StsXPolygonAction.onFirstPick)
            {
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.WHITE, gl, 8);
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.BLACK, gl, 4);
            }
            else if(type == StsXPolygonAction.TYPE_VERTEX)
            {
                int pickIndex = StsXPolygonAction.pickIndex;
                StsGLDraw.drawPoint(points[pickIndex], StsColor.BLACK, gl, 8);
                StsGLDraw.drawPoint(points[pickIndex], StsColor.WHITE, gl, 6);

                if(pickIndex != nPoints-1)
                    StsGLDraw.drawPoint(points[nPoints-1], StsColor.WHITE, gl, 4);
            }
            else if(mode == StsXPolygonAction.MODE_CREATE)
            {
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.BLACK, gl, 8);
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.WHITE, gl, 6);
            }
            else
            {
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.BLACK, gl, 6);
                StsGLDraw.drawPoint(points[nPoints-1], StsColor.WHITE, gl, 4);
            }
        }
    }

    // Move these routines to GLDraw for general use after debugging
    private static void drawDottedEdge(GL gl, StsColor color, float edgeWidth, double[][] points, int nPoints, boolean closed)
	{
        gl.glLineWidth(edgeWidth);

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
       	drawLineStrip(gl, color, points, nPoints, closed);

        gl.glLineStipple(1, StsGraphicParameters.altDottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
       	drawLineStrip(gl, StsColor.BLACK, points, nPoints, closed);

        gl.glDisable(GL.GL_LINE_STIPPLE);
//        drawLineStrip(gl, color, points);
    }

	private static void drawLineStrip(GL gl, StsColor color, double[][] points, int nPoints, boolean closed)
    {
        if(gl == null || color == null || points == null) return;

       /** IMPORTANT: turn LIGHTING off and then back on for lines */

        try
        {
        if(nPoints < 2) return;

            gl.glDisable(GL.GL_LIGHTING);
            color.setGLColor(gl);

            if(closed)
                gl.glBegin(GL.GL_LINE_LOOP);
            else
                gl.glBegin(GL.GL_LINE_STRIP);

            for(int n = 0; n < nPoints; n++)
                gl.glVertex2dv(points[n], 0);
        }
        catch(Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
	}

    private boolean drawConcave(GL gl, GLU glu)
    {
        libraryType.getStsColor().setGLColor(gl);
        drawingConcave = true;
		if(concavePolygon == null) concavePolygon = new StsConcavePolygon(gl, glu);
		GLUtessellator tesselator = glu.gluNewTess();
		concavePolygon.initialize(tesselator, this, debug, true, normal, name);
        glu.gluTessBeginPolygon(tesselator, null);
        glu.gluTessNormal(tesselator, 0, 0, 1);

        for(int n = 0; n < nPoints; n++)
        {
//            float[] pn = pntNrmls[n];
            double[] vertex = new double[] { points[n][0], points[n][1] };
            glu.gluTessVertex(tesselator, vertex, 0, points[n]);
        }

        glu.gluTessEndPolygon(tesselator);
        return drawingConcave;
    }

    public String getLabel()
    {
        return new String("name: " + getName() + " nPoints: " + nPoints);
    }

    public void drawConcaveFailed(int error)
    {
        drawingConcave = false;
    }

	public void drawConcaveFailed(String errorString)
    {
        drawingConcave = false;
    }

    public void checkDeletePolygon(StsObjectRefList polygons)
    {
        if(closed) return;
        this.delete();
        polygons.delete(this);
    }

    public boolean isClosed() { return closed; }

    public int getID() { return ID; }

    public StsBoundingBox getBoundingBox()
    {
        if(boundingBox == null)
        {
            boundingBox = new StsBoundingBox(false);
            boundingBox.addPoints(points);
        }
        return boundingBox;
    }
}
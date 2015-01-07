//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;

public class StsConcavePolygon extends GLUtessellatorCallbackAdapter
{
	GL gl;
	GLU glu;
	IPolygon polygon;
	boolean isPlanar;
	float[] normal;
	boolean debug = false;
    String name;
    static final boolean debugAll = false;

    public StsConcavePolygon(GL gl, GLU glu)
	{
		this.gl = gl;
		this.glu = glu;
	}

	public void initialize(GLUtessellator tesselator, IPolygon polygon, boolean debug, boolean isPlanar, float[] normal, String name)
	{
		this.polygon = polygon;
		this.debug = debug || debugAll;
		this.isPlanar = isPlanar;
		this.normal = normal;
        this.name = name;
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_BEGIN, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_END, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_VERTEX, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_VERTEX_DATA, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_COMBINE, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_EDGE_FLAG, this);
		glu.gluTessCallback(tesselator, GLU.GLU_TESS_ERROR, this);
        glu.gluTessNormal(tesselator, (double)normal[0], (double)normal[1], (double)normal[2]);
        if (this.debug) debugMessage("initialize", "Polygon " + name);
    }

    public void setNormal(float[] normal)
	{
        this.normal = normal;
	}

	/** Called by the begin callback */
	public void begin(int mode)
	{
		gl.glBegin(mode);
		if (debug) debugMessage("begin", "Polygon " + name + "\n Concave Polygon mode: " + getModeString(mode));
		if(isPlanar) gl.glNormal3fv(normal, 0);
	}

    public void debugMessage(String methodName, String message)
    {
        StsException.systemDebug(this, methodName, message);
    }

    public void vertex(Object dataObject)
	{
		try
		{
			if (dataObject == null)return;
			if (dataObject instanceof float[])
			{
				float[] data = (float[]) dataObject;
				gl.glVertex3fv(data, 0);
				if (debug) debugMessage("vertex", " Concave poly vertex: " + data[0] + " " + data[1] + " " + data[2]);
			}
			else
			{
				double[] data = (double[]) dataObject;
				gl.glVertex3dv(data, 0);
				if (debug) debugMessage("vertex", " Concave Poly vertex: " + data[0] + " " + data[1] + " " + data[2]);
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsConcavePolygon.vertex() failed.", e, StsException.WARNING);
		}
	}

	public void vertexData(Object coordObject, Object vertexObject)
	{
		try
		{
			if (coordObject == null)return;
			double[] coords = (double[]) coordObject;
		/*
			if (isPlanar)
			{
				gl.glNormal3fv(normal);
				debugMessage(" normal: " + normal[0] + " " + normal[1] + " " + normal[2]);
			}
			else
		*/
	        if(!isPlanar)
			{
				gl.glNormal3d(coords[3], coords[4], coords[5]);
//				debugMessage(" normal: " + coords[3] + " " + coords[4] + " " + coords[5]);
			}
			gl.glVertex3dv(coords, 0);
            if (debug) debugMessage("vertexData", " vertex: " + coords[0] + " " + coords[1] + " " + coords[2]);
		}
		catch (Exception e)
		{
			StsException.outputException("StsConcavePolygon.vertex() failed.", e, StsException.WARNING);
		}
	}

	/** Called by the combine callback */
	public void combine(double[] coords, Object[] data, float[] weight, Object[] outData)
	{
        if(coords == null || data == null || weight == null || outData == null) return;
        
        double[] vertex;
        try
		{
			if(data[0] instanceof double[])
			{
                int nPoints = data.length;
                int nElements = ((double[])data[0]).length;
                vertex = new double[nElements];
                vertex[0] = coords[0];
			    vertex[1] = coords[1];
			    vertex[2] = coords[2];
                for(int n = 0; n < nPoints; n++)
                {
                    double[] values = (double[])data[n];
                    if(values == null) continue;
                    double wt = weight[n];
                    if(wt == 0.0) continue;
                    for(int i = 3; i < nElements; i++)
                        vertex[i] += wt * values[i];
                }
			}
			else
			{
				int nPoints = data.length;
                int nElements = ((float[])data[0]).length;
                vertex = new double[nElements];
			    vertex[0] = coords[0];
			    vertex[1] = coords[1];
			    vertex[2] = coords[2];
                for(int n = 0; n < nPoints; n++)
                {
                    float[] values = (float[])data[n];
                    if(values == null) continue;
                    double wt = weight[n];
                    if(wt == 0.0) continue;
                    for(int i = 3; i < nElements; i++)
                        vertex[i] += wt * values[i];
                }
			}
			outData[0] = vertex;
            if (debug) debugMessage("combine", " vertex: " + coords[0] + " " + coords[1] + " " + coords[2]);
        }
		catch(Exception e)
		{
			StsException.systemError("StsConcavePolygon.combine() failed.");
		}
	}

	/** Called by the end callback */
	public void end()
	{
		gl.glEnd();
		if (debug) debugMessage("end", " polygon end.");
	}

	/** Called by the error callback */
	public void error(int error)
	{
		String errorString = glu.gluErrorString(error);
		StsException.systemError("StsConcavePolygon GLU error for polygon: " + polygon.getLabel() + " error: " +
								 errorString);
		polygon.drawConcaveFailed(error);
	}

	/** Called by the edge flag callback */
	public void edgeFlag(boolean flag)
	{
		gl.glEdgeFlag(flag ? true : false);
		if (debug) debugMessage("edgeFlag", " edgeFlag set: " + flag);
	}

	private String getModeString(int mode)
	{
		switch (mode)
		{
			case 4:
				return new String("GL_TRIANGLES");
			case 5:
				return new String("GL_TRIANGLE_STRIP");
			case 6:
				return new String("GL_TRIANGLE_FAN");
			default:
				return new String("UNKNOWN MODE");
		}
	}

	private String getErrorString(int error)
	{
		switch (error)
		{
			case 100109:
				return new String("Data error.");
			case 100151:
				return new String("Missing gluBeginPolygon.");
			case 100152:
				return new String("Missing gluBeginContour.");
			case 100153:
				return new String("Missing gluEndPolgyon.");
			case 100154:
				return new String("Missing gluEndContour.");
			case 100155:
				return new String("Misoriented or self-intersecting loops, coordinate too large.");
			case 100156:
				return new String("Coincident vertices: need combine callback.");
			case 100157:
				return new String("Colinear vertices.");
			case 100158:
				return new String("Intersecting edges.");
			default:
				return new String("Unknown error: " + error);
		}
	}
}

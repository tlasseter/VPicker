package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsFractureDisplayable;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.Utilities.StsConcavePolygon;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;

import javax.media.opengl.GL;
import javax.media.opengl.glu.GLU;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 8:17:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsGolderFracture extends StsSerialize
{
    int fractureNumber;
	int nNodes;
	public float[][] vertices;
	public float[] properties;
    float[] normal;
	private float[] center;
    long time;
    public boolean highlight = false;
    boolean connected = false;
    StsObject correlatedObject = null;
    int objectIndex = -1;

	/** zDomain currently being isVisible. Changing domains requires building new display lists and textures;
	 *  in which case zDomainDisplayed is set to none, display() method deletes display lists, rebuilds displays
	 *  for current project zDomain and sets ZDomainDisplayed to this zDomain.
	 */
	transient protected byte zDomainDisplayed = StsParameters.TD_NONE;
 /** original domain this object was built in. If and when velocity model is rebuilt, this domain is used as the coordinate source. */
    transient protected byte zDomainOriginal = StsParameters.TD_DEPTH;
	transient public int setNumber;
	transient public int nVertices;
	/** Display lists should be used (controlled by View:Display Options) */
	transient boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	transient boolean usingDisplayLists = true;

    public StsGolderFracture()
    {
    }


	static public StsGolderFracture fractureConstructor(String[] tokens)
	{
		return new StsGolderFracture(tokens, 3);
	}

	private StsGolderFracture(String[] tokens, int nSkip)
	{
		setProperties(tokens, nSkip);
	}

	static public StsGolderFracture tessFractureConstructor(String[] tokens, float[][] allVertices)
	{
		return new StsGolderFracture(tokens, allVertices);
	}

	private StsGolderFracture(String[] tokens, float[][] allVertices)
	{
		int nVertex = -1;
		try
		{
			fractureNumber = Integer.parseInt(tokens[0]);
			nNodes = Integer.parseInt(tokens[1]);
			vertices = new float[nNodes][];
			for(int n = 0; n < nNodes; n++)
			{
				nVertex = Integer.parseInt(tokens[n+2]) - 1;
				vertices[n] = allVertices[nVertex];
			}
			computeNormal();
			setProperties(tokens, nNodes + 2);
			zDomainDisplayed = currentModel.getProject().getZDomain();
		}
		catch(Exception e)
		{
			StsException.systemError(this, "constructor", "nNodes " + nNodes + " nVertex " + nVertex + " nTotalVertices " + allVertices.length);
		}
	}

	private void computeNormal()
	{
		int n00, n01, n10, n11;
		int nVertices = vertices.length;
		n00 = 0;
		n01 = nVertices/2;
		n10 = n01/2;
		n11 = (n01 + nVertices)/2;
		float[] v0 = StsMath.subtract(vertices[n01],  vertices[n00]);
		float[] v1 = StsMath.subtract(vertices[n11],  vertices[n10]);
		normal = StsMath.leftCrossProduct(v0, v1);
		StsMath.normalize(normal);
	}

	private void setProperties(String[] tokens, int nSkip)
	{
		int nProperties = tokens.length - nSkip;
		properties = new float[nProperties];
		for(int n = nSkip, i = 0; n < tokens.length; n++, i++)
			properties[i] = Float.parseFloat(tokens[n]);
	}

    public void setHighlight(StsObject object, int index)
    {
        this.highlight = true;
        correlatedObject = object;
        objectIndex = index;
    }

    public void setHighlight(StsFractureDisplayable object, int index)
    {
        this.highlight = true;
        correlatedObject = (StsObject)object;
        objectIndex = index;
    }

    public void setConnection(boolean val)
    {
        this.connected = val;
    }

    public void clearHighlight()
    {
        highlight = false;
        correlatedObject = null;
        objectIndex = -1;
    }

	public void addVertices(float[][] vertices)
	{
		this.vertices = vertices;
		computeCenter();
	}

	private void computeCenter()
	{
		center = StsMath.average(vertices, 3);
	}

	public float[] getCenter()
	{
		if(center != null) return center;
		computeCenter();
		return center;
	}

	public void addNormal(float[] normal)
	{
		this.normal = normal;
		//debugCheckNormal();
	}

    public float getAzimuth()
    {
        double x = normal[1];
		double y = normal[0];
		return (float)StsMath.atan2d(y, x);
    }

	private void debugCheckNormal()
	{
		int n00, n01, n10, n11;
		int nVertices = vertices.length;
		n00 = 0;
		n01 = nVertices/2;
		n10 = n01/2;
		n11 = (n01 + nVertices)/2;
		float[] v0 = StsMath.subtract(vertices[n01],  vertices[n00]);
		float[] v1 = StsMath.subtract(vertices[n11],  vertices[n10]);
		float[] debugNormal = StsMath.leftCrossProduct(v0, v1);
		StsMath.normalize(normal);
		StsMath.normalize(debugNormal);
		float dot = StsMath.dot(normal, debugNormal);
		// StsException.systemDebug(this, "debugCheckNormal", "dot: " + dot);
		if(dot < 0.0f)
			StsMath.scale(normal, -1.0f);
		float azimuth = getAzimuth();
		StsException.systemDebug("normal " + normal[0] +  " " + normal[1]+  " " + normal[2]+  " azimuth " + azimuth);
	}

    public long getTime()
    {
        if(correlatedObject != null)
        {
            if(correlatedObject instanceof StsDynamicSensor)
                return ((StsTimeCurve)((StsDynamicSensor)correlatedObject).getTimeCurves().getElement(0)).getTimeVector().getValue(objectIndex);
        }
        return time;
    }

    public boolean isHighlighted()
    {
        return highlight;
    }
    public boolean isConnected() { return connected; }

	public void display(StsGLPanel3d glPanel3d, boolean displayEdges, StsColor color)
	{
		byte projectZDomain = currentModel.getProject().getZDomain();
		if (projectZDomain != zDomainDisplayed)
		{
            if(!canDraw()) return;
            zDomainDisplayed = projectZDomain;
			checkAddTime();
		}

        GL gl = glPanel3d.getGL();
        color.setGLColor(gl);
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);
        gl.glBegin(GL.GL_POLYGON);
		gl.glNormal3fv(normal, 0);
		for (float[] vertex : vertices)
		{
		    gl.glNormal3fv(normal, 0);
		    gl.glVertex3fv(getXYZorT(vertex), 0);
		}
        gl.glEnd();
	}

    public float[] getXYZorT(float[] xyzmt)
    {
        if (isDepth || xyzmt.length < 5)
            return xyzmt;
        else
            return new float[] { xyzmt[0],  xyzmt[1], xyzmt[4] };
    }

    public int getObjectIndex() { return objectIndex; }
    public StsObject getObject() { return correlatedObject; }
	private void checkAddTime()
	{
		StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
        if(velocityModel == null) return;
		int nVertices = vertices.length;
		for(int n = 0; n < nVertices; n++)
		{
			float[] vertex = vertices[n];
			try
			{
				double t = velocityModel.getT(vertex);
				if(vertex.length < 5)
				{
					float[] newVertex = new float[5];
					System.arraycopy(vertex, 0, newVertex, 0, 3);
					vertices[n] = newVertex;
					vertex = newVertex;
				}
				vertex[4] = (float)t;
			}
			catch(Exception e)
			{
				StsException.outputWarningException(this, "checkAddTime", e);
				return;
			}
		}
	}

   /** section can be drawn if the original domain in which it was constructed is the same as the current zDomain or
     *  a velocityModel exists which can convert from this original zDomain.
     * @return true if it can or can't be drawn for reasons defined above.
     */
    public boolean canDraw()
    {
        if(zDomainOriginal == currentModel.getProject().zDomain) return true;
        return currentModel.getProject().velocityModel != null;
    }

	public void displayX(StsGLPanel3d glPanel3d, boolean displayEdges, StsColor color)
	{
        color = overrideFractureColor(color);
        GL gl = glPanel3d.getGL();
        color.setGLColor(gl);
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);
        gl.glBegin(GL.GL_POLYGON);
		gl.glNormal3fv(normal, 0);
		for (float[] vertex : vertices)
		{
		    gl.glNormal3fv(normal, 0);
		    gl.glVertex3fv(vertex, 0);
		}
        gl.glEnd();
     }

    public StsColor overrideFractureColor(StsColor color)
    {
        if(correlatedObject != null)
        {
            if(correlatedObject instanceof StsDynamicSensor)
            {
                float colorFloat = 0.0f;
			    float[] colorFloats = ((StsDynamicSensor)correlatedObject).getColorFloats();
                if(colorFloats != null)
                    colorFloat = colorFloats[objectIndex];
                return ((StsDynamicSensor)correlatedObject).defineColor(colorFloat);
            }
            else
            {
                return ((StsFracture)correlatedObject).getStsColor();
            }
        }
        return color;
    }

	public float getPropertyValue(int index)  { return properties[index]; }

	public void highlightIfInsideLimits(StsDynamicSensor sensor, int index, float maxDistance, boolean azimuthLimit, float minAzimuth, float maxAzimuth)
	{
		if(center == null) return;
		float fractureDistance = StsMath.distance(sensor.getXYZ(index), center, 3);
		if(azimuthLimit)
		{
			float azimuth = getAzimuth();
			if((azimuth > maxAzimuth) || (azimuth < minAzimuth))
				 return;
		}
		if(fractureDistance < maxDistance)
			setHighlight(sensor, index);
	}

	public void setCenter(float[] center)
	{
		this.center = center;
	}

	/** for all vertex points, find maxZ and project XY coordinates to vertical plane thru center which is orthogonal to normal.
	 *  take the extreme xy range of this points and use to define a line thru the max Z */

	public StsPoint[][] getTopAndBottomLines()
	{
		float minZ = StsParameters.largeFloat;
		float maxZ = -StsParameters.largeFloat;
		float dotMin = StsParameters.largeFloat;
		float dotMax = -StsParameters.largeFloat;
		float lineX = normal[1];
		float lineY = -normal[0];
		getCenter();
		for(float[] vertex : vertices)
		{
			minZ = Math.min(minZ, vertex[2]);
			maxZ = Math.max(maxZ, vertex[2]);
			float dx = vertex[0] - center[0];
			float dy = vertex[1] - center[1];

			float dot = dx*lineX + dy*lineY;
			if(dot < dotMin) dotMin = dot;
			if(dot > dotMax) dotMax = dot;
		}
		float xMin = center[0] + dotMin*lineX;
		float yMin = center[1] + dotMin*lineY;
		float xMax = center[0] + dotMax*lineX;
		float yMax = center[1] + dotMax*lineY;
		StsPoint topMinPoint = new StsPoint(xMin, yMin, minZ);
		StsPoint topMaxPoint = new StsPoint(xMax, yMax, minZ);
		StsPoint botMinPoint = new StsPoint(xMin, yMin, maxZ);
		StsPoint botMaxPoint = new StsPoint(xMax, yMax, maxZ);
		return new StsPoint[][] { { topMinPoint, topMaxPoint }, { botMinPoint, botMaxPoint} };
	}
}

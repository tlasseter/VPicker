package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 8:17:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsMVFracture extends StsSerialize
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 6768992662538283926L;
	static final int UP_RIGHT = 0;
    static final int UP_LEFT = 1;
    static final int LOW_LEFT = 2;
    static final int LOW_RIGHT = 3;

    StsMVFractureSet fractureSet;

	/** Display lists should be used (controlled by View:Display Options) */
	transient boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	transient boolean usingDisplayLists = true;
    
    float[] center;
    float[] heightAxis;
    float[] lengthAxis;
    float[] normal;
  
	float aperture;
    float length;
    float azimuth;
    float dip;
    float aspectRatio;
    long time;

    float[] attributes;
    String[] attributeNames;

    public static final int NONE = 0;
    public static final int APERTURE = 1;
    public static final int LENGTH = 2;
    public static final int DIP = 3;
    public static final int AZIMUTH = 4;
    public static final int ASPECTRATIO = 5;
    public static int[] staticAtts = new int[] { NONE, APERTURE, LENGTH, DIP, AZIMUTH, ASPECTRATIO };
    public static String S_NONE = "None";
    public static String S_APERTURE = "Aperture";
    public static String S_LENGTH = "Length";
    public static String S_DIP = "Dip";
    public static String S_AZIMUTH = "Azimuth";
    public static String S_ASPECTRATIO = "Aspect Ratio";

    public StsMVFracture()
    {
    }

    public StsMVFracture(StsMVFractureSet fractureSet, float[] center, float azimuth, float dip, float length, float aspectRatio, float aperture, long time, float[] attributes, String[] attributeNames )
	{
        this(fractureSet, center, azimuth, dip, length, aspectRatio, aperture);
        this.time = time;
        this.attributes = attributes;
        this.attributeNames = attributeNames;
    }

    public StsMVFracture(StsMVFractureSet fractureSet, float[] center, float azimuth, float dip, float length, float aspectRatio, float aperture )
	{
        this.fractureSet = fractureSet;
        this.center = center;
        this.length = length;
        this.aperture = aperture;
        this.azimuth = azimuth;
        this.dip = dip;
        this.aspectRatio = aspectRatio;
        initialize(azimuth, dip, aspectRatio );
    }
    
    /**
    *    private initialization function that initializes values
    *
    */  
    void initialize( float azimuth, float dip, float aspectRatio )
    {
     	azimuth = (float)Math.toRadians((float)azimuth);
    	dip = (float)Math.toRadians((float)dip);
    	normal = new float[] { (float) (Math.sin(azimuth) * Math.sin(dip)), 
    			(float) (Math.cos(azimuth) * Math.sin(dip)),
    			(float)  Math.cos(dip) };


        // The minor axis is the y direction rotated anti clockwise by az degrees times the 
        // the width
    	lengthAxis = new float[] { -(float)Math.cos(azimuth),  (float)Math.sin(azimuth) , 0.0f };
    	
    	lengthAxis[0] *= (length/2.0);
    	lengthAxis[1] *= (length/2.0);
    	lengthAxis[2] *= (length/2.0);
        

        /// The major axis is at right angles to both the minor axis and the normal - multiplied by length
    	heightAxis = StsMath.crossProduct(lengthAxis, normal );
    	StsMath.normalizeVector(heightAxis);
        float height = length * aspectRatio;
        heightAxis[0] *= (height/2.0);
        heightAxis[1] *= (height/2.0);
        heightAxis[2] *= (height/2.0);
    }

    /**
     *
     */
    public String[] getAttributeNames() { return attributeNames; };
    public long getTime() { return time; }
    /**
    *
    * @return the area of this fracture
    *
    */
    public double area()
    {
    	double lengthAxisLength = Math.sqrt(lengthAxis[0]*lengthAxis[0] + lengthAxis[1]*lengthAxis[1] + lengthAxis[2]*lengthAxis[2]);
    	double heightAxisLength = Math.sqrt(heightAxis[0]*heightAxis[0] + heightAxis[1]*heightAxis[1] + heightAxis[2]*heightAxis[2]);
        return 4.0 * lengthAxisLength * heightAxisLength;
    }

    /**
    * @return a corner of the fracture 
    * 
    *
    */
    public float[]  getCorner(int index)
    {
    	float[] tmpVec = new float[3];

    	switch(index)
    	{
    	case UP_RIGHT:						//upper right
    		tmpVec[0] =  center[0] + heightAxis[0] * fractureSet.getVerticalScale() + lengthAxis[0] ;
    		tmpVec[1] =  center[1] + heightAxis[1] * fractureSet.getVerticalScale()+ lengthAxis[1] ;
    		tmpVec[2] =  center[2] + heightAxis[2] * fractureSet.getVerticalScale()+ lengthAxis[2] ;
    		break;
    	case UP_LEFT:						//upper left
    		tmpVec[0] =  center[0] + heightAxis[0] * fractureSet.getVerticalScale() - lengthAxis[0] ;
    		tmpVec[1] =  center[1] + heightAxis[1] * fractureSet.getVerticalScale() - lengthAxis[1] ;
    		tmpVec[2] =  center[2] + heightAxis[2] * fractureSet.getVerticalScale() - lengthAxis[2] ;
    		break;
    	case LOW_LEFT:						// lower left
    		tmpVec[0] =  center[0] - heightAxis[0] * fractureSet.getVerticalScale() - lengthAxis[0] ;
    		tmpVec[1] =  center[1] - heightAxis[1] * fractureSet.getVerticalScale() - lengthAxis[1] ;
    		tmpVec[2] =  center[2] - heightAxis[2] * fractureSet.getVerticalScale() - lengthAxis[2] ;
    		break;
    	case LOW_RIGHT:						// lower right
    		tmpVec[0] =  center[0] - heightAxis[0] * fractureSet.getVerticalScale() + lengthAxis[0] ;
    		tmpVec[1] =  center[1] - heightAxis[1] * fractureSet.getVerticalScale() + lengthAxis[1] ;
    		tmpVec[2] =  center[2] - heightAxis[2] * fractureSet.getVerticalScale() + lengthAxis[2] ;
    		break;
    	}
    	return tmpVec;
    }

        
	public void addToSet()
	{
		fractureSet.addFracture(this);
	}
	
    
	public void display(StsGLPanel3d glPanel3d, boolean displaySectionEdges, StsColor color)
	{
        float[][] xy = new float[4][3];
        xy[0] = getCorner(LOW_LEFT);
        xy[1] = getCorner(UP_LEFT);
        xy[2] = getCorner(LOW_RIGHT);
        xy[3] = getCorner(UP_RIGHT);

		GL gl = glPanel3d.getGL();

        color.setGLColor(gl);
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);

        gl.glBegin(GL.GL_QUAD_STRIP);
        gl.glNormal3fv(normal, 0);

        gl.glVertex3fv(xy[0], 0);
        gl.glVertex3fv(xy[1], 0);

        gl.glVertex3fv(xy[2], 0);
        gl.glVertex3fv(xy[3], 0);
        gl.glEnd();
        if(displaySectionEdges)
        {
            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE);
            StsColor.BLACK.setGLColor(gl);
            gl.glLineWidth(2.0f);
            gl.glBegin(gl.GL_POLYGON);
        
            gl.glVertex3fv(xy[0], 0);
            gl.glVertex3fv(xy[1], 0);

            gl.glVertex3fv(xy[3], 0);
            gl.glVertex3fv(xy[2], 0);

            gl.glVertex3fv(xy[0], 0);
            gl.glEnd();
        }
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);      // Reset to default or sensors and others are not filled.
    }

    public float getAttributeValue(int attIdx)
    {
        if(attIdx > staticAtts.length -1)
             return attributes[attIdx-staticAtts.length]; 
        switch(attIdx)
         {
             case NONE:
                 return -1;
             case APERTURE:
                 return aperture;
             case LENGTH:
                 return length;
             case DIP:
                 return dip;
             case AZIMUTH:
                 return azimuth;
             case ASPECTRATIO:
                 return aspectRatio;
             default:
                 return -1;
         }
    }
}


//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTimeEnabled;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Sounds.StsSound;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;

public class StsPerforationMarker extends StsWellMarker implements StsSelectable
{
    // instance fields
    protected float length = 1.0f;
    protected boolean highlighted = false;
    protected int numShots = 1;

    transient boolean firstAlive = false;
    transient boolean accentPerf = false;
    transient boolean initialized = false;

    /** DB constructor */
    public StsPerforationMarker()
    {
    }

    /** constructor */
    private StsPerforationMarker(String name, StsWell well, byte type, StsPoint location, float length, int nShots, long time) throws StsException
    {
        if (well==null)
            throw new StsException(StsException.WARNING, "StsPerforationMarker:" + " Cannot create a marker for a null well.");
        if (!StsMarker.isTypeValid(type))
            throw new StsException(StsException.WARNING, "StsPerforationMarker:" + " Marker type is invalid.");

        marker = new StsMarker(name, StsMarker.PERFORATION, StsColor.colors32[0]);
        marker.addWellMarker(this);
	    this.well = well;
	    setType(type);
	    setLocation(location);
	    setLength(length);
        setNumShots(nShots);
        setBornDate(time);
        well.addMarker(this);	    
    }
    
    static public StsPerforationMarker constructor(String name, StsWell well, byte type, float mdepth, float length)
    {
        try
        {
            StsPoint location = well.getPointAtMDepth(mdepth, false);
            return constructor(name, well, type, location, length);
        }
        catch(Exception e)
        {
            StsMessageFiles.errorMessage(e.getMessage());
            return null;
        }
    }
    
    static public StsPerforationMarker constructor(String name, StsWell well, byte type, StsPoint location, float length)
    {
        try
        {
            return new StsPerforationMarker(name, well, type, location, length, 1, 0l);
        }
        catch(Exception e)
        {
            StsMessageFiles.errorMessage(e.getMessage());
            return null;
        }
    }

    static public StsPerforationMarker constructor(String name, StsWell well, byte type, StsPoint location, float length, int nShots, long time)
    {
        try
        {
            return new StsPerforationMarker(name, well, type, location, length, nShots, time);
        }
        catch(Exception e)
        {
            StsMessageFiles.errorMessage(e.getMessage());
            return null;
        }
    }
	private void initializeColor()
	{
		StsColor color = StsColor.colors32[colorIndex++];
        setStsColor(color);
	}
	
    public float getLength() { return length; }
    public void setLength(float len)
    {
    	length = len;
    }

    public int getNumShots() { return numShots; }
    public void setNumShots(int nshots)
    {
    	numShots = nshots;
    }

    public void setHighlighted(boolean highlight)
    {
    	highlighted = highlight;
        dbFieldChanged("highlighted", highlighted);        
        currentModel.viewObjectRepaint(well, well);
    }

    public boolean delete()
    {
        well.getMarkers().delete(this);
        marker.getWellMarkers().delete(this);
        return super.delete();
    }
    public void display2d(StsGLPanel3d glPanel3d, int dirNo, boolean displayName)
    {
            display2d(glPanel3d, dirNo, displayName, getStsColor(), 1.0);
    }     
    public void display2d(StsGLPanel3d glPanel3d, int dirNo, boolean displayName, StsColor color, double viewshift)
    {
        boolean isAlive = isAlive(currentModel.getProject().getProjectTime());
        if(!well.isVisible || !well.getDrawPerfMarkers() || !well.getWellClass().getDisplayPerfMarkers())
            return;

        accentPerf = false;
        if(well.getWellClass().getEnableTime())
        {
            if(!isAlive)
            {
                firstAlive = false;
                return;
            }

            if((firstAlive == false) && isAlive)
            {
                firstAlive = true;
                accentPerf = true;
            }
        }

        float[] xyz = location.getXYZorT();
        float[] xy = new float[2];
        float[][] points = new float[2][];
        int xidx = 0, yidx = 1;
        switch(dirNo)
        {
        	case StsCursor3d.ZDIR:
        		xy[0] = xyz[0];
        		xy[1] = xyz[1];
        		xidx = 0;
        		yidx = 1;
        		break;
        	case StsCursor3d.YDIR:
        		xy[0] = xyz[0];
        		xy[1] = xyz[2];
        		xidx = 0;
        		yidx = 2;
        		break;
        	case StsCursor3d.XDIR:
        		xy[0] = xyz[1];
        		xy[1] = xyz[2];
        		xidx = 1;
        		yidx = 2;
        		break;         		
        }
        final StsColor currentColor = StsWell.getWellClass().getDefaultColor(getIndex());
    	StsColor.setGLJavaColor(glPanel3d.getGL(), currentColor.getColor());
	    float mdepth = location.getM();
	    
	    //System.out.println("Mdepth=" + mdepth + " Length=" + length);
	    points[0] = well.getPointAtMDepth(mdepth-length/2, true).getXYZorT();
	    points[1] = well.getPointAtMDepth(mdepth+length/2, true).getXYZorT();


        GL gl = glPanel3d.getGL();
        //StsPoint[] stsPts = new StsPoint[2];
        //stsPts[0] = new StsPoint(points[0]);
        //stsPts[1] = new StsPoint(points[1]);
        StsColor aColor = new StsColor(currentColor, 0.5f);
        //StsGLDraw.drawLineStrip2d(gl, aColor, stsPts, 20);
        int scale = (int)well.getWellClass().getPerfScale();
        if(highlighted)
        {
            StsPoint[] stsPoints = new StsPoint[2];
            stsPoints[0] = new StsPoint(points[0]);
            stsPoints[1] = new StsPoint(points[1]);
            StsGLDraw.drawLine2d(gl, currentColor, scale*2, stsPoints, 2);
        }
        else
        {
            StsGLDraw.drawDottedLine2d(gl, currentColor, StsColor.BLACK, scale, points[0][xidx], points[0][yidx], points[1][xidx], points[1][yidx]);
        }
        accent2DPerforation(currentColor, location.getXYZorT());

        if(!displayName) return;
        
    	StsColor.setGLJavaColor(glPanel3d.getGL(), currentColor.getColor());
        // Display the name
        gl.glDisable(GL.GL_LIGHTING);
        StsGLDraw.fontOutput(gl, xy[0]+5, xy[1]+5, marker.getName(), GLHelvetica12BitmapFont.getInstance(gl));         
        gl.glEnable(GL.GL_LIGHTING);    	 
    }

   public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean isDrawingCurtain)
   {
          display(glPanel3d, displayName, isDrawingCurtain, 20, 60, getStsColor(), 1.0);
   }

   public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean isDrawingCurtain, int width, int height,
       StsColor color, double viewshift)
    {
        if (!currentModel.getProject().canDisplayZDomain(well.zDomainSupported))
            return;

        boolean isAlive = isAlive(currentModel.getProject().getProjectTime());
        if(!well.isVisible || !well.getDrawPerfMarkers() || !well.getWellClass().getDisplayPerfMarkers())
            return;

        accentPerf = false;
        if(well.getWellClass().getEnableTime())
        {
            if(!isAlive)
            {
                firstAlive = false;
                return;
            }

            if((firstAlive == false) && isAlive)
            {
                firstAlive = true;
                accentPerf = true;
            }
        }

	    float[][] points = new float[2][];
        if(isDrawingCurtain) viewshift += 2.0;

	    float[] xyz = location.getXYZorT();
	    float mdepth = location.getM();
	    
	    //System.out.println("Mdepth=" + mdepth + " Length=" + length);
        int scale = (int)well.getWellClass().getPerfScale();
        GL gl = glPanel3d.getGL();
        StsColor aColor = new StsColor(color, 0.5f);
        points[0] = well.getPointAtMDepth(mdepth-length/2, true).getXYZorT();
        if(points[0] == null) return;
	    points[1] = well.getPointAtMDepth(mdepth+length/2, true).getXYZorT();
        double[] screenPoint0 = glPanel3d.getScreenCoordinates(points[0]);
        double[] screenPoint1 = glPanel3d.getScreenCoordinates(points[1]);
        double pixelDistance = StsMath.distance(screenPoint0, screenPoint1, 3);

            if(pixelDistance < 2)
            {
                double stretchFactor = 2/pixelDistance;
                screenPoint1 = StsMath.interpolate(screenPoint0, screenPoint1, stretchFactor);
                points[1] = glPanel3d.getWorldCoordinateFloats(screenPoint1);
                if(highlighted)
                    StsGLDraw.drawLineStrip(gl, aColor, points, scale);
                else
                    StsGLDraw.drawDottedLineStrip(gl, aColor, points, scale);
            }
            else
            {
                if(highlighted)
                    StsGLDraw.drawLineStrip(gl, aColor, points, scale);
                else
                    StsGLDraw.drawDottedLineStrip(gl, aColor, points, scale); 
            }

            // Perforation Shots
            float start = mdepth-length/2;
            float step = length/(numShots-1);
            for(int jj=0; jj<numShots; jj++)
            {
                float md = start + (jj*step);
                StsGLDraw.drawSphere(glPanel3d, well.getPointAtMDepth(md, true).getXYZorT(), StsColor.GREY, 5);
            }

        accent3DPerforation(glPanel3d, aColor, mdepth-length/2, scale);
        if(!displayName) return;

        color.setGLColor(gl);
        glPanel3d.setViewShift(gl, viewshift + 2.0);
        gl.glDisable(GL.GL_LIGHTING);
        StsGLDraw.fontHelvetica12(gl, xyz, marker.getName());
        gl.glEnable(GL.GL_LIGHTING);
        glPanel3d.resetViewShift(gl);
    }

    public void accent3DPerforation(StsGLPanel3d glPanel3d, StsColor color, float start, int scale)
    {
        float[] point = null;
        if(accentPerf)
        {
//            if(well.getWellClass().getEnableSound())
//                StsSound.play(well.getWellClass().getDefaultPerfSound());
        }
    }
    public void accent2DPerforation(StsColor color, float[] xyz)
    {
        if(accentPerf)
        {
            ;
        }
    }

    public StsColor getStsColor()
    {
        return StsWell.getWellClass().getDefaultColor(getIndex());
    }
}














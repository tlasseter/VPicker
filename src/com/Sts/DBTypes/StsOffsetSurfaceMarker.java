
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsProject;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.Types.StsPoint;
import com.Sts.UI.StsSelectable;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsGLDraw;
import com.Sts.Utilities.StsMath;
import com.magician.fonts.GLHelvetica12BitmapFont;

import javax.media.opengl.GL;


public class StsOffsetSurfaceMarker extends StsWellMarker implements StsSelectable
{
    // instance fields
    protected StsPoint offsetLocation = null;

    /** DB constructor */
    public StsOffsetSurfaceMarker()
    {
    }

    /** constructor */

    public StsOffsetSurfaceMarker(StsWell well, StsMarker marker, float mdepth, float zoffset) throws StsException
    {
        if (well==null)
        {
            throw new StsException(StsException.WARNING, "StsOffsetSurfaceMarker.StsOffsetSurfaceMarker: Cannot create a marker for a null well.");
        }
        if (marker==null)
        {
            throw new StsException(StsException.WARNING, "StsOffsetSurfaceMarker.StsOffsetSurfaceMarker:Marker is null.");
        }

        this.well = well;
        this.marker = marker;
        this.setName(marker.getName());

        location = well.getPointAtMDepth(mdepth, false);
        float md = well.getMDepthFromDepth(location.getZ() + zoffset);
        float time = well.getTimeFromMDepth(md);
        offsetLocation = location.copy();
        offsetLocation.setZ(offsetLocation.getZ() + zoffset);
        offsetLocation.setT(time);

        marker.addWellMarker(this);
        well.addMarker(this);
    }


     static public StsOffsetSurfaceMarker constructor(StsWell well, StsMarker marker, float mdepth, float zoffset)
     {
         try
         {
             return new StsOffsetSurfaceMarker(well, marker, mdepth, zoffset);
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

    public StsPoint getOffsetLocation() { return offsetLocation; }
    public void setOffsetLocation(StsPoint pt)
    {
    	offsetLocation = pt;
    }
    public boolean delete()
    {
        well.getMarkers().delete(this);
        marker.getWellMarkers().delete(this);
        return super.delete();
    }
     public void display2d(StsGLPanel3d glPanel3d, int dirNo, boolean displayName, boolean drawDifferent)
     {
         if(drawDifferent)
         {
             display2d(glPanel3d, dirNo, false, 8, 7, StsColor.BLACK, 2.0);
             display2d(glPanel3d, dirNo, displayName, 6, 3, getStsColor(), 3.0);
         }
         else
             display2d(glPanel3d, dirNo, displayName, 6, 3, getStsColor(), 1.0);
     }
     public void display2d(StsGLPanel3d glPanel3d, int dirNo, boolean displayName, int width, int height, StsColor color, double viewshift)
     {
         float[] xyz = offsetLocation.getXYZorT();
         float[] xy = new float[2];
         switch(dirNo)
         {
         	case StsCursor3d.ZDIR:
         		xy[0] = xyz[0];
         		xy[1] = xyz[1];
         		break;
         	case StsCursor3d.YDIR:
         		xy[0] = xyz[0];
         		xy[1] = xyz[2];
         		break;
         	case StsCursor3d.XDIR:
         		xy[0] = xyz[1];
         		xy[1] = xyz[2];
         		break;
         }
     	 StsColor.setGLJavaColor(glPanel3d.getGL(), color.getColor());
         StsGLDraw.drawPoint2d(xy, glPanel3d.getGL(), width);
         if(!displayName) return;
         GL gl = glPanel3d.getGL();

         // Display the name
         gl.glDisable(GL.GL_LIGHTING);
         StsGLDraw.fontOutput(gl, xy[0]+width, xy[1]+width, getName(), GLHelvetica12BitmapFont.getInstance(gl));
         gl.glEnable(GL.GL_LIGHTING);
     }

     public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean isDrawingCurtain, StsColor color)
     {
            display(glPanel3d, displayName, isDrawingCurtain, 10, 3, color, 1.0);
     }
     public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean isDrawingCurtain, boolean drawDifferent)
    {
        if(drawDifferent)
        {
            display(glPanel3d, displayName, isDrawingCurtain, 14, 7, StsColor.BLACK, 2.0);
            display(glPanel3d, displayName, isDrawingCurtain, 10, 3, getStsColor(), 3.0);
        }
        else
            display(glPanel3d, displayName, isDrawingCurtain, 10, 3, getStsColor(), 1.0);
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean isDrawingCurtain, int width, int height,
        StsColor color, double viewshift)
     {
         if(isDrawingCurtain) viewshift += 2.0;
         float[] xyz = new float[] {location.getX(), location.getY(), location.getZorT()};
         StsPoint first = new StsPoint(xyz);

         xyz = new float[] {offsetLocation.getX(), offsetLocation.getY(), offsetLocation.getZorT()};
         StsPoint last = new StsPoint(xyz);
         StsPoint[] points = new StsPoint[] {first, last};

         GL gl = glPanel3d.getGL();
         StsGLDraw.drawPoint(xyz, color, glPanel3d, width, height, viewshift);
         StsGLDraw.drawDottedLine(gl, color, false, points, 2);
         if(!displayName) return;
         glPanel3d.setViewShift(gl, viewshift + 2.0);
         gl.glDisable(GL.GL_LIGHTING);
         StsGLDraw.fontHelvetica12(glPanel3d.getGL(), xyz, marker.getName());
         gl.glEnable(GL.GL_LIGHTING);
         glPanel3d.resetViewShift(gl);
     }

     // Override to include the zOffset
     public StsPoint getLocation()
     {
         return offsetLocation;
     }
     public float getZ() { return offsetLocation.getZ(); }
}

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;

public class StsWellMarker extends StsMainTimeObject implements StsSelectable
 {
     // display fields
     static public StsFieldBean[] displayFields = null;
     static public StsFieldBean[] propertyFields = null;

     // constants
     static public final byte UNSPECIFIED_TYPE = 0;
     static public final byte TOP = 1;
     static public final byte BASE = 2;
     static public final byte FAULT = 3;
     static public final byte UNCONFORMITY = 4;

     // instance fields
     protected StsMarker marker = null;
     protected StsWell well = null;
     protected StsPoint location = null;

//     protected float mdepth = nullValue;
     protected int domain = 0;
 //    protected int domain = HDFWellMarker.DEPTH_DOMAIN;
 //    protected boolean isVisible = true;

     static int colorIndex = 0;
     static final float nullValue = StsParameters.nullValue;
     /** DB constructor */
     public StsWellMarker()
     {
     }

     public StsWellMarker(boolean persistent)
     {
         super(persistent);
     }

     /** constructor */
     private StsWellMarker(String name, StsWell well, byte type, StsPoint location) throws StsException
     {
         if (well==null)
         {
             throw new StsException(StsException.WARNING, "StsWellMarker.StsWellMarker:"
                       + " Cannot create a marker for a null well.");
         }
         if (name==null)
         {
             throw new StsException(StsException.WARNING, "StsWellMarker.StsWellMarker:"
                       + " Marker name is null.");
         }
         if (!StsMarker.isTypeValid(type))
         {
             throw new StsException(StsException.WARNING, "StsWellMarker.StsWellMarker:"
                       + " Marker type is invalid.");
         }
         marker = (StsMarker)currentModel.getObjectWithName(StsMarker.class, name);
         if (marker == null)  // create a new marker
         {
             marker = new StsMarker(name, type);
          }
         else if (type != marker.getType())  // same name but different type
         {
             // create a new marker with a unique name
 //             marker = new StsMarker(name, type);
             type = marker.getType();  // change wellMarker to same type
         }

         marker.addWellMarker(this);
         well.addMarker(this);
         this.well = well;
         this.location = location;
     }

     static public StsWellMarker constructor(String name, StsWell well, byte type, StsPoint location)
     {
         try
         {
             return new StsWellMarker(name, well, type, location);
         }
         catch(Exception e)
         {
             StsMessageFiles.errorMessage(e.getMessage());
             return null;
         }
     }

     public StsWellMarker(StsWell well, StsMarker marker, float mdepth) throws StsException
     {
         if (well==null)
         {
             throw new StsException(StsException.WARNING, "StsWellMarker.StsWellMarker:"
                       + " Cannot create a marker for a null well.");
         }
         if (marker==null)
         {
             throw new StsException(StsException.WARNING, "StsWellMarker.StsWellMarker:"
                       + " Marker is null.");
         }

         this.well = well;
         this.marker = marker;
         location = well.getPointAtMDepth(mdepth, false);

         marker.addWellMarker(this);
         well.addMarker(this);
     }

     static public StsWellMarker constructor(StsWell well, StsMarker marker, float mdepth)
     {
         try
         {
             return new StsWellMarker(well, marker, mdepth);
         }
         catch(Exception e)
         {
             StsMessageFiles.errorMessage(e.getMessage());
             return null;
         }
     }

     /** constructor for surface/well intersection */
     private StsWellMarker(StsModelSurface surface, StsWell well) throws StsException
     {
         super(false);
         this.well = well;
         float mdepth = computeWellSurfaceIntersection(surface, well);
         location.setM(mdepth);
         if(mdepth == nullValue) throw new StsException("Failed to construct well-surface intersection marker.");

         String name = surface.getName();
         marker = (StsMarker)currentModel.getObjectWithName(StsMarker.class, name);
         if (marker == null)  // create a new marker
         {
             marker = new StsMarker(false);
             marker.setName(name);
             marker.setType(StsMarker.SURFACE);
             marker.setStsColor(surface.getStsColor());
             marker.setModelSurface(surface);
             surface.setMarker(marker);
             currentModel.add(marker);
        }

         currentModel.add(this);
         well.addMarker(this); // can't add until marker is persistent in line above
         marker.addWellMarker(this);
     }

     public static StsWellMarker constructor(StsModelSurface surface, StsWell well)
     {
         if (well == null) return null;
         if (surface==null) return null;
         try
         {
              StsWellMarker marker = new StsWellMarker(surface, well);
              return marker;
         }
         catch(Exception e)
         {
             StsMessageFiles.errorMessage(e.getMessage());
             return null;
         }
     }

     static public float computeWellSurfaceIntersection(StsModelSurface surface, StsWell well)
     {
         if(surface == null) return nullValue;
         StsGridPoint gridPoint = well.computeGridIntersect(surface);
         if(gridPoint == null) return nullValue;
         if(isDepth)
             return well.getMDepthFromDepth(gridPoint.getZorT());
         else
             return well.getMDepthFromTime(gridPoint.getZorT());
     }

     public void resetLocation()
     {
         float mdepth = location.getM();
         location = well.getPointAtMDepth(mdepth, false);
     }

     public boolean initialize(StsModel model)
     {
//         resetLocation();
         StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
         if (velocityModel != null)  adjustTime();
         return true;
     }

     // Accessors
     public String getName() {  return marker.getName(); }
     public String toString() { return marker.getName(); }
     public StsWell getWell(){ return well; }
     public void setLocation(StsPoint xyz) { location = xyz; }
     public StsPoint getLocation() { return location; }
     public void setZ(float z) { location.setZ(z); }
     public float getZ() { return location.getZ(); }
     public float getMDepth() { return location.getM(); }
 //    public void setMDepth(float value) { mdepth = value; }
     public boolean setDomain(int domain)
     {
 /*
         switch (domain)
         {
             case HDFWellMarker.DEPTH_DOMAIN:
             case HDFWellMarker.TIME_DOMAIN:
                 this.domain = domain;
                 return true;
         }
         return false;
 */
         return true;
     }
     public int getDomain() { return domain; }
     public void setMarker(StsMarker marker) throws StsException
     {
         if (marker==null) throw new StsException(StsException.WARNING,
                 "StsWellMarker.setMarker:  marker cannot be set to null.");
         this.marker = marker;
     }
     public StsMarker getMarker() { return marker; }

     /** convenience methods */
     public void setStsColor(StsColor color) { marker.setStsColor(color); }
     public StsColor getStsColor() { return marker.getStsColor(); }

     public void setModelSurface(StsModelSurface s) { marker.setModelSurface(s); }
     public StsModelSurface getModelSurface() { return marker.getModelSurface(); }

     /** get well marker z value to use in ordering markers */
     public float getOrderingValue() { return this.getZ(); }

     public void adjustTime()
     {
         float t = well.getTimeFromMDepth(location.getM());
         location.setT(t);
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
         float[] xyz = location.getXYZorT();
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
         float[] xyz = location.getXYZorT();
         StsGLDraw.drawPoint(xyz, color, glPanel3d, width, height, viewshift);
         if(!displayName) return;
         GL gl = glPanel3d.getGL();
         glPanel3d.setViewShift(gl, viewshift + 2.0);
         gl.glDisable(GL.GL_LIGHTING);
         StsGLDraw.fontHelvetica12(glPanel3d.getGL(), xyz, marker.getName());
         gl.glEnable(GL.GL_LIGHTING);
         glPanel3d.resetViewShift(gl);
     }

     public StsFieldBean[] getDisplayFields()
     {
         if(displayFields == null)
         {
             displayFields = new StsFieldBean[]
             {
                 new StsBooleanFieldBean(StsWellMarker.class, "isVisible", "Enable"),
                 new StsColorComboBoxFieldBean(StsWellMarker.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
             };
         }
         return displayFields;
     }

     public StsFieldBean[] getPropertyFields()
     {
         if(propertyFields == null)
         {
             propertyFields = new StsFieldBean[]
             {
                 new StsFloatFieldBean(StsWellMarker.class, "z", "Z"),
                 new StsFloatFieldBean(StsWellMarker.class, "mDepth", "Measured")
             };
         }
         return propertyFields;
     }

     public XMLobject getXMLobject()
     {
         return new XMLobject();
     }

     public class XMLobject
     {
         String markerName;
         float depth;

         public XMLobject()
         {
             markerName = marker.name;
             depth = location.getM();
         }
     }
 }














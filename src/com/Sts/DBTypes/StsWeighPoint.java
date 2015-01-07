
/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie Class instantiated by the movie wizard.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.Iterator;

public class StsWeighPoint extends StsMainObject implements StsTreeObjectI
{
    private float azimuth = 0.0f;
    private float elevation = 0.0f;
    private float distance = 0.0f;
    private float zscale = 1.0f;
    private long time = 0l;
    private StsPoint center = null;
    private StsPoint cursor = null;
    private byte domain = StsProject.TD_TIME;
    private StsObjectRefList notes = null;
    public float[][] axisRanges = new float[2][0];
    public int direction = StsCursor3d.ZDIR;

    transient StsWeighPointDialog wpd = null;
    transient private StsModel model = null;
    transient private StsView view = null;
    static protected StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    /**
     * Default constructor
     */
    public StsWeighPoint()
    {

    }

    /**
     * constructor
     */
    public StsWeighPoint(StsModel model, StsWindowFamily family)
    {
        super(false);
        this.model = model;
        StsWeighPoint[] wps = (StsWeighPoint[])model.getCastObjectList(StsWeighPoint.class);
        setName("WayPoint" + wps.length);
        Iterator<StsView> windowViewIterator = family.getWindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = windowViewIterator.next();
            StsCursor3d cursor = view.glPanel3d.getCursor3d();
            direction = view.glPanel3d.getCursor3d().getCurrentDirNo();
            setCursor(cursor.getCurrentDirCoordinate(StsCursor3d.XDIR),
                      cursor.getCurrentDirCoordinate(StsCursor3d.YDIR),
                      cursor.getCurrentDirCoordinate(StsCursor3d.ZDIR));
            if(view instanceof StsView3d)
            {
                float[] parms = ((StsView3d)view).getCenterAndViewParameters();
                center = new StsPoint(parms[0], parms[1], parms[2]);
                azimuth = parms[4];
                distance = parms[3];
                elevation = parms[5];
                zscale = parms[6];
            }
            else if(view instanceof StsView2d)
            {
                float[][] viewRanges = ((StsView2d)view).getAxisRanges();
                axisRanges[0] = new float[] {viewRanges[0][0], viewRanges[0][1]};
                axisRanges[1] = new float[] {viewRanges[1][0], viewRanges[1][1]};
            }
		}
        domain = model.getProject().getZDomain();
        time = model.getProject().getProjectTime();

        addToModel();
//		refreshObjectPanel();
//        model.getActionManager().fireChangeEvent();
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            this.model = model;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWeighPoint.classInitialize() failed.", e, StsException.WARNING);
            return false;
        }
    }

    /**
     * Does this weighpoint have notes
     */
    public String getHasNotes()
    {
        if(notes == null)
            return new String("No");
        else
            return new String("Yes");
    }

    public boolean delete()
    {
        if(notes != null)
            notes.deleteAll();
        return super.delete();    
    }

    /**
     * Set the azimuth off North the user is facing
     * @param azimuth is degrees off North
     */
    public void setAzimuth(float azimuth) 
    { 
    	this.azimuth = azimuth; 
    	treeObjectSelected(); 
    	dbFieldChanged("azimuth", azimuth);
    	}
    
    public float getAzimuth() { return azimuth; }

    /**
     * Set the elevation the user is above horizontal
     * @param elevation is degrees off horizontal
     */
    public void setElevation(float elevation) { this.elevation = elevation; treeObjectSelected(); dbFieldChanged("elevation", elevation);}
    public float getElevation() { return elevation; }

    /**
     * Set the distance the user is from the COG
     * @param distance is units from Center of Gravity
     */
    public void setDistance(float distance) { this.distance = distance; treeObjectSelected(); dbFieldChanged("distance", distance);}
    public float getDistance() { return distance; }

    /**
     * Set the z scale
     * @param scale multiple of z
     */
    public void setZscale(float scale) { this.zscale = scale; treeObjectSelected(); dbFieldChanged("zscale", zscale);}
    public float getZscale() { return zscale; }

    /**
     * Set the center of the scene
     * @param xpos is the x position
     * @param ypos is the y position
     * @param zpos is the z position
     */
    public void setCenter(float xpos, float ypos, float zpos)
    {
        center.setX(xpos);
        center.setY(ypos);
        center.setZ(zpos);
    }
    public float[] getCenter() { return center.getPointValues();}
    public void setXCenter(float xpos) {  center.setX(xpos); treeObjectSelected(); dbFieldChanged("center", center);}
    public void setYCenter(float ypos) {  center.setY(ypos); treeObjectSelected(); dbFieldChanged("center", center);}
    public void setZCenter(float zpos) {  center.setZ(zpos); treeObjectSelected(); dbFieldChanged("center", center);}
    public float getXCenter() {  return center.getX(); }
    public float getYCenter() {  return center.getY(); }
    public float getZCenter() {  return center.getZ(); }

    /**
     * Set the cursor position
     * @param xpos is the x position
     * @param ypos is the y position
     * @param zpos is the z position
     */
    public void setCursor(float xpos, float ypos, float zpos)
    {
        cursor = new StsPoint(xpos, ypos, zpos);
        cursor.setX(xpos);
        cursor.setY(ypos);
        cursor.setZ(zpos);
    }

    public float[] getCursor() { return cursor.getPointValues();}
    public void setXCursor(float xpos) {  cursor.setX(xpos); treeObjectSelected(); dbFieldChanged("cursor", cursor);}
    public void setYCursor(float ypos) {  cursor.setY(ypos); treeObjectSelected(); dbFieldChanged("cursor", cursor);}
    public void setZCursor(float zpos) {  cursor.setZ(zpos); treeObjectSelected(); dbFieldChanged("cursor", cursor);}
    public float getXCursor() {  return cursor.getX(); }
    public float getYCursor() {  return cursor.getY(); }
    public float getZCursor() {  return cursor.getZ(); }

    /**
     * Display method for a weighPoint
     * @param glPanel3d the graphics context
     */
    public void display(StsGLPanel3d glPanel3d)
    {
        if (glPanel3d == null)
            return;
        /* Debating about the value of showing the weighPoints */
    }

    /**
     * Add a note to this weighPoint
     * @param jrnl - The text based note
     */
    public void addNote(String jrnl)
    {
        if(notes == null)
		  {
			  notes = StsObjectRefList.constructor(5, 5, "notes", this);
			  this.dbFieldChanged("notes", notes);
		  }
        StsNote note = new StsNote(jrnl);
//        note.addToModel();
        notes.add(note);
    }

    /**
     * Get all the notes associated with this WeighPoint
     */
    public StsObjectRefList getNotes()
    {
        return notes;
    }

    /**
     * WeighPoint selected on the Object tree
     */
     public void treeObjectSelected()
     {
        if(currentModel.win3d == null)
            return;
        //view = (StsView3d)currentModel.win3d.getView(StsView3d.class);
        getWeighPointClass().currentObject = this;
        currentModel.getProject().setZDomain(domain);

         changeModelView();

        if(getWeighPointClass().getEnableTimeJumps())
             currentModel.getProject().setProjectTime(time);

        currentModel.win3dDisplayAll();

                  // Open the note dialog
        if(getWeighPointClass().getDisplayWeighPointDialog())
        {
            if(wpd != null)
            {
                wpd.hide();
                wpd = null;
            }
            wpd = new StsWeighPointDialog(model, this, false);
            wpd.setVisible(true);
        }
    }
    public boolean changeModelView()
    {
        float[] parms = new float[] {center.getX(), center.getY(), center.getZ(), distance, azimuth, elevation, zscale};

        Iterator<StsView> windowViewIterator = currentModel.getMainWindowFamilyViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = windowViewIterator.next();

            // ToDo: Does not work when sliders are actual coordinates.
            view.glPanel3d.getCursor3d().setCoordinates(cursor);
            view.glPanel3d.getCursor3d().setCurrentDirNo(direction);           
            if(view instanceof StsView3d)
                view.changeModelView3d(parms);
            else if(view instanceof StsView2d)
            {
                if(axisRanges != null)
                    view.changeModelView2d(axisRanges);
                //System.out.println("Not enabled yet for 2D....");
            }
		}
        return true;
    }

    static public StsWeighPointClass getWeighPointClass()
    {
       return (StsWeighPointClass)currentModel.getCreateStsClass(StsWeighPoint.class);
    }

         public boolean anyDependencies()
                             {
                                 return false;
                             }
         public StsFieldBean[] getDisplayFields()
         {
                           if(displayFields == null)
                                 {
                                     displayFields = new StsFieldBean[]
                                     {
                                         new StsFloatFieldBean(StsWeighPoint.class, "xCenter", true, "Center X:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "yCenter", true, "Center Y:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "zCenter", true, "Center Z:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "azimuth", true, "Azimuth:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "elevation", true, "Elevation:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "distance", true, "Direction:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "zscale", true, "Z Scale:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "xCursor", true, "Cursor X:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "yCursor", true, "Cursor Y:"),
                                         new StsFloatFieldBean(StsWeighPoint.class, "zCursor", true, "Cursor Z:")
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
                                         new StsStringFieldBean(StsWeighPoint.class, "name", true, "Name:"),
                                         new StsStringFieldBean(StsWeighPoint.class, "hasNotes", false, "Has Notes:"),
                                         new StsStringFieldBean(StsWeighPoint.class, "zDomainString", false, "Z Domain:"),
                                     };
                                 }
                                 return propertyFields;
                             }

         public Object[] getChildren() { return new Object[0]; }

         public StsObjectPanel getObjectPanel()
                             {
                                 if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
                                 return objectPanel;
                             }

         public String getZDomainString()
                   {
                       if(domain == StsProject.TD_TIME) return "Time";
                       else return "Depth";
                   }
     }
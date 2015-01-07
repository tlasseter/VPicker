package com.Sts.Actions.Wizards.EditWellMarkers;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsAddWellMarker extends StsAction
{
    StsWellViewModel wellViewModel;
    StsFieldBeanPanel addMarkerPanel = new StsFieldBeanPanel();
    StsGroupBox selectGroupBox = new StsGroupBox("1. Select/Create");
    StsGroupBox mdepthGroupBox = new StsGroupBox("2. Pick/Edit");
    StsFieldBeanPanel selectOperationsPanel;
    StsFieldBeanPanel mdepthPanel;
    StsComboBoxFieldBean markerTypeBean;
    StsButtonListFieldBean oldOrNewRadioButtonBean;
    String newOrOld = oldString;
    StsComboBoxFieldBean markersComboBoxBean;
    StsStringFieldBean newMarkerStringBean;
    StsStringFieldBean surfaceStringBean;
    StsComboBoxFieldBean surfacesComboBoxBean;
    StsFloatFieldBean mdepthFloatBean;
    StsFloatFieldBean lengthFloatBean;
    StsFloatFieldBean offsetFloatBean;
    StsDateFieldBean dateBean;
    StsButton acceptBtn;

    StsComboBoxFieldBean subTypeComboBoxBean;
    byte markerSubType = StsEquipmentMarker.SENSOR;
    StsMarker[] markers = null;
    StsMarker marker = null;
    StsWellMarker wellMarker = null;
    StsSurface surface = null;
    String markerSurfaceName = null;
    float mdepth;
    float length = 10.0f;
    float zOffset = 0.0f;
    long perfTime = 0l;
    byte subType = StsEquipmentMarker.SENSOR;
    byte markerType = StsMarker.GENERAL;
    String newMarkerName = "";
    GridBagConstraints gbcMarkerBean;

    static String newString = "new";
    static String oldString = "existing";
    static String[] newOrOldStrings = new String[] { oldString, newString };
    static String noSurfaceString = "none";

    static StsMarker nullMarker = null;

    public StsAddWellMarker(StsModel model, StsWellViewModel wellViewModel)
    {
        super(model);
        this.wellViewModel = wellViewModel;
    }

    public boolean start()
    {
        return initializeAddMarkerPanel();
    }

    /** Permanently at the top of the addMarkerPanel, a comboBox of markerTypes is displayed.
     *  Selecting type controls what and how will be displayed below
     */

    private boolean initializeAddMarkerPanel()
    {
        try
        {
        	addMarkerPanel.gbc.gridwidth = 2;
            addMarkerPanel.addEndRow(selectGroupBox);
            addMarkerPanel.addEndRow(mdepthGroupBox);

            markerTypeBean = new StsComboBoxFieldBean(this, "markerType", "Marker Type:", StsMarker.markerTypes);
            selectGroupBox.add(markerTypeBean);

            buildSelectOperationsPanel();
            buildDepthPanel();

            StsButton acceptBtn = new StsButton("Accept", "Press to save the current edited or new marker.", this, "acceptMarker");

            addMarkerPanel.gbc.gridwidth = 1;
            addMarkerPanel.addEndRow(acceptBtn);
            addMarkerPanel.add(acceptBtn);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsAddWellMarker.initializeAddMarkerPanel() failed.",
                                         e, StsException.WARNING);
            return false;
        }
   }

    /** If we have existing markers, get a list of those not already on this well.
     *  If there are any, provide two radio buttons (Existing and New) which allow
     *  the user to create a new wellMarker from an existing marker or a completely
     *  new marker.
     *
     *  If "Existing" radioButton is selected, a comboBox of existing markers is displayed below.
     *  If "New" radioButton is selected, an empty textBox is displayed for user to fill in.
     */
    private void buildSelectOperationsPanel()
    {
        if(selectOperationsPanel != null)
            selectGroupBox.remove(selectOperationsPanel);

        selectOperationsPanel = new StsFieldBeanPanel();
        selectGroupBox.gbc.gridwidth = 2;
        selectGroupBox.add(selectOperationsPanel);

        newMarkerStringBean = new StsStringFieldBean(this, "newMarkerName", "New:");
        newMarkerStringBean.setSize(30, 15);

        markers = getMarkerList(wellViewModel.getWell());
        if(markers != null)
        {
            oldOrNewRadioButtonBean = new StsButtonListFieldBean(this, "newOrOldMarker", "Select:", newOrOldStrings, true);
            oldOrNewRadioButtonBean.setValueObject(newOrOld);
            selectOperationsPanel.add(oldOrNewRadioButtonBean);

            markersComboBoxBean = new StsComboBoxFieldBean(this, "marker", "Existing:", markers);
            if(newOrOld == oldString) // select marker from list
            {
                selectOperationsPanel.add(markersComboBoxBean);
                if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
                {
                    surfaceStringBean = new StsStringFieldBean(this, "markerSurfaceName", "Surface:");
                    selectOperationsPanel.add(surfaceStringBean);
                }
            }
            else // create new marker
            {
                selectOperationsPanel.add(newMarkerStringBean);
                if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
                {
                    StsObject[] surfaces = model.getObjectList(StsSurface.class);
                    if(surfaces.length > 0)
                    {
                        surfacesComboBoxBean = new StsComboBoxFieldBean(this, "surface", "Surface:", surfaces);
                    }
                    else
                    {
                        surfacesComboBoxBean = new StsComboBoxFieldBean(this, "noSurface", "Surface:", new String[]
                            { noSurfaceString } );
                        surface = null;
                    }
                    surfacesComboBoxBean.setSelectedIndex(0);
                    selectOperationsPanel.add(surfacesComboBoxBean);
                }
            }
            markersComboBoxBean.setSelectedIndex(0);
        }
        else // create new marker
        {
            selectOperationsPanel.add(newMarkerStringBean);
            if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
        	{
        		StsObject[] surfaces = model.getObjectList(StsSurface.class);
        		if(surfaces.length > 0)
        		{
        			surfacesComboBoxBean = new StsComboBoxFieldBean(this, "surface", "Surface:", surfaces);
        		}
        		else
        		{
        			surfacesComboBoxBean = new StsComboBoxFieldBean(this, "noSurface", "Surface:", new String[]
        			                                            { noSurfaceString } );
        			surface = null;
        		}
        		surfacesComboBoxBean.setSelectedIndex(0);
        		selectOperationsPanel.add(surfacesComboBoxBean);
        	}
        }
    }

    public String getMarkerSubType()
    {
        return StsEquipmentMarker.typeToString(markerSubType);
    }

    public void setMarkerSubType(String subType)
    {
    	byte newMarkerSubType = StsEquipmentMarker.stringToType(subType);
        if(newMarkerSubType == markerSubType) return;
        markerSubType = newMarkerSubType;
    }

    private void buildDepthPanel()
    {
        if(mdepthPanel != null)
             mdepthGroupBox.remove(mdepthPanel);

         mdepthPanel = new StsFieldBeanPanel();
         mdepthGroupBox.add(mdepthPanel);

         StsWell well = wellViewModel.getWell();
         float mdepthMin = well.getMinMDepth();
         float mdepthMax = well.getMaxMDepth();
         mdepthFloatBean = new StsFloatFieldBean(this, "mdepth", mdepthMin, mdepthMax, "Measure Depth:");
         mdepthPanel.add(mdepthFloatBean);
         if(markerType == StsMarker.PERFORATION)
         {
             lengthFloatBean = new StsFloatFieldBean(this, "length", 1, 1000, "Perforation Length:");
             mdepthPanel.add(lengthFloatBean);
             dateBean = new StsDateFieldBean(this, "perfTime", true, "Time:");
             dateBean.setValue(StsDateFieldBean.convertToString(System.currentTimeMillis()));
             mdepthPanel.add(dateBean);
         }
         else if(markerType == StsMarker.OFFSET_SURFACE)
         {
             offsetFloatBean = new StsFloatFieldBean(this, "zOffset", -500, 500, "Z Offset (- is Up)");
             mdepthPanel.add(offsetFloatBean);
         }
         else if(markerType == StsMarker.EQUIPMENT)
         {
        	 subTypeComboBoxBean = new StsComboBoxFieldBean(this, "markerSubType", "SubType:", StsEquipmentMarker.subTypeStrings);
        	 mdepthPanel.add(subTypeComboBoxBean);
         }
    }

    private void rebuildAddMarkerPanel()
    {
        buildSelectOperationsPanel();
        buildDepthPanel();
        wellViewModel.rebuild();
    }

    private StsMarker[] getMarkerList(StsWell well)
    {
        StsObjectRefList existingMarkersList = well.getMarkers();
        StsObjectList existingMarkers = null;
        if(existingMarkersList != null) existingMarkers = existingMarkersList.getList();

        StsClass markerClass = model.getStsClass(StsMarker.class);
        if(markerClass == null) return null;
        //StsObjectList markersList = markerClass.getObjectListOfTypeExcluding(markerType, existingMarkers);    // Need to allow multiple markers on a horizontal well
        StsObjectList markersList = null;
        if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
        {
            markersList = markerClass.getObjectListOfTypeExcluding(StsMarker.SURFACE, null);
            markersList.addList(markerClass.getObjectListOfTypeExcluding(StsMarker.OFFSET_SURFACE, null));
        }
        else
            markersList = markerClass.getObjectListOfTypeExcluding(markerType, null);

        if(markersList == null) return null;
        return (StsMarker[])markersList.getCastListCopy();
    }

    public void setNewOrOldMarker(String newOrOld)
    {
//         if(this.newOrOld == newOrOld) return;
        this.newOrOld = newOrOld;
        rebuildAddMarkerPanel();
    }

    public String getNewOrOldMarker()
    {
        return newOrOld;
    }

    public void setMarker(StsMarker marker)
    {
        this.marker = marker;
        wellMarker = null;

        //markerType = marker.getType();   This is not required since the list only contains markers of the selected type. It also causes problems with offset_surface markers since it does not allow the selection of a surface marker to create an offset surface marker.
        if(markerTypeBean != null) markerTypeBean.doSetValueObject(StsMarker.markerTypes[markerType]);

        if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
        {
            if(surfaceStringBean != null) 
                surfaceStringBean.setValue(marker.getSurface().getName());
        }
        else if(markerType == StsMarker.PERFORATION) // Nothing to fill since a perforation is unrelated between wells
        {
        	;
        }
        else if(markerType == StsMarker.EQUIPMENT) // Nothing to fill since equipment is unrelated between wells
        {
        	;
        }
    }

    public StsMarker getMarker()
    {
        return marker;
    }

    public void setNewMarkerName(String name)
    {
        newMarkerName = name;
        marker = null;
        wellMarker = null;
    }

    public String getNewMarkerName()
    {
        return newMarkerName;
    }

    public void setMarkerType(String stringType)
    {
        byte newMarkerType = StsMarker.stringToType(stringType);
        if((model.getNObjects(StsModelSurface.class) == 0) && (model.getNObjects(StsSurface.class) == 0) && ((newMarkerType == StsMarker.SURFACE) || (newMarkerType == StsMarker.OFFSET_SURFACE)))
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "Unable to create surface markers. Project has no surfaces.");
            markerTypeBean.setSelectedItem(StsMarker.typeToString(markerType));
            return;
        }
        if(newMarkerType == markerType) return;
        markerType = newMarkerType;
        rebuildAddMarkerPanel();
    }

    public String getMarkerType()
    {
        return StsMarker.typeToString(markerType);
    }

    public void setSurface(StsObject surface)
    {
        this.surface = (StsSurface)surface;
    }
    public StsObject getSurface() { return surface; }

    public void setMarkerSurfaceName(String name)
     {
         markerSurfaceName = name;
     }
     public String getMarkerSurfaceName()
     {
         return markerSurfaceName;
     }

    public void setNoSurface(String surface) { this.surface = null; }
    public String getNoSurface() { return null; }

    public void setMdepth(float mdepth)
    {
        this.mdepth = mdepth;
    }
    public float getMdepth() { return mdepth; }

    public void setZOffset(float zoffset)
    {
        this.zOffset = zoffset;
    }
    public float getZOffset() { return zOffset; }

    public void setLength(float len)
    {
        this.length = len;
    }
    public float getLength() { return length; }

    public String getPerfTime()
    {
        return StsDateFieldBean.convertToString(perfTime);
    }
    public void setPerfTime(String time)
    {
         this.perfTime = StsDateFieldBean.convertToLong(time);
    }
    public StsJPanel getAddMarkerPanel()
    {
        return addMarkerPanel;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        try
        {
            if (mouse.getCurrentButton() == StsMouse.LEFT)
            {
                {

                    double depth = wellViewModel.getMdepthFromMouseY(mouse.getMousePoint().y, glPanel);
                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
                    if (buttonState == StsMouse.PRESSED || buttonState == StsMouse.DRAGGED)
                    {
                        wellViewModel.cursorPicked = true;
                        wellViewModel.moveCursor(mouse.getMousePoint().y);
                        wellViewModel.repaint();
                        this.mdepth = (float)depth;
                        mdepthFloatBean.setValue(mdepth);
                    }
                    else if (buttonState == StsMouse.RELEASED)
                    {
                        wellViewModel.cursorPicked = false;
                        this.mdepth = (float)depth;
                        mdepthFloatBean.setValue(mdepth);
                    }

                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "performMouseAction", e);
            return false;
        }
    }

    public void acceptMarker()
    {
    	createChangeWellMarker();
    }

    private void createChangeWellMarker()
    {
        StsWell well = wellViewModel.getWell();
        if(newOrOld == newString)
        {
            String newMarkerName = getNewMarkerName();
            // Horizontal wells may intersect the same bed multiple times.
            /*
            if(well.getMarker(getNewMarkerName()) != null)
            {
                new StsMessage(addMarkerPanel, StsMessage.ERROR,
                        "Failed to create marker (" + newMarkerName + "). it already exists for this well (" + well.getName() + ").");
                return;
            }
            */
            if(getNewMarkerName().length() == 0)
            {
                new StsMessage(addMarkerPanel, StsMessage.ERROR, "Invalid marker name, specify valid name first.");
                return;
            }
        }
        if(wellMarker == null)
        {
            if(marker == null)
            {
               StsMarker newMarker;
               if(markerType == StsMarker.SURFACE)
               {
                   newMarker = new StsMarker(newMarkerName, markerType, surface);
                   wellMarker = StsWellMarker.constructor(well, newMarker, mdepth);
                   StsPoint location = well.getPointAtMDepth(mdepth, false);
                   wellMarker.setLocation(location);
               }
               else if (markerType == StsMarker.OFFSET_SURFACE)
               {
                   newMarker = new StsMarker(newMarkerName, markerType, surface);
                   wellMarker = StsOffsetSurfaceMarker.constructor(well, newMarker, mdepth, zOffset);
                   StsPoint location = well.getPointAtMDepth(mdepth, false);
                   wellMarker.setLocation(location);
               }
               else if(markerType == StsMarker.PERFORATION)
               {
            	   wellMarker = StsPerforationMarker.constructor(newMarkerName, well, StsMarker.PERFORATION, well.getPointAtMDepth(mdepth, false), length, 1, System.currentTimeMillis());
               }
               else if(markerType == StsMarker.EQUIPMENT)
               {
            	   wellMarker = StsEquipmentMarker.constructor(newMarkerName, well, mdepth, subType);
               }
               else
               {
                   newMarker = new StsMarker(newMarkerName, markerType);
                   wellMarker = StsWellMarker.constructor(well, newMarker, mdepth);
               }
            }
            else
            {
               if(markerType == StsMarker.SURFACE)
               {
                   wellMarker = StsWellMarker.constructor(well, marker, mdepth);
               }
               else if (markerType == StsMarker.OFFSET_SURFACE)
               {
            	   wellMarker = StsOffsetSurfaceMarker.constructor(well, marker, mdepth, zOffset);
               }
               StsPoint location = well.getPointAtMDepth(mdepth, false);
               wellMarker.setLocation(location);
            }
        }
        markerTypeBean.setSelectedItem(StsMarker.typeToString(StsMarker.GENERAL));
        setMdepth(0.0f);
        setNewMarkerName("");
        rebuildAddMarkerPanel();
        wellViewModel.repaint();
        model.win3dDisplayAll();
    }
}

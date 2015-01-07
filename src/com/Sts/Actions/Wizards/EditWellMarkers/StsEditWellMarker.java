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

public class StsEditWellMarker extends StsAction
{
    StsWellViewModel wellViewModel;
    StsJPanel editMarkerPanel = StsJPanel.addInsets();
    StsGroupBox selectGroupBox = new StsGroupBox("Step 1 - Select & Edit Marker");
    StsGroupBox mdepthGroupBox = new StsGroupBox("Step 2 - Move/Delete Marker");
    StsFieldBeanPanel selectOperationsPanel;
    StsFieldBeanPanel mdepthPanel;
    StsComboBoxFieldBean markerTypeBean;
    StsComboBoxFieldBean markersComboBoxBean;
    StsComboBoxFieldBean surfacesComboBoxBean;
    StsFloatFieldBean mdepthFloatBean;
    StsFloatFieldBean lengthFloatBean;
    StsIntFieldBean numShotsBean;
    StsComboBoxFieldBean subTypeComboBoxBean;
    StsDateFieldBean bornField;

    StsButton acceptBtn;

    StsWell well = null;
    StsWellMarker[] wellMarkers = null;
    StsMarker marker = null;
    StsWellMarker wellMarker = null;
    StsSurface surface = null;
    String markerSurfaceName = null;
    float mdepth, length = 10.0f;
    int nShots = 1;
    long time = -1l;
    byte markerSubType = StsEquipmentMarker.SENSOR;
    byte markerType = StsMarker.GENERAL;
    String newMarkerName = "";
    GridBagConstraints gbcMarkerBean;

    static String newString = "New";
    static String oldString = "Existing";
    static String[] newOrOldStrings = new String[] { oldString, newString };
    static String noSurfaceString = "None";

    static StsMarker nullMarker = null;

    public StsEditWellMarker(StsModel model, StsWellViewModel wellViewModel)
    {
        super(model);
        this.wellViewModel = wellViewModel;
        well = wellViewModel.getWell();
    }

    public boolean start()
    {
        return initializeEditMarkerPanel();
    }

    /** Permanently at the top of the addMarkerPanel, a comboBox of markerTypes is displayed.
     *  Selecting type controls what and how will be displayed below
     */

    private boolean initializeEditMarkerPanel()
    {
        try
        {
        	editMarkerPanel.gbc.gridwidth = 2;
            editMarkerPanel.addEndRow(selectGroupBox);
            editMarkerPanel.addEndRow(mdepthGroupBox);

            wellMarkers = getMarkerList();
            markersComboBoxBean = new StsComboBoxFieldBean(this, "wellMarker", "Well Markers:", wellMarkers);
            selectGroupBox.add(markersComboBoxBean);
            // the following will classInitialize the selectOperationsPanel which will be added to the selectGroupBox
            markersComboBoxBean.setSelectedIndex(0);
 //           buildSelectOperationsPanel();
            buildDepthPanel();

            StsButton acceptBtn = new StsButton("Accept", "Press to save the current edited or new marker.", this, "acceptMarker");
            StsButton deleteButton = new StsButton("Delete", "Click to delete current marker.", this, "deleteWellMarker");

            editMarkerPanel.gbc.gridwidth = 1;
            editMarkerPanel.addToRow(deleteButton);
            editMarkerPanel.addEndRow(acceptBtn);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEditWellMarker.initializeEditMarkerPanel() failed.",
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

        markerTypeBean = new StsComboBoxFieldBean(this, "markerType", "Marker Type:", StsMarker.markerTypes);
        markerTypeBean.setSelectedItem(StsMarker.typeToString(markerType));
        markerTypeBean.setEditable(false);
        selectOperationsPanel.add(markerTypeBean);

        if((markerType == StsMarker.SURFACE) || (markerType == StsMarker.OFFSET_SURFACE))
        {
            getConstructSurfacesComboBox();
            selectOperationsPanel.add(surfacesComboBoxBean);
        }
    }

    private void getConstructSurfacesComboBox()
    {
        if(surfacesComboBoxBean == null)
        {
            StsObject[] surfaces = model.getObjectList(StsSurface.class);
            int nSurfaces = surfaces.length;
            Object[] surfaceItems = new Object[nSurfaces + 1];
            surfaceItems[0] = noSurfaceString;
            for (int n = 0; n < nSurfaces; n++)
                surfaceItems[n + 1] = surfaces[n];
            surfacesComboBoxBean = new StsComboBoxFieldBean(this, "surface", "Surface:", surfaceItems);
       }
        surface = marker.getSurface();
        if(surface != null)
            surfacesComboBoxBean.setSelectedItem(surface);
        else
            surfacesComboBoxBean.setSelectedItem(noSurfaceString);
    }

    private void buildDepthPanel()
    {
        if(mdepthPanel != null)
             mdepthGroupBox.remove(mdepthPanel);

         mdepthPanel = new StsFieldBeanPanel();
         mdepthGroupBox.add(mdepthPanel);

         float mdepthMin = well.getMinMDepth();
         float mdepthMax = well.getMaxMDepth();
         mdepthFloatBean = new StsFloatFieldBean(this, "mdepth", mdepthMin, mdepthMax, "Measured Depth:");
         mdepthPanel.add(mdepthFloatBean);
         if(markerType == StsMarker.PERFORATION)
         {
             lengthFloatBean = new StsFloatFieldBean(this, "length", 1, 1000, "Perforation Length:");
             numShotsBean = new StsIntFieldBean(this, "numShots", 1, 100, "Number of Shots:");
             bornField = new StsDateFieldBean(this, "timeString", true, "Perforation Date:");
             bornField.setFormat(model.getProject().getTimeDateFormat());
             mdepthPanel.add(lengthFloatBean);
             mdepthPanel.add(numShotsBean);
             mdepthPanel.add(bornField);
         }
         else if(markerType == StsMarker.EQUIPMENT)
         {
        	 subTypeComboBoxBean = new StsComboBoxFieldBean(this, "markerSubType", "SubType:", StsEquipmentMarker.subTypeStrings);
        	 mdepthPanel.add(subTypeComboBoxBean);
         }
     }

    public void deleteWellMarker()
    {
        well.deleteMarker(wellMarker);
        wellViewModel.repaint();
        markersComboBoxBean.getComboBox().removeItem(wellMarker);
    }

    private void rebuildAddMarkerPanel()
    {
        buildSelectOperationsPanel();
        buildDepthPanel();
        wellViewModel.rebuild();
    }

    private StsWellMarker[] getMarkerList()
    {
        StsObjectRefList wellMarkers = well.getMarkers();
        return (StsWellMarker[])wellMarkers.getCastList();
    }

    public void setWellMarker(StsWellMarker wellMarker)
    {
        this.wellMarker = wellMarker;
        marker = wellMarker.getMarker();
        mdepth = wellMarker.getMDepth();
		if(wellMarker instanceof StsPerforationMarker)
		{
			markerType = StsMarker.PERFORATION;
	        length = ((StsPerforationMarker)wellMarker).getLength();
            nShots = ((StsPerforationMarker)wellMarker).getNumShots();
            time = ((StsPerforationMarker)wellMarker).getBornDateLong();
		}
		else if(wellMarker instanceof StsEquipmentMarker)
		{
			markerType = StsMarker.EQUIPMENT;
	        markerSubType = ((StsEquipmentMarker)wellMarker).getSubType();
		}
        rebuildAddMarkerPanel();

    }

    public StsWellMarker getWellMarker()
    {
        return wellMarker;
    }

    public void setMarkerType(String stringType)
    {
        byte newMarkerType = StsMarker.stringToType(stringType);
        if(newMarkerType == markerType) return;
        markerType = newMarkerType;
        marker.fieldChanged("type", markerType);
        rebuildAddMarkerPanel();
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

    public String getMarkerType()
    {
        return StsMarker.typeToString(markerType);
    }

    public void setSurface(StsObject surface)
    {
        this.surface = (StsSurface)surface;
        marker.fieldChanged("surface", this.surface);
		this.surface.fieldChanged("marker", marker);
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

    public StsJPanel getEditMarkerPanel()
    {
        return editMarkerPanel;
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
                        this.mdepth = (float)depth;
                        mdepthFloatBean.setValue(mdepth);
                        wellViewModel.repaint();
                    }
                    else if (buttonState == StsMouse.RELEASED)
                    {
                        wellViewModel.cursorPicked = false;
                        setMdepth((float)depth);
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

    public void setTimeString(String timeString)
    {
        time = model.getProject().getTime(timeString);
    }
    public String getTimeString()
    {
        return model.getProject().getTimeString(time);
    }

    public void setNumShots(int nShots)
    {
        this.nShots = nShots;
    }
    public int getNumShots() { return nShots; }

    public void setLength(float len)
    {
        this.length = len;
    }
    public float getLength() { return length; }

    public void acceptMarker()
    {
        StsPoint location = well.getPointAtMDepth(mdepth, false);
        wellMarker.fieldChanged("location", location);
        if(wellMarker instanceof StsPerforationMarker)
        {
        	((StsPerforationMarker)wellMarker).setLength(length);
            ((StsPerforationMarker)wellMarker).setNumShots(nShots);
            ((StsPerforationMarker)wellMarker).setBornDate(time);
        }
        else if(wellMarker instanceof StsEquipmentMarker)
        {
        	((StsEquipmentMarker)wellMarker).setSubType(markerSubType);
        }
        wellViewModel.repaint();
        model.win3dDisplayAll();
    }
}


//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PerforationAttributes;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsCorrelateMarkers extends StsWizardStep
	implements ActionListener, ListSelectionListener
{
    StsHeaderPanel header = new StsHeaderPanel();
    StsObject[] sensors = null;
	JPanel mainPanel = new JPanel(new BorderLayout());
	StsDoubleListPanel panel = new StsDoubleListPanel();

    public StsCorrelateMarkers(StsWizard wizard)
    {
		super(wizard);
		setPanels(mainPanel, header);
        mainPanel.add(panel);
    	panel.setLeftListRenderer(new StsColorListRenderer());
    	panel.setRightListRenderer(new StsColorListRenderer());
        panel.setLeftTitle("Available Sensors:");
        panel.setRightTitle("Available Markers:");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Correlate Perforation Markers with Sensors");
        panel.getRightButton().setToolTipText("Remove Perforation Markers with Sensors");

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Perforation Correlation with Sensors");
        header.setSubtitle("Correlate with Well Markers");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PerforationAttributes");        
        header.setInfoText(wizardDialog,"(1) Select a Perforation Marker from left list.\n" +
                                  "(2) Select the related Sensor from the right list.\n" +
                                  "   ***** Perforation marker list is generated from from selected well *****\n" +
                                  "(3) Press the > button to assign the marker to the sensor or\n" +
                                  "    Press the < button to unassign the marker from the selected sensor.\n" +
                                  "(4) Press the Finish Button to dismiss wizard.");
    }
	public void setSensors(StsDynamicSensor[] sensors) { this.sensors = sensors; }

    public void initSensorList()
    {
        sensors = model.getObjectList(StsDynamicSensor.class);
    	if( sensors != null )
        {
            int nSensors = sensors == null ? 0 : sensors.length;
            if( nSensors > 0 )
            {
                StsColorListItem[] items = new StsColorListItem[nSensors];
                StsColor color = new StsColor(StsColor.BLACK);
                for( int i=0; i<nSensors; i++ )
                {
                    String sensorName = new String(sensors[i].getName());

                    StsPerforationMarker m = ((StsDynamicSensor)sensors[i]).getPerforationMarker();
                    if (m == null) sensorName += " ( )";
                    else sensorName += " (" + m.getName() + ")";

                    items[i] = new StsColorListItem(color, sensorName);
                }
                int[] selected = panel.getSelectedLeftIndices();

                // check to see if any of the currently selected have been removed
                if( selected != null )
                {
                	int nSelected = selected.length;
                	for( int i=0; i<selected.length; i++ )
                    	if( selected[i] >= nSensors ) nSelected--;
                    if( nSelected == 0 ) selected = null;
                    else
                    {
                    	int[] newSelected = new int[nSelected];
                        int iSelected = 0;
                        for( int i=0; i<nSelected; i++ )
	                    	if( selected[i] < nSensors )
                            {
                            	newSelected[iSelected] = selected[i];
                                iSelected++;
                            }
                        selected = new int[nSelected];
                        System.arraycopy(newSelected, 0, selected, 0, nSelected);
                    }

                }
                panel.setLeftItems(items);
                if( selected != null ) panel.setSelectedLeftIndices(selected);
            }
        }
    }

    private void initMarkerList()
    {
        StsWell well = ((StsPerforationAttributesWizard)wizard).getSelectedWell();
		StsPerforationMarker[] markers = well.getPerforationMarkers();
        if(markers == null) return;

        int nMarkers = markers.length;
        if(nMarkers == 0) return;

        StsColorListItem[] items = new StsColorListItem[nMarkers];
        for( int i=0; i<nMarkers; i++ )
        {
            if(markers[i].getName() == null)
            	continue;

            String markerName = new String(markers[i].getName());
            markerName = well.getName() + " - " + markerName;

            items[i] = new StsColorListItem(markers[i].getStsColor(), markerName);
        }

        items = (StsColorListItem[])StsMath.arraycopy(items, nMarkers);
        panel.setRightItems(items);
    }

    public boolean start()
    {
		refreshLists();
        // select all of the surfaces
        if(sensors != null)
        {
        	int[] selected = new int[sensors.length];
            for( int i=0; i<sensors.length; i++ ) selected[i] = i;
            panel.setSelectedLeftIndices(selected);
        }
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public void refreshLists()
    {
    	initSensorList();
    	initMarkerList();
        refreshButtons();
    }

    public String trimMarkerName(String fullString)
    {
    	int index = fullString.indexOf(" - ");
        if( index > 0 )
			return fullString.substring(index+3, fullString.length());
        return null;
    }

    public String trimSensorName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index-1);
        return null;
    }

    public StsPerforationMarker[] getSelectedMarkers()
    {
        StsPerforationMarker[] markers = null;

        try
        {
            Object[] items = panel.getSelectedRightItems();
            int nItems = items == null ? 0 : items.length;
            if( nItems == 0) return null;

            markers = new StsPerforationMarker[nItems];
            for( int i=0; i<nItems; i++ )
            {
                StsColorListItem item = (StsColorListItem)items[i];
                markers[i] = (StsPerforationMarker)model.getObjectWithName(StsPerforationMarker.class, trimMarkerName(item.getName()));
            }
    	    return markers;
        }
        catch(Exception e)
        {
            StsException.outputException("StsCorrelateMarkers.getSelectedMarkers() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	refreshButtons();
    }

    public void refreshButtons()
    {
		Object[] sensorItems = panel.getSelectedLeftItems();
		Object[] markerItems = panel.getSelectedRightItems();
    	int nSensorItems = sensorItems == null ? 0 : sensorItems.length;
    	int nMarkerItems = markerItems == null ? 0 : markerItems.length;
		panel.enableLeftButton(nSensorItems > 0 && nSensorItems == nMarkerItems);

		if(nSensorItems > 0 && nMarkerItems == 0)
        {
			for( int i=0; i<nSensorItems; i++ )
            {
            	StsColorListItem item = (StsColorListItem) sensorItems[i];
            	StsDynamicSensor sensor = null;
                try { sensor = (StsDynamicSensor)model.getObjectWithName(StsDynamicSensor.class, trimSensorName(item.getName())); }
                catch(Exception ex) { }
                if( sensor != null && sensor.getPerforationMarker() != null )
                {
					panel.enableRightButton(true);
                    return;
                }
            }
        }
		panel.enableRightButton(false);
    }
    public void actionPerformed(ActionEvent e)
    {
    	if( e.getActionCommand().equals("<") )
        {
			Object[] sensorItems = panel.getSelectedLeftItems();
			Object[] markerItems = panel.getSelectedRightItems();
	    	int nSensorItems = sensorItems == null ? 0 : sensorItems.length;
    		int nMarkerItems = markerItems == null ? 0 : markerItems.length;
            if( nSensorItems == nMarkerItems )
            {
                StsWell well = ((StsPerforationAttributesWizard)wizard).getSelectedWell();
                for( int i=0; i<nSensorItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem) sensorItems[i];
                    StsDynamicSensor sensor = null;
                    try { sensor = (StsDynamicSensor)model.getObjectWithName(StsDynamicSensor.class, trimSensorName(item.getName())); }
                    catch(Exception ex) { }

                    item = (StsColorListItem) markerItems[i];
                    StsPerforationMarker marker = (StsPerforationMarker)well.getMarker(trimMarkerName(item.getName()));
                    if( sensor != null && marker != null )
                    {
                    	sensor.setPerforationMarker(marker);
                        if(!sensor.computeDistanceToPerforation())
                            sensor.setPerforationMarker(null);
                    }
                }
                refreshLists();
            }
        }
        else if( e.getActionCommand().equals(">") )
        {
			Object[] sensorItems = panel.getSelectedLeftItems();
	    	int nSensorItems = sensorItems == null ? 0 : sensorItems.length;
            for( int i=0; i<nSensorItems; i++ )
            {
                StsColorListItem item = (StsColorListItem) sensorItems[i];
                StsDynamicSensor sensor = null;
                try { sensor = (StsDynamicSensor)model.getObjectWithName(StsDynamicSensor.class, trimSensorName(item.getName())); }
                catch(Exception ex) { }
                if( sensor != null )
                {
                	StsPerforationMarker oldMarker = sensor.getPerforationMarker();
                    if( oldMarker != null ) oldMarker.setModelSurface(null);
                    sensor.setPerforationMarker(null);
                }
            }
            refreshLists();
        }
    }
}
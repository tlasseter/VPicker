
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

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
    StsModelSurface[] modelSurfaces = null;
	JPanel mainPanel = new JPanel(new BorderLayout());
	StsDoubleListPanel panel = new StsDoubleListPanel();

    public StsCorrelateMarkers(StsWizard wizard)
    {
		super(wizard);
		setPanels(mainPanel, header);
        mainPanel.add(panel);
    	panel.setLeftListRenderer(new StsColorListRenderer());
    	panel.setRightListRenderer(new StsColorListRenderer());
  //      panel.setTitle("Correlate markers with selected surfaces:");
        panel.setLeftTitle("Selected Surfaces:");
        panel.setRightTitle("Available Markers:");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Correlate selected markers/surfaces");
        panel.getRightButton().setToolTipText("Remove marker correlation with selected surfaces");
//        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Horizon Construction");
        header.setSubtitle("Correlate with Well Markers");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");                                        
        header.setInfoText(wizardDialog,"(1) Select a horizon from left list.\n" +
                                  "(2) Select the related well marker from the right list.\n" +
                                  "   ***** Marker list is generated from all unique marker names from all wells *****\n" +
                                  "   ***** If same marker has multiple names, use marker editor to fix first. *****\n" +
                                  "(3) Press the > button to assign the marker to the horizon or\n" +
                                  "    Press the < button to unassign the marker from the selected horizon.\n" +
                                  "(4) Press the Next>> Button to complete horizon construction.");
    }
	public void setModelSurfaces(StsModelSurface[] modelSurfaces) { this.modelSurfaces = modelSurfaces; }

    public void initSurfaceList()
    {
    	if( modelSurfaces != null )
        {
            int nSurfs = modelSurfaces == null ? 0 : modelSurfaces.length;
            if( nSurfs > 0 )
            {
                StsColorListItem[] items = new StsColorListItem[nSurfs];

                for( int i=0; i<nSurfs; i++ )
                {
                    String surfName = new String(modelSurfaces[i].getName());

                    StsMarker m = modelSurfaces[i].getMarker();
                    if (m == null) surfName += " ( )";
                    else surfName += " (" + m.getName() + ")";

                    StsZone z = modelSurfaces[i].getZoneAbove();
                    if (z != null) surfName += "  :  base of " + z.getName();
                    z = modelSurfaces[i].getZoneBelow();
                    if (z != null) surfName += "  :  top of " + z.getName();

                    items[i] = new StsColorListItem(modelSurfaces[i].getStsColor(), surfName);
                }
                int[] selected = panel.getSelectedLeftIndices();

                // check to see if any of the currently selected have been removed
                if( selected != null )
                {
                	int nSelected = selected.length;
                	for( int i=0; i<selected.length; i++ )
                    	if( selected[i] >= nSurfs ) nSelected--;
                    if( nSelected == 0 ) selected = null;
                    else
                    {
                    	int[] newSelected = new int[nSelected];
                        int iSelected = 0;
                        for( int i=0; i<nSelected; i++ )
	                    	if( selected[i] < nSurfs )
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
		StsClass markers = model.getCreateStsClass(StsMarker.class);
        if(markers == null) return;

        int nMarkers = markers.getSize();
        if(nMarkers == 0) return;

        int nMarkersWithoutSurface = 0;
        StsColorListItem[] items = new StsColorListItem[nMarkers];
        for( int i=0; i<nMarkers; i++ )
        {
            StsMarker marker = (StsMarker)markers.getElement(i);
            if(marker.getName() == null)
            	continue;
            if((marker.getType() == StsMarker.FMI) || (marker.getType() == StsMarker.EQUIPMENT) ||
            		(marker.getType() == StsMarker.PERFORATION))
            	continue;
            items[nMarkersWithoutSurface++] = new StsColorListItem(marker.getStsColor(), marker.getName());
        }

        if(nMarkersWithoutSurface == 0) return;
        items = (StsColorListItem[])StsMath.arraycopy(items, nMarkersWithoutSurface);
        panel.setRightItems(items);
    }

    public boolean start()
    {
		refreshLists();
        // select all of the surfaces
        if(modelSurfaces != null)
        {
        	int[] selected = new int[modelSurfaces.length];
            for( int i=0; i<modelSurfaces.length; i++ ) selected[i] = i;
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
    	initSurfaceList();
    	initMarkerList();
        refreshButtons();
    }

    public String trimName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index-1);
        return null;
    }

    public StsMarker[] getSelectedMarkers()
    {
        StsMarker[] markers = null;

        try
        {
            Object[] items = panel.getSelectedRightItems();
            int nItems = items == null ? 0 : items.length;
            if( nItems == 0) return null;

            markers = new StsMarker[nItems];
            for( int i=0; i<nItems; i++ )
            {
                StsColorListItem item = (StsColorListItem)items[i];
                markers[i] = (StsMarker)model.getObjectWithName(StsMarker.class, item.getName());
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
		Object[] surfaceItems = panel.getSelectedLeftItems();
		Object[] markerItems = panel.getSelectedRightItems();
    	int nSurfaceItems = surfaceItems == null ? 0 : surfaceItems.length;
    	int nMarkerItems = markerItems == null ? 0 : markerItems.length;
		panel.enableLeftButton(nSurfaceItems > 0 && nSurfaceItems == nMarkerItems);

		if(nSurfaceItems > 0 && nMarkerItems == 0)
        {
			for( int i=0; i<nSurfaceItems; i++ )
            {
            	StsColorListItem item = (StsColorListItem) surfaceItems[i];
            	StsModelSurface surface = null;
                try { surface = (StsModelSurface)model.getObjectWithName(StsModelSurface.class, trimName(item.getName())); }
                catch(Exception ex) { }
                if( surface != null && surface.getMarker() != null )
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
			Object[] surfaceItems = panel.getSelectedLeftItems();
			Object[] markerItems = panel.getSelectedRightItems();
	    	int nSurfaceItems = surfaceItems == null ? 0 : surfaceItems.length;
    		int nMarkerItems = markerItems == null ? 0 : markerItems.length;
            if( nSurfaceItems == nMarkerItems )
            {
//		        System.out.println("Correlating " + nSurfaceItems + " surfaces");
                for( int i=0; i<nSurfaceItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem) surfaceItems[i];
                    StsModelSurface surface = null;
                    try { surface = (StsModelSurface)model.getObjectWithName(StsModelSurface.class, trimName(item.getName())); }
                    catch(Exception ex) { }

                    item = (StsColorListItem) markerItems[i];
                    StsMarker marker = (StsMarker)model.getObjectWithName(StsMarker.class, item.getName());
                    if( surface != null && marker != null )
                    {
                    	surface.setMarker(marker);
                    }
                    // Might have multiple surface markers in a single well
                    StsObjectList.ObjectIterator iter = model.getObjectIterator(StsMarker.class);
                    while(iter.hasNext())
                    {
                        marker = (StsMarker)iter.next();
                        if(marker.getName().equals(item.getName()))
                            marker.setModelSurface(surface);
                    }
                }
                refreshLists();
            }
        }
        else if( e.getActionCommand().equals(">") )
        {
			Object[] surfaceItems = panel.getSelectedLeftItems();
	    	int nSurfaceItems = surfaceItems == null ? 0 : surfaceItems.length;
//            System.out.println("Uncorrelating " + nSurfaceItems + " surfaces");
            for( int i=0; i<nSurfaceItems; i++ )
            {
                StsColorListItem item = (StsColorListItem) surfaceItems[i];
                StsModelSurface surface = null;
                try { surface = (StsModelSurface)model.getObjectWithName(StsModelSurface.class, trimName(item.getName())); }
                catch(Exception ex) { }
                if( surface != null )
                {
                	StsMarker oldMarker = surface.getMarker();
                    if( oldMarker != null ) oldMarker.setModelSurface(null);
                    surface.setMarker(null);
                }
            }
            refreshLists();
        }
    }
}

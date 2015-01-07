
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Zones;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsDefineZones extends StsWizardStep
{
    StsTablePanel panel;
    StsHeaderPanel header;

    StsModelSurface[] surfaces = null;
    StsZone[] zones = null;

    public StsDefineZones(StsWizard wizard)
    {
        super(wizard, new StsTablePanel(), null, new StsHeaderPanel());
        panel = (StsTablePanel) getContainer();
        panel.setTitle("Select horizon pairs to define units:");
        panel.addColumns(new Object[] { "Zone", "Top", "Base" });
        panel.setPreferredSize(new Dimension(400,300));
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Model Construction");
        header.setSubtitle("Select Zones");
        header.setInfoText(wizardDialog,"(1) Select the top-base horizon pairs to include in model.\n" +
                                 "   ***** The list is automatically generated from all horizons in database ***** \n" +
                                 "(2) Press the Next>> Button");
    }

    public boolean start()
    {
        surfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
        if(surfaces == null || surfaces.length < 2)
        {
            logMessage("Can't build model.  Less than two model horizons available."
                    +  "  Terminating action...");
            return false;
        }
        int nSurfaces = surfaces.length;

        int nZones = nSurfaces - 1;
        int[] selected = new int[nZones];
        Object[] row = new Object[3];
        String botSurfaceName = surfaces[0].getName();
        for( int n = 0; n < nZones; n++)
        {
            String topSurfaceName = botSurfaceName;
            botSurfaceName = surfaces[n+1].getName();
            row[0] = Integer.toString(n);
            row[1] = topSurfaceName;
            row[2] = botSurfaceName;
            panel.addRow(row);
            selected[n] = n;
        }
        panel.setSelectedIndices(selected);
        return true;
    }

    public boolean end() { return getSelectedZones() != null; }

    public StsZone[] getSelectedZones()
    {
        try
        {
            if(zones == null)
            {
                int[] selected = panel.getSelectedIndices();
                int nItems = selected == null ? 0 : selected.length;
                if(nItems > 0)
                {
                    model.deleteStsClass(StsZone.class);

                    zones = new StsZone[nItems];
                    for(int n = 0; n < nItems; n++)
                    {
                        StsModelSurface top = surfaces[selected[n]];
                        StsModelSurface base = surfaces[selected[n]+1];
					/*
                        StsZone existingZone = zoneClass.getExistingZone(top, base);
                        if(existingZone != null)
                            zones[n] = existingZone;
	                */
			//		 else
                        {
                            String name = top.getOriginalSurface().getName();
                            zones[n] = createZone(name, top, base);
                        }
                    }
                }
                panel.setEnabled(false);
            }
            return zones;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSelectZones.getSelectedZones() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	private StsZone createZone(String name, StsModelSurface topModelSurface, StsModelSurface baseModelSurface)
    {
    	// create the zone from the imported surfaces ??
    	StsZone zone = null;
    	try
        {
            zone = getMatchingZone(name, topModelSurface, baseModelSurface);
            if(zone != null) return zone;
            zone = new StsZone(name, topModelSurface, baseModelSurface);
			/*
            zone.setTopModelSurface(topModelSurface);
            topModelSurface.setZoneBelow(zone);
            zone.setBaseModelSurface(baseModelSurface);
            baseModelSurface.setZoneAbove(zone);
	*/
        /*
            zone.buildWellZoneSet(model);
            StsWellZoneSet wzs = zone.getWellZoneSet();
            if (wzs != null)
            {
                StsWellZoneSet.buildWellZones(model, wzs.getTopMarker(), wzs.getBaseMarker());
            }
        */
            StsMarker.tryToLinkIntervalBelow(topModelSurface.getMarker(), baseModelSurface);

        }
        catch(Exception e) { e.printStackTrace(); }
        return zone;
    }

    private StsZone getMatchingZone(String name, StsModelSurface topModelSurface, StsModelSurface baseModelSurface)
    {
        StsObject[] zoneObjects = model.getObjectList(StsZone.class);
        int nZones = zoneObjects.length;
        StsZone matchingZone = null;
        for(int n = 0; n < nZones; n++)
        {
            StsZone zone = (StsZone)zoneObjects[n];
            if(zone.matches(name, topModelSurface, baseModelSurface))
            {
                matchingZone = zone;
                break;
            }
        }
        return matchingZone;
    }
}
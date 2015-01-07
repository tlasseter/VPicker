
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;

public class StsBuildZonesFromWellZones extends StsAction
{
 	public StsBuildZonesFromWellZones(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            // retrieve a well zone sets for zone creation
            StsClass wellZoneSets = model.getCreateStsClass(StsWellZoneSet.class);
            if (wellZoneSets==null)
            {
                logMessage("No well zone sets found.");
                return true;
            }
            int nWellZoneSets = wellZoneSets.getSize();
            if (nWellZoneSets==0)
            {
                logMessage("No well zone sets found.");
                return true;
            }

            String[] names = new String[nWellZoneSets];
            StsColor[] colors = new StsColor[nWellZoneSets];
            for (int i=0; i<nWellZoneSets; i++)
            {
                StsWellZoneSet wellZoneSet = (StsWellZoneSet)wellZoneSets.getElement(i);
                names[i] = wellZoneSet.getName();
                StsObject[] zoneLinks = { wellZoneSet.getParentZone() };
                if (zoneLinks[0]==null) names[i] += "   ( )";
                else names[i] += "   (" + zoneLinks[0].getName() + ")";
                colors[i] = wellZoneSet.getStsColor();
            }
            StsColorListSelector selector = new StsColorListSelector(new JFrame(),
                    model.getName(), "Well Zones  (with zone links)", colors, names);
            selector.setSingleSelectionMode(false);
            selector.setVisible(true);
            int[] selectedIndices = selector.getSelectedIndices();
            if (selectedIndices==null)
            {
                logMessage("No well zone sets selected.");
    	        actionManager.endCurrentAction();
                return true;
            }

            int nZonesCreated = 0;
            for (int i=0; i<selectedIndices.length; i++)
            {
                StsWellZoneSet wellZoneSet =
                        (StsWellZoneSet)wellZoneSets.getElement(selectedIndices[i]);
                StsZone zone = wellZoneSet.getParentZone();
                if (zone != null)
                {
                    logMessage("Well zone set " + wellZoneSet.getName() +
                            " is already linked with zone " + zone.getName() +
                            ".  Continuing...");
                    continue;
                }
                zone = new StsZone(wellZoneSet);
                logMessage("Well zone set " + wellZoneSet.getName() +
                            " has created zone " + zone.getName());
                nZonesCreated++;
            }

            if (nZonesCreated == 0)
            {
                StsMessageFiles.logMessage("No new zones built.");
            }

            else if (nZonesCreated==1)
            {
                StsMessageFiles.logMessage("Successfully built 1 new zone.");
            }
            else
            {
                StsMessageFiles.logMessage("Successfully built " + nZonesCreated
                        + " new zones.");
            }
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildZonesFromWellZones failed.", e, StsException.WARNING);
            return false;
        }
        return true;
    }
}

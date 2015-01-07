
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
import com.Sts.RFUI.*;
import com.Sts.UI.*;

public class StsLinkZoneWithWellZoneSet extends StsAction
{
 	public StsLinkZoneWithWellZoneSet(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            // check special cases
            StsClass zoneList = model.getCreateStsClass(StsZone.class);
            int nZones = zoneList.getSize();
            if (nZones<1)
            {
                logMessage("No zones found to link.");
                return true;
            }
            StsClass wellZoneSetList = model.getCreateStsClass(StsWellZoneSet.class);
            int nWellZoneSets = wellZoneSetList.getSize();
            if (nWellZoneSets<1)
            {
                logMessage("No well zones sets found to link.");
                return true;
            }

            while (true)  // keep asking until user cancels
            {
                // retrieve zones and select one
                StsWellZoneLinkDialog selector = StsWellZoneLinkDialog.constructor(model);
				if(selector == null) return false;
                StsZone zone = (StsZone)selector.selectObject();
                if (zone==null)   // user cancelled
                {
                	actionManager.endCurrentAction();
                    return true;
                }

                logMessage("Selected zone:  " + zone.getName());
                StsWellZoneSet currentWellZoneSet = zone.getWellZoneSet();

                // retrieve well zone sets and select one
                StsZoneLinkDialog selector2 = StsZoneLinkDialog.constructor(model,
                        "Zone: " + zone.getName());
				if(selector2 == null) return false;
                if (currentWellZoneSet!=null)
                {
                    StsListSelector listSelector = selector2.getSelector();
                    if (listSelector!=null)
                    {
                        listSelector.setSelectedItem(currentWellZoneSet.getName());
                    }
                }
                StsWellZoneSet newWellZoneSet = (StsWellZoneSet)selector2.selectObject();
                if (newWellZoneSet==null) continue;

                logMessage("Selected well zone set:  " + newWellZoneSet.getName());
                if (newWellZoneSet!=currentWellZoneSet)
                {
                    zone.setWellZoneSet(newWellZoneSet);
                    newWellZoneSet.setParentZone(zone);
                    logMessage("Linked zone: " + zone.getName() +
                            " with well zone set:  " + newWellZoneSet.getName());
                }
                else
                {
                    logMessage("Previously linked zone: " + zone.getName() +
                            " with well zone set:  " + newWellZoneSet.getName());
                }
                if (nZones==1 || nWellZoneSets==1)  // user can't cancel!
                {
                	actionManager.endCurrentAction();
                    return true;
                }
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsLinkZoneWithWellZoneSet.start()\n" + e);
			actionManager.endCurrentAction();
            return false;
        }
    }
}

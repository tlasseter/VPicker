
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

public class StsBuildWellZonesFromMarkers extends StsAction
{
 	public StsBuildWellZonesFromMarkers(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            StsClass markers = model.getCreateStsClass(StsMarker.class);
            if (markers==null)
            {
                logMessage("No markers found." +
                        "  At least two are required to build a well zone.");
    	        actionManager.endCurrentAction();
                return true;
            }
            int nMarkers = markers.getSize();
            if (nMarkers<2)
            {
                logMessage("At least two markers are required to build a well zone.");
    	        actionManager.endCurrentAction();
                return true;
            }
            String[] markerNames = null;
            while (true)
            {
    		    StsWellZoneBuilderDialog dialog = new StsWellZoneBuilderDialog(model.win3d);
                dialog.setLocationRelativeTo(model.win3d);
                dialog.setVisible(true);
                markerNames = dialog.getSelectedItems();
                if (markerNames==null)  // user cancelled
                {
                    logMessage("Marker selection cancelled.");
			        actionManager.endCurrentAction();
                    return true;
                }
                if (markerNames.length==1)
                {
                    logMessage("Only 1 marker chosen.  Please reselect.");
                    dialog.dispose();
                    continue;
                }
                break;
            }
            markerNames = StsOrderedListDialog.parseNames(markerNames);
            int nWellZoneSets = buildWellZoneSets(markerNames);
            if (nWellZoneSets<1)
            {
                logMessage("Error! Unable to build well zone sets.");
    	        actionManager.endCurrentAction();
                return false;
            }
            if (nWellZoneSets==1)
            {
                logMessage("Successfully built 1 well zone set.");
            }
            else
            {
                logMessage("Successfully built " + nWellZoneSets
                        + " well zone sets.");
            }

			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsBuildWellZonesFromMarkers.start()\n" + e);
            return false;
        }
        return true;
    }

    /** build well zone sets from markers (bottom to top) */
    public int buildWellZoneSets(String[] markerNames)
    {
        if (markerNames == null || markerNames.length < 2) return 0;
        int nWellZoneSetsBuilt = 0;

        StsMarker top = null;
        StsMarker base = (StsMarker)model.getObjectWithName(StsMarker.class, markerNames[0]);
        for (int i=1; i<markerNames.length; i++)
        {
            top = base;
            base = (StsMarker)model.getObjectWithName(StsMarker.class, markerNames[i]);
            if (!StsWellZoneSet.isNewWellZoneSet(model, top, base)) continue;
            StsWellZoneSet wzs = null;
            try { wzs = new StsWellZoneSet(StsWellZoneSet.STRAT, top, base); }
            catch (Exception e) { return nWellZoneSetsBuilt; }
            wzs.buildZone(model);
            StsMarker.tryToLinkIntervalBelow(top, top.getModelSurface());
            if (!StsWellZoneSet.buildWellZones(model, top, base))
            {
                return nWellZoneSetsBuilt;
            }
            nWellZoneSetsBuilt++;
        }
        return nWellZoneSetsBuilt;
    }


}


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
import com.Sts.Utilities.*;

public class StsEdgesFromImportedEdges extends StsAction implements Runnable
{
    private boolean actionOK = false;

   	/** convert REFERENCE StsSectionEdges to FAULT StsSectionEdges */
    public StsEdgesFromImportedEdges(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

	public void run()
    {

		if(!actionOK()) return;

        try
        {
            model.disableDisplay();

            statusArea.setTitle("Building Fault Cuts from Imported cuts: ");

            // since we may be deleting and adding to the sectionEdges, get a list first
            StsClass sectionEdges = model.getCreateStsClass(StsSectionEdge.class);
            StsObject[] edgesList = sectionEdges.getElements();
            int nSectionEdges = edgesList.length;

            int nEdgesConverted = 0;
            for(int n = 0; n < nSectionEdges; n++)
            {
                StsSectionEdge edge = (StsSectionEdge)edgesList[n];
                if(edge.getType() == StsParameters.REFERENCE) edge.convertFromRefToFault();
            }
            logMessage(nEdgesConverted + " imported edges converted to faults.");

            model.enableDisplay();
            model.win3d.win3dDisplay();
            actionManager.endCurrentAction();
//            model.glPanel.getActionManager().fireChangeEvent();
        }
        catch(Exception e)
        {
            actionOK = false;
            return;
        }
    }

	private boolean actionOK()
    {
    	try
        {
            actionOK = true;

            if(model.hasObjectsOfType(StsModelSurface.class, StsModelSurface.MODEL))
            {
                logMessage("No horizons have been built."
                        +  "  Terminating action.");
                actionOK = false;
            }
	        return actionOK;
        }
        catch(Exception e)
        {
            actionOK = false;
            return actionOK;
        }
    }

    public boolean end()
    {
        if (!actionOK) logMessage("Fault cuts from imported cuts failed.");
        return actionOK;
    }
}


//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;

public class StsImportFaultCuts extends StsAction implements Runnable
{
    private static final byte STS = 1;
    private static final byte CPS3 = 2;
    private static final byte OPENWORKS = 3;
    private static final String[] descriptions = new String[] { "STS", "CPS3", "OpenWorks" };

	boolean success = false;

 	public StsImportFaultCuts(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            statusArea.setTitle("Read Fault Cuts:");
			success = getFaultCuts();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsImportFaultCuts.start()\n" + e);
            success = false;
        }
        actionManager.endCurrentAction();
    }

    public boolean end()
    {
        if (!success) logMessage("Reading of fault cuts failed!");
    	return success;
    }

    public boolean getFaultCuts()
    {
        try
        {
            // get the data directory
            String dirPath = ".";
            try { dirPath = model.getProject().getDataFullDirString(); }
            catch (Exception e) { }

            // build the file selector
            StsFileChooser fileChooser = new StsFileChooser(model.win3d, "Select fault cut ASCII file(s):", dirPath);

            // pop up the selector and see what was selected
	    	if (!fileChooser.show())
            {
		        logMessage("No fault cuts file(s) selected.  Terminating...");
    			actionManager.endCurrentAction();
                return true;
            }

            // get filenames
            String[] filenames = fileChooser.getFilenames();

            // get dirPath (could have been changed by user)
            dirPath = fileChooser.getDirectoryPath();

            return getFaultCuts(dirPath, filenames, fileChooser.getFileFilterIndex());
        }
        catch (Exception e)
        {
            logMessage("Fault cut import failed!");
            StsException.outputException("StsImportFaultCuts.getFaultCuts: ",
                    e, StsException.FATAL);
            return false;
        }
    }

    public boolean getFaultCuts(String dirPath, String[] filenames, int type)
    {
        if (model==null) return false;
        StsProject project = model.getProject();
        if (project == null) return false;
        if (!model.hasObjectsOfType(StsModelSurface.class, StsModelSurface.MODEL)) return false;

        if (dirPath==null || filenames==null)
        {
		    logMessage("No fault cut file(s) specified.  Terminating...");
			actionManager.endCurrentAction();
            return true;
        }

        int nFilesLoaded = 0;
        int nCutsLoaded = 0;
        try
        {
            // set up status area
            statusArea.setMaximum(filenames.length);
            statusArea.setProgress(0);
            statusArea.addProgress();
            statusArea.sleep(10);

            // go thru files
            for (int i=0; i<filenames.length; i++)
            {
        	    // open the fault cut file (use Open Works fault polygon format for now)
                File f = new File(dirPath, filenames[i]);
                StsFaultCutFile fcf = null;
                switch(type)
                {
                case OPENWORKS:
                    fcf = new StsOWFaultCutFile(true, f.getAbsolutePath());
                    break;

                case CPS3:
                    fcf = new CPS3FaultCutFile(true, f.getAbsolutePath());
                    break;
                }


                if (fcf == null)
                {
	        	    logMessage("File " + filenames[i] +
                            ": not found.  Continuing...");
                    statusArea.setProgress(i+1); // update progress
                    continue;
                }

                // read the cuts
                StsFaultCut[] fcs = null;
                try
                {
                    fcs = fcf.getFaultCuts();
                    if (fcs == null || fcs.length == 0)
                    {
                        logMessage("No fault cuts read from file " + filenames[i] +
                                ".  Continuing...");
                        statusArea.setProgress(i+1); // update progress
                        continue;
                    }
                }
                catch(StsException Stse)
                {
                	Toolkit.getDefaultToolkit().beep();
                    logMessage(Stse.getMessage() + " Continuing...");
                    continue;
                }

                // select a horizon to associate cuts with
                StsSelectStsObjects selector = StsSelectStsObjects.constructor(model,
                    model.getObjectListOfType(StsModelSurface.class, StsModelSurface.MODEL), filenames[i],
                    "Select the corresponding horizon:", true);
				if(selector == null) return false;
                StsModelSurface horizon = (StsModelSurface)selector.selectObject();
                if (horizon == null)
                {
                    logMessage("No horizon selected: skipping fault cuts...");
                    statusArea.setProgress(i+1); // update progress
                    continue;
                }
//                model.setViewProperty("displaySections", false); // turn sections off until we have built new edges
                StsSectionEdge[] newEdges = new StsSectionEdge[fcs.length];
                byte edgeType = type == STS ? StsSection.FAULT : StsSectionEdge.REFERENCE;
                for (int j=0; j<fcs.length; j++)
                {
                    newEdges[j] = new StsSectionEdge(edgeType, null, horizon, 1);
                    logMessage("Setting vertices for fault cut " + j);
                    newEdges[j].constructFromFaultCut(fcs[j], horizon);
                }

                logMessage("Added " + fcs.length + " fault cuts to surface "
                        + horizon.getName());
                statusArea.setProgress(i+1); // update progress

                nCutsLoaded += fcs.length;
                nFilesLoaded++;
            }
            logMessage("Loaded:  " + nCutsLoaded + " fault cuts " +
                    "from " + nFilesLoaded + descriptions[type-1] + " ASCII files.");
            statusArea.textOnly();
        }
        catch(Exception e)
        {
            statusArea.removeProgress();
            logMessage("Loaded:  " + nCutsLoaded + " fault cuts " +
                    "from " + nFilesLoaded + " OpenWorks ASCII files.");
            StsException.outputException("StsImportFaultCuts.getFaultCuts: ",
                    e, StsException.FATAL);
        }
        ((StsSectionClass)model.getCreateStsClass(StsSection.class)).setDisplaySections(true); // turn sections off until we have built new edges
        model.win3d.win3dDisplay();
        return true;
    }
}

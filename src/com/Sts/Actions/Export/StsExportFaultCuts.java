
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Export;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsExportFaultCuts extends StsAction implements Runnable
{
	boolean success = false;

 	public StsExportFaultCuts(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            statusArea.setTitle("Write Fault Cuts:");
			success = putFaultCuts();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsExportFaultCuts.start()\n" + e);
            success = false;
        }
        actionManager.endCurrentAction();
    }

    public boolean end()
    {
        if (!success) logMessage("Writing of fault cuts failed!");
    	return success;
    }

    private boolean putFaultCuts()
    {
        try
        {
            // get the data directory
            String path = ".";
            try { path = model.getProject().getDataFullDirString(); }
            catch (Exception e) { }

            // get a file to save to
	        StsFileChooser chooser = StsFileChooser.createFileChooserPostfix(model.win3d, "Save fault cuts to Ascii file", path, ".*");
            String filename = null;
			if (chooser.showSave())
            {
                filename = chooser.getFilePath();
                File file = new File(filename);
                if (file.exists())
                {
                    StsYesNoPanel yesNo = new StsYesNoPanel("Overwrite existing file?",
                        "Yes", "No");
                    yesNo.setVisible(true);
                    if (yesNo.isNo())
                    {
                        statusArea.textOnly();
            		    logMessage("Previous file not replaced.  Terminating...");
//			            glPanel.getActionManager().endCurrentAction();
                        return true;
                    }
                }
            }
            else
            {
                statusArea.textOnly();
            	logMessage("No fault cut file specified.  Terminating...");
//			    glPanel.getActionManager().endCurrentAction();
                return true;
            }

            return putFaultCuts(filename);
        }
        catch (Exception e)
        {
            logMessage("Fault cut Export failed!");
            StsException.outputException("StsExportFaultCuts.getFaultCuts: ",
                    e, StsException.FATAL);
            return false;
        }
    }

    private boolean putFaultCuts(String filename)
    {
		int nCuts = 0;
        StsFaultCutFile fcf = null;

		StsObjectList.ObjectIterator faults = model.getObjectOfTypeIterator(StsSection.class, StsParameters.FAULT);
		int nFaults = faults.getSize();
		if(nFaults == 0) return false;

        try
        {
            fcf = new StsOWFaultCutFile(false, filename);
			fcf.open();

            // set up status area
            statusArea.setMaximum(nFaults);
            statusArea.setProgress(0);
            statusArea.addProgress();
            statusArea.sleep(10);

            // go thru section edges and write fault cuts not in boundaries
            sectionEdgeLoop:
            while(faults.hasNext())
            for (int i=0; i < nFaults; i++)
            {
                StsSection fault = (StsSection)faults.next();
                StsObjectRefList sectionEdges = fault.getSectionEdges();
				int nEdges = sectionEdges.getSize();
				for(int n = 0; n < nEdges; n++)
				{
                    StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(n);

					// write the cuts
					if (!fcf.writeFaultCut(edge.convertToFaultCut()))
					{
						logMessage("Unable to write fault cut to file " + filename +
								".  Continuing...");
					}
					else nCuts++;

					statusArea.setProgress(i+1); // update progress
				}
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsExportFaultCuts.putFaultCuts() failed. ",
                    e, StsException.WARNING);
        }
		finally
		{
			try {  fcf.close(); } catch(Exception e) { }
            logMessage("Saved:  " + nCuts + " fault cuts " +
                    "to OpenWorks ASCII file:  " + filename);
            statusArea.textOnly();
		}
        return true;
    }
}

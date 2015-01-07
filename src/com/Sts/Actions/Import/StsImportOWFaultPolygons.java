
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

import java.io.*;

public class StsImportOWFaultPolygons extends StsAction implements Runnable
{
	boolean success = false;

 	public StsImportOWFaultPolygons(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            if(!model.hasObjectsOfType(StsModelSurface.class, StsModelSurface.IMPORTED))
            {
                logMessage("Unable to import fault polygons before "
                        + "surfaces are imported.  Terminating action.");
                success = false;
                actionManager.endCurrentAction();
                return;
            }

            statusArea.setTitle("Read Polygon File:");
			success = getPolygons();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsImportOWFaultPolygons.start()\n" + e);
            success = false;
        }
        actionManager.endCurrentAction();
    }

    public boolean end()
    {
        if (!success) logMessage("Reading of fault polygons failed!");
    	return success;
    }

    public boolean getPolygons()
    {
        try
        {
            // get the data directory
            String path = ".";
            try { path = model.getProject().getDataFullDirString(); }
            catch (Exception e) { }

            // build the file selector
            StsFileChooser fileChooser = StsFileChooser.createMultiFileChooserPrefix(model.win3d, "Select OpenWorks fault polygon ASCII file(s):", path, null);
//            StsFileChooser fileChooser = new StsFileChooser(model.win3d, "Select OpenWorks fault polygon ASCII file(s):", path, null, true, true);

            // pop up the selector and see what was selected
	    	if (!fileChooser.show())
            {
		        logMessage("No fault polygon file selected.  Terminating...");
    			actionManager.endCurrentAction();
                return true;
            }

            // get filenames
            String[] filenames = fileChooser.getFilenames();

            // get path (could have been changed by user)
            path = fileChooser.getDirectoryPath();

            return getPolygons(path, filenames);
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("OpenWorks fault polygon import failed!");
            StsException.outputException("StsOWImportFaultPolygons.getPolygons: ",
                    e, StsException.FATAL);
            return false;
        }
    }

    public boolean getPolygons(String path, String[] filenames)
    {
        if (model==null) return false;
        StsProject project = model.getProject();
        if (project == null) return false;
        if (!model.hasObjectsOfType(StsModelSurface.class, StsModelSurface.IMPORTED)) return false;

        if (path==null || filenames==null)
        {
		    StsMessageFiles.logMessage("No fault polygon files specified.  Terminating...");
			actionManager.endCurrentAction();
            return true;
        }

        int nFilesLoaded = 0;
        int nPolygonsLoaded = 0;
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
        	    // open the polygon file
                File f = new File(path, filenames[i]);
                StsFaultPolygonFile fpf = new StsFaultPolygonFile(f.getAbsolutePath(),
                        StsFaultPolygonFile.OPEN_WORKS);
                if (fpf == null)
                {
	        	    logMessage("File " + filenames[i] +
                            ": not found.  Continuing...");
                    statusArea.setProgress(i+1); // update progress
                    continue;
                }

                // read the polygons
                StsFaultPolygon[] fps = fpf.getFaultPolygons();
                if (fps == null)
                {
	        	    logMessage("No fault polygons read from file " + filenames[i] +
                            ".  Continuing...");
                    statusArea.setProgress(i+1); // update progress
                    continue;
                }

                // select a surface to associate polygons with
                StsSelectStsObjects selector = StsSelectStsObjects.constructor(model,
                    model.getObjectListOfType(StsModelSurface.class, StsModelSurface.IMPORTED), filenames[i],
                    "Select the corresponding surface:", true);
				if(selector == null) return false;
                StsModelSurface s = (StsModelSurface)selector.selectObject();
                if (s == null)
                {
                    logMessage("No surface selected: skipping fault polygons...");
                    continue;
                }
                s.addFaultPolygons(fps);
                logMessage("Added " + fps.length + " fault polygons to surface "
                        + s.getName());
                statusArea.setProgress(i+1); // update progress

                nPolygonsLoaded += fps.length;
                nFilesLoaded++;
            }
            logMessage("Loaded:  " + nPolygonsLoaded + " fault polygons " +
                    "from " + nFilesLoaded + " OpenWorks ASCII files.");
            statusArea.textOnly();
        }
        catch (Exception e)
        {
            statusArea.removeProgress();
            logMessage("Loaded:  " + nPolygonsLoaded + " fault polygons " +
                    "from " + nFilesLoaded + " OpenWorks ASCII files.");
            StsException.outputException("StsImportOWFaultPolygons.getPolygons: ",
                    e, StsException.FATAL);
        }

        model.win3d.win3dDisplay();
        return true;
    }
}

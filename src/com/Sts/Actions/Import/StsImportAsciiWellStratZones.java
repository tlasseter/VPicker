
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsImportAsciiWellStratZones extends StsAction
{

   	/** import well zones from an ascii file */
 	public StsImportAsciiWellStratZones(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
			create();
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            StsException.outputException("StsImportAsciiWellsStratZones.start() failed.",
                e, StsException.WARNING);
            return false;
        }
        return true;
    }

    /** read in multiple well zones from a list of Ascii files chosen in a dialog */
    public void create()
    {

        if (model.getCreateStsClass(StsWell.class).getSize() == 0)
        {
            logMessage("Unable to load strat zones:  no wells found.");
            return;
        }

        // get the data directory
        String path = ".";
        try { path = model.getProject().getDataFullDirString(); }
        catch (Exception e) { }

        // build the file selector
        StsFileChooser fileChooser = StsFileChooser.createMultiFileChooserPrefix(model.win3d, "Select 1 or more well zone files:", path, "zone.strat");

        // pop up the selector and see what was selected
		if (!fileChooser.show()) return;  // nothing selected

        // get filenames
        String[] filenames = fileChooser.getFilenames();

        // get path (could have been changed by user)
        path = fileChooser.getDirectoryPath();

        // turn on the wait cursor

        logMessage("Preparing to load " + filenames.length
                + " wells with Ascii strat zones ...");

        // read the files and build the well entities
        int nLoaded = create(path, filenames);

        logMessage("Loaded " + nLoaded + " Ascii strat zones ...");
    }

    /** read in multiple well zones from a list of Ascii files */
    public int create(String path, String[] filenames)
    {
        int totalLoaded = 0;

        try
        {
            if (path==null || filenames==null) return totalLoaded;

            for (int i=0; i<filenames.length; i++)
            {
                // extract the well name
                StsWellKeywordIO.setParseFilename(filenames[i], "zone", "strat");
                String wellname = StsWellKeywordIO.getWellName();

                // be sure we've already got this well before proceeding
                StsWell well = (StsWell)model.getObjectWithName(StsWell.class, wellname);
                if (well==null) continue;

                // read well strat zone file
                int nLoaded = StsWellKeywordIO.constructWellStratZones(well, path, filenames[i]);
                totalLoaded += nLoaded;
            }
            return totalLoaded;
        }
        catch(Exception e)
        {
            StsException.outputException("StsImportAsciiWellStratZones.create() failed.",
                e, StsException.WARNING);
            return totalLoaded;
        }
    }
}

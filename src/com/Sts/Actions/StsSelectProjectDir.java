
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsSelectProjectDir extends StsAction
	implements Runnable
{
    private boolean success = false;

   	/** create a new section and attach a mouse listener to it */
 	public StsSelectProjectDir(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public void run()
    {
        File projectDir = new File(".");
        try
        {
            StsFileChooser f = new StsFileChooser(model.win3d, "Choose a root project directory:",
                                                  projectDir.getPath());
            boolean selected = f.show();

            String[] filenames = f.getFilenames();
            String filename = null;
            if (filenames != null) filename = filenames[0];

            if( filename != null )
                projectDir = new File(f.getDirectoryPath(), filename);
            else
                projectDir = new File(f.getDirectoryPath());

            if (projectDir == null || !projectDir.exists())
            {
                logMessage("Invalid project directory selected.");
                return;
            }

            model.getProject().setRootDirectory(projectDir);
            success = true;
        }
        catch(Exception exc)
        {
            StsException.outputException("StsSelectProjectDir.run() failed.", exc, StsException.FATAL);
        }
   	}

    public boolean end()
    {
        if (success) logMessage("Project directory selected.");
        else logMessage("Project directory selection failed.");
        return success;
    }
}




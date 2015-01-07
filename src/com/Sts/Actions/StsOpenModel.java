
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsOpenModel extends StsAction // implements Runnable
{
    private boolean openSuccess = true;
    private String dbPathname = null;

     public StsOpenModel(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

     public StsOpenModel(StsActionManager actionManager, String dbPathname)
    {
        super(actionManager, true);
        this.dbPathname = dbPathname;
        model.win3d.rebuildMenubar();
    }

    public void run()
    {
        start();
    }

    public boolean start()
    {
        final StsModel oldModel = model;
        StsAbstractFile file;

        try
        {
            if(oldModel != null)
            {
                actionManager.endCurrentAction();
                oldModel.stopTime();
                oldModel.close();
                String message = "Closed old model : " + oldModel.getName();
                StsMessageFiles.logMessage(message);
                if(Main.isDbCmdDebug) StsException.systemDebug(this, "start", message);
            }

            if(dbPathname == null)
            {
                dbPathname = openExistingDB("." + File.separator);
                if(Main.isDbCmdDebug) StsException.systemDebug(this, "start", "opening db " + dbPathname);
                if(dbPathname == null) return false;
            }

            file = StsFile.constructor(dbPathname);
            if(file == null) return false;
            if(!StsDBFileModel.fileOK(file, null)) return false;

            openSuccess = false;
            StsModel newModel = StsModel.constructor(file);
            if(newModel == null) return false;

            openSuccess = true;
            return openSuccess;
        }
        catch(Exception e)
        {
            StsException.outputException("StsOpenModel.run() failed.", e, StsException.WARNING);
            return false;
        }
        finally
        {
            statusArea.textOnly();
        }
   	}

    private String openExistingDB(String dirPath)
    {
        String filename = "";

        try
        {
            StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(null, "Open model file", dirPath, "db.");
            while(true)
            {
                if(!chooser.show())
                {
                    model.mainWindowActionManager.endCurrentAction();
                    return null;
                }
                else
                {
                    File f = new File(chooser.getFilePath());
                    if( !f.exists() )
                        new StsMessage(null, StsMessage.ERROR, "Must select an existing db file.");
                    else if( f.isDirectory() )
                        new StsMessage(null, StsMessage.ERROR, "Selected file " + f.getPath() + " is a directory.");
                    else
                    {
                        try
                        {
                            return chooser.getFilePath();
                        }
                        catch(Exception e)
                        {
                            new StsMessage(null, StsMessage.WARNING, "Couldn't open db file: " + filename);
                            return null;
                        }
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsOpenExistingDB() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean end()
    {
        return openSuccess;
    }

}




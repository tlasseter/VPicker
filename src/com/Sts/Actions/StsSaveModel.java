
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
import com.Sts.Utilities.StsException;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class StsSaveModel extends StsAction implements Runnable
{
	public static final int SAVE = 0;
	public static final int SAVE_AS = 1;

    private boolean saveSuccess = false;
    private int type;
    private String dbStemname;

    static final String DB_GROUP = StsNewModel.DB_GROUP;

 	public StsSaveModel(StsActionManager actionManager, Integer type)
    {
        super(actionManager, true);
        this.type = type.intValue();
    }

   	public void run()
    {
       	// save as
        String dirPath;
        String filePath = null;
        StsDBFileModel db = model.getDatabase();
        StsDBFileModel newDb = null;

        StsProject project = model.getProject();
        do
        {

            dirPath = project.getModelFullDirString();
            StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(model.win3d, "Specify Database name ", dirPath, StsNewModel.DB_GROUP);
            int option = chooser.showReturnOption();
            if(option != JFileChooser.APPROVE_OPTION) return;
            filePath = selectDatabasePathname(chooser);
        }
        while (filePath == null);

        if(filePath == null) return;

        File f = new File(filePath);
        if(f.exists())
        {
            if(filePath.endsWith(model.getName()))
            {
                new StsMessage(model.win3d, StsMessage.WARNING,
                    "Project is same as current project,\n\n" +
                    "     Projects are automatically saved as they are worked on.\n" +
                    "     There is no need to save under same name.\n");
                    return;
            }
            else
            if (StsYesNoDialog.questionValue(model.win3d, "Project already exist,\n\n" +
                "        Do you want to continue and overwrite \n" +
                "        " + filePath + "?\n"))
                f.delete();
            else
                return;

        }
        try
        {
            StsFile file = StsFile.constructor(filePath);
            newDb = StsDBFileModel.openWrite(file, null);
            if(db != null) db.copyTo(newDb);
            model.setDatabase(newDb);
            project.setDatabaseInfo(newDb);
         }
        catch(Exception e)
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "Couldn't open db file: " + filePath);
            return;
        }
        if(newDb != null)
        {
            try
            {
                // change model name to match output file
                model.setName(buildName(newDb.getFilename()));
                 model.resetWindowTitles(dbStemname);
                statusArea.addProgress();
             /*
	            StsTimer overallTimer = new StsTimer();
	            StsTimer timer = new StsTimer();
                overallTimer.start();
        	    overallTimer.stopPrint(" Model saved in " );
        	*/
                saveSuccess = true;
           	}
           	catch(Exception exc)
           	{
				System.out.println("Exception in save database\n" + exc);
           	}
        }
        actionManager.endCurrentAction();
   	}

    public String selectDatabasePathname(StsFileChooser chooser)
    {
        try
        {
            String filePathname = chooser.getFilePath();
            if(filePathname == null) return null;

            String filename = chooser.getFilename();
            if(!filename.startsWith(StsNewModel.DB_GROUP))
            {
                dbStemname = filename;
                filename = DB_GROUP + filename;
            }
            if(filename.endsWith(".db"))
            {
                dbStemname = filename.substring(0,filename.length()-3);
                filename =  DB_GROUP + filename;
            }
            else
            {
                dbStemname = filename.substring(3, filename.length());
            }
            String directory = chooser.getDirectoryPath();
            filePathname = directory + File.separator + filename;

            File f = new File(filePathname);
            if( f.exists())	{
				boolean ok = StsYesNoDialog.questionValue(model.win3d, "Database exists, really overwrite it ?");
		        if (ok)
					f.delete();
				else
					return null;
			}


//            project.setDatabase(filePathname);
            return filePathname;
        }
        catch(Exception e)
        {
            StsException.outputException("StsNewModel.selectDatabasePathname() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    public boolean end()
    {
        logMessage("Model saved to:  " + "db." + model.getName());
        StsMessageFiles.writeTimeStampToFile(StsMessageFiles.LOG);
        return saveSuccess;
    }

    /* build a model name from the db filename */
    private String buildName(String filename)
    {
        String name = null;
        StringTokenizer stok = new StringTokenizer(filename, ".");
        int nTokens = stok.countTokens();
        for (int i=0; i<nTokens-1; i++) stok.nextToken();
        name = new String(stok.nextToken());
        if (name==null) name = new String("Default"); // give it something
        return name;
    }
/*
    class EmptyDB implements Runnable
    {
        StsDBFile db;
        EmptyDB(StsDBFile db) { this.db = db; }
        public void run() { db.empty(); }
    }
*/
}


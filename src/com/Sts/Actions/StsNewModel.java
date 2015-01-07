
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.io.*;

public class StsNewModel extends StsAction
{
    boolean success = false;
    String projectDirname = null;
    String dbStemname;

    transient static public final String title = "Select project file";
    transient static public final String filterString = "Project Files";
    transient static public final String filter = "proj";

    static public final String DB_GROUP = "db.";
    static public final String jarFilename = "projects.jar";

    public StsNewModel(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public StsNewModel(StsActionManager actionManager, String projectDirname)
    {
        super(actionManager, true);
        this.projectDirname = projectDirname;
    }

	public boolean start()
    {
        final StsModel oldModel = model;
        StsWin3d newWin3d = null;
		StsGLPanel newglPanel;
        StsProjectParms parms = null;
        File file = null;

        try
        {
            oldModel.disableDisplay();

            // because we are changing models, we need to end this transaction,
            // close the db, change models in actionManager and start a new db and transaction
            actionManager.endCurrentAction();
            oldModel.commit();
            oldModel.stopTime();
            oldModel.close();
			String dbPathname = null;
			StsProject project = null;
			do
            {
				// if user hasn't selected a preexisting project directory, select one here
				if(projectDirname == null)
				{
					projectDirname = selectProjectDirectory();
					if(projectDirname == null)return false;
				}

				project = new StsProject(projectDirname);
                String modelFolderPathname = project.getModelFullDirString();
                StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(model.win3d, "Specify Database name ", modelFolderPathname, DB_GROUP);
                int option = chooser.showReturnOption();
                if(option != JFileChooser.APPROVE_OPTION)
                    return false;
                dbPathname = selectDatabasePathname(chooser);
			}
            while (dbPathname == null);

            if(dbPathname == null) return false;

            if(!setupProject(project)) return false;

            model = StsModel.constructor(project, dbPathname);
            logMessage("Model now switched from:  " + oldModel.getName() + " to " + model.getName());
            success = true;
            return success;
        }
        catch (Exception e)
        {
			StsException.outputException("StsNewModel.start() failed.", e, StsException.WARNING);
            if(newWin3d != null) logMessage("Unable to classInitialize display.");
            return success;
        }
    }

    String selectProjectDirectory()
    {
        String currentDirectory = model.getProject().getRootDirString();
        JFileChooser chooser = new JFileChooser(currentDirectory);
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setDialogTitle("Select or Enter Desired Directory and Press Open Directory");
        chooser.setApproveButtonText("Open Directory");
        chooser.showOpenDialog(null);
        File selectedFile = chooser.getSelectedFile();
        if(selectedFile == null)
            return null;
        String projectDirname = selectedFile.getPath();
        if(projectDirname != null) StsProject.addProjectDirectory(projectDirname);
        return projectDirname;
    }

    public String selectDatabasePathname(StsFileChooser chooser)
    {
        try
        {
            String filePathname = chooser.getFilePath();
            if(filePathname == null) return null;

            String filename = chooser.getFilename();
            if(!filename.startsWith(DB_GROUP))
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
            if( f.exists())
            {
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

    private boolean setupProject(StsProject project)
    {
        try
        {
            StsNewProjectDialog propertyDialog = new StsNewProjectDialog(project);
            propertyDialog.setVisible(true);
            if(propertyDialog.wasCanceled())
                return false;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsNewModel.setupProject() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }


    public boolean end()
    {
        logMessage(" ");
        if(success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return success;
    }
}

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Actions.*;
import com.Sts.Collaboration.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsCreateCollaboration extends StsAction // implements Runnable
{
	StsCollaboration collaboration = null;
    private boolean openSuccess = false;

    public StsCreateCollaboration(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public boolean start()
//    public void run()
    {
//		collaboration = StsCollaboration.getInstance(model, entry);
//		SessionManager sessionManager = new SessionManager();
//		if(!sessionManager.createSession()) return false;

//	    if(sessionManager.isLeader()) return true;
        final StsModel oldModel = model;
        String filename = null;
        StsDBFile db = null;
        String directory;

        try
        {
            FilenameFilter filenameFilter = new DBFilenameFilter();

            directory = "." + File.separator;
            if (oldModel != null)
            {
                oldModel.commit();
                db = oldModel.getDatabase();
                if (db != null)
                {
                    db.close();
                    db.getURLDirectory();
                    //                        directory = oldModel.getProject().getModelFullDirString();
                    //                    StsModel.setCurrentModel(null);
                }
            }
			BufferedInputStream bis = new BufferedInputStream(collaboration.pc.getInputStream());
			DataInputStream dis = new DataInputStream(bis);
            StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { oldModel.win3d.closeWindows(); }});

            logMessage("Model now switched from:  " + oldModel.getName() + " to collaboration db from leader.");

	        StsModel newModel = StsModel.constructor(dis);

            openSuccess = true;
            actionManager.endCurrentAction();
//            newModel.win3d.glPanel3d.getActionManager().fireChangeEvent();
            newModel.initializeActionStatus();
			newModel.refreshObjectPanel();
            System.out.println("Peer reading db from leader");

            return openSuccess;
        }
        catch (Exception e)
        {
            StsException.outputException("StsCreateCollaboration.start() failed.", e, StsException.WARNING);
            return false;
        }
        finally
        {
            statusArea.textOnly();
        }
    }

    /**
     * Method generateKernelConfigContainer.
     * @return KernelConfigContainer
     */
    /*
         private KernelConfigContainer generateKernelConfigContainer() {

        KernelConfigContainer kcc = new KernelConfigContainer();

        // add a UserID ---------------
        kcc.add(new UserConfigContainer());

        // create the nengine config ------------ just use defaults.
        NetworkEngineConfigContainer necc = new NetworkEngineConfigContainer();
        kcc.add(necc);

        // add the BasicConnectivityManager configuration to the config
        BCMConfigContainer bcm = new BCMConfigContainer();

        bcm.setID(NetworkEngineConfigContainer.CM_ID);
        bcm.setGuaranteedPort(4600);
        bcm.setUDPPort(4000);

        necc.setConnectivityManagerClassName(
            "ar.connmanager.basic.BasicConnectivityManager");
        necc.add(bcm);

        // create the console configs ------------
        ConsoleConfigContainer ccc = new ConsoleConfigContainer();

        // create the primary console
        ccc.setConsole(new TextConsole());
        ccc.setID(KernelConfigContainer.PRIMARY_CONSOLE_ID);

        kcc.add(ccc);

        // ensure everything is configured properly
        try {
            kcc.validate();
        } catch (Exception ve) {
            System.err.println("Invalid configuration");
            ve.printStackTrace();
            return null;
        }

        return kcc;
         }
     */

    private StsDBFile openExistingDB(String dirPath)
    {
        String filename = "";

        try
        {
            StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(model.win3d, "Open model file", dirPath, "db.");

            while (true)
            {
                if (!chooser.show())
                {
                    actionManager.endCurrentAction();
                    return null;
                }
                else
                {
                    File f = new File(chooser.getFilePath());
                    if (!f.exists())
                    {
                        logMessage("File " + f.getPath() + " not found.");
                    }
                    else if (f.isDirectory())
                    {
                        logMessage("Selected file " + f.getPath() + " is a directory.");
                    }
                    else
                    {
                        try
                        {
                            filename = chooser.getFilename();
                            dirPath = chooser.getDirectoryPath();
                            StsFile file = StsFile.constructor(dirPath, filename);
                            return StsDBFileModel.openRead(file, null, true);
                        }
                        catch (Exception e)
                        {
                            new StsMessage(model.win3d, StsMessage.WARNING, "Couldn't open db file: " + filename);
                            return null;
                        }
                    }
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsOpenExistingDB() failed.", e, StsException.WARNING);
            return null;
        }
    }

    private StsDBFile openNewDB(String dirPath)
    {
        String filename = "null";

        try
        {
            StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(model.win3d, "New Database (use db. as prefix)", dirPath, "db.");
            if (!chooser.showSave())
            {
                return null;
            }

            String filePath = chooser.getFilePath();
            File f = new File(filePath);
            if (f.exists())
            {
                f.delete();
            }
            try
            {
                filename = chooser.getFilename();
                StsFile file = StsFile.constructor(dirPath, filename);
                return StsDBFileModel.openWrite(file, null);
            }
            catch (Exception e)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Couldn't open new db file: " + filename);
                return null;
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsOpenNewDB() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static final class DBFilenameFilter implements FilenameFilter
    {
        public DBFilenameFilter()
        {
        }

        public boolean accept(File dir, String name)
        {
            System.out.println("StsOpenModel.DBFilenameFilter.accept():" +
                               " dirNo: " + dir.toString() + " filename: " + name);
            return name.startsWith("db.");
        }
    }

    public boolean end()
    {
        if (openSuccess)
        {
            logMessage("Model opened successfully.");
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        else
        {
            logMessage("Model open failed.");
        }
        return openSuccess;
    }
}

//class PawnWriter implements TransactionListener {
//
//    private Pawn pawn;

// constructor
//    PawnWriter(Pawn pawn) {
//        this.pawn = pawn;
//    }

/**
 * @see com.Sts.DB.TransactionListener#transactionOccured(byte[])
 */
//    public void transactionOccured(byte[] transaction) {

//        System.out.println("************ Transaction occured.");
// apply the transaction to the pawn
//        this.pawn.executeCommand(new ApplyTransactionCmd(transaction));
//    }
//}

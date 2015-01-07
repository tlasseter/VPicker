package com.Sts.MVC;
/**
 * <p>Title:        StsMessageFile</p>
 * <p>Description:  Class used to output messages to the message panel and to log files</p>
 * <p>Copyright:    Copyright (c) 2001</p>
 * <p>Company:      4D Systems LLC</p>
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class StsMessageFiles
{
    /** Default file folder to write logs and message files */
	static public final String MESSAGE_FOLDER = "messageFiles";
	static File messageFolder = null;

	static JTabbedPane tabbedTextPanels;
	static StsRandomAccessFileListModel[] listModels;

    /** current message string */
    static String message = null;
    /** current message type */
    static int type;

    /** Static LOG constant = 0 */
	static public final int LOG = 0;
    /** Static INFO constant = 1 */
	static public final int INFO = 1;
    /** Static ERROR constant = 2 */
	static public final int ERROR = 2;
    /** Static JOURNAL constant = 3 */
	static public final int JOURNAL = 3;
    /** Available file prefixes - logFile, infoFile, errorFile, and journalFile */
	static public final String[] filePrefixes = new String[] { "logFile", "infoFile", "errorFile", "journalFile" };
    /** Available message labels - Log, Info, Error and Journal */
    static public final String[] labels = new String[] { "Log", "Info", "Error", "Journal" };
    /** Number of message types */
	static public final int nTypes = 3;
    /** Maximum number of lines per block */
	static public final int maxLinesPerBlock = 200;

    /**
     * Create the message tab panel
     * @params _tabbedTextPanels the tab pane for messages
     **/
    static public void create(JTabbedPane _tabbedTextPanels)
    {
        tabbedTextPanels = _tabbedTextPanels;
        for(int n = 0; n < nTypes; n++)
        {
            StsListScrollPanel scrollPanel = new StsListScrollPanel(null, n);
            tabbedTextPanels.add(labels[n], scrollPanel);
        }
    }

    /**
     * Initialize the tabbed pane
     * @param model the current model
     * @param _tabbedTextPanels the tabbed Panel for messages
     * @param messageFolderDirectory the folder directory where message files are written
     */
    static public void initialize(StsModel model, JTabbedPane _tabbedTextPanels, String messageFolderDirectory)
    {
        tabbedTextPanels = _tabbedTextPanels;

        if(listModels != null)
        {
            for(int n = 0; n < nTypes; n++)
                if(listModels[n] != null) listModels[n].closeFile();
        }

        messageFolder = new File(messageFolderDirectory + File.separator + MESSAGE_FOLDER);
        if(!messageFolder.exists()) messageFolder.mkdir();

        listModels = new StsRandomAccessFileListModel[nTypes];
        String name = model.getName();
        boolean deleteExistingFile = true;

        for(int n = 0; n < 3; n++)
        {
            listModels[n] = getCreateListModel(messageFolder, filePrefixes[n] + "." + name + ".txt", deleteExistingFile, n);
            StsListScrollPanel scrollPanel = new StsListScrollPanel(listModels[n], n);
            tabbedTextPanels.add(labels[n], scrollPanel);
        }
    }

    /**
     * Create the message tab panel
     * @param tabbedTextPanels the tabbed pane where the panel will be added
     * @param label the label for the new panel
     * @param type the type of panel
     */
    private void createMessagePanel(JTabbedPane tabbedTextPanels, String label, int type)
    {
        StsListScrollPanel scrollPanel = new StsListScrollPanel(null, type);
        tabbedTextPanels.add(label, scrollPanel);
	}

    /**
     * classInitialize the Tabbed Pane
     * @param model the current model
     * @param _tabbedTextPanels the tabbed Panel for messages
     */
    static public void initialize(StsModel model, JTabbedPane _tabbedTextPanels)
    {
        String rootDirString = model.getProject().getRootDirString();
        initialize(model, _tabbedTextPanels, rootDirString);
    }

	static private StsRandomAccessFileListModel getCreateListModel(File directory, String name, boolean deleteExistingFile, int type)
	{
		try
		{
			File file = new File(directory, name);
			if(file.exists()) file.delete();
			file = new File(directory, name);
			// hack!  Still manages to read old files even when we want to delete them!
		/*
			if(file.exists())
			{
				if(deleteExistingFile)
				{
					file.delete();
					file = new File(directory, name);
				}
			}
//			else
//				file.createNewFile();
        */

			RandomAccessFile raf = new RandomAccessFile(file, "rw");
			return new StsRandomAccessFileListModel(raf, maxLinesPerBlock, type);
		}
		catch(Exception e)
		{
			StsException.outputException("StsMessageFiles.getCreateMessageFile() for file: " + name + " failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    /**
     * Get the list model associated with one of the three message types
     * @param type the desired list model
     * @return the random access list model associated with a particular message type
     */
	static public StsRandomAccessFileListModel getListModel(int type)
	{
	    if(type < 0 || type >= nTypes) return null;
		return listModels[type];
	}

    /**
     * Write a time stamp to all the log files.
     */
    static public void writeTimeStampToFiles()
    {
		String msg =  DateFormat.getDateTimeInstance(DateFormat.LONG,
							DateFormat.LONG).format(new Date());
		for(int n = 0; n < nTypes; n++) message(msg, n);
    }

    /**
     * Write a time stamp to the specified log file
     * @param type the log file type to time stamp
     */
    static public void writeTimeStampToFile(int type)
    {
		String msg =  DateFormat.getDateTimeInstance(DateFormat.LONG,
							DateFormat.LONG).format(new Date());
		message(msg, type);
    }

    /**
     * Write out a log message to the log message file
          * @param msg the message
          */
         static public boolean logMessage(String msg)
         {
             message = (new SimpleDateFormat("dd-MMM-yy kk:mm:ss.S")).format(new Date(System.currentTimeMillis())) + ":" + msg;
             type = LOG;
             StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { doLogMessage(); }});
             return true;
         }

         static public boolean doLogMessage()
         {
             return message(message, type);
         }

         /**
          * Output a message to the message panel if it exists, if not print to the system
          * @param msg the message
          * @param type the type (log, info, error or journal)
          */
         static private boolean message(String msg, int type)
         {
             if(tabbedTextPanels == null)
             {
                 return false;
             }

             int nTabs = tabbedTextPanels.getTabCount();
             if(type > nTabs-1)
             {
                 return false;
             }

             if (SwingUtilities.isEventDispatchThread())
             {
                 tabbedTextPanels.setSelectedIndex(type);
             }
             else
             {
                 try
                 {
                     final int _type = type;
                     SwingUtilities.invokeAndWait
                     (
                         new Runnable()
                         {
                             public void run()
                             {
                                 tabbedTextPanels.setSelectedIndex(_type);
                             }
                         }
                     );
                 }
                 catch (Exception e)
                 {
                     e.printStackTrace();
                     return false;
                 }
             }
             if(listModels == null || listModels[type] == null) return true;
             listModels[type].writeLine(msg);
             return true;
         }

         /**
          * Write out an informational message to the info message file
          * @param msg the message
          */
         static public boolean infoMessage(String msg)
         {
             msg = (new SimpleDateFormat("dd-MMM-yy kk:mm:ss.S")).format(new Date(System.currentTimeMillis())) + ":" + msg;
             return message(msg, INFO);
         }

         /**
          * Write out a error message to the error message file
          * Return true if written to System.out.
          * @param msg the message
          */
         static public boolean errorMessage(String msg)
         {
             msg = (new SimpleDateFormat("dd-MMM-yy kk:mm:ss.S")).format(new Date(System.currentTimeMillis())) + ":" + msg;
             return message(msg, ERROR);
         }

         static public void errorStackTrace(Exception e)
         {
             message(e.toString(), ERROR);
             StackTraceElement[] elements = e.getStackTrace();
             for(StackTraceElement element : elements)
                message("\t" + element.toString(), ERROR);
         }
     }

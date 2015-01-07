
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.io.*;

public class StsFileChooser implements ListSelectionListener
{
    private Frame parent = null;
    private FileFilter[] filters = null;

    private int[] currentSelected = null;
    private int[] lastSelected = null;
    private JList fileList = null;
    private boolean multipleSelection = false;
    private int fileSelectionMode = JFileChooser.FILES_AND_DIRECTORIES;
    String title = "test";
    String path;
    JFileChooser chooser = null;

    /** constructor for multiple selection file with prefix filter */
    private StsFileChooser(Frame parent, String title, String path, String filter)
    {
        this(parent, title, path, filter, true, true);
    }

    private StsFileChooser(Frame parent, String title, String path, String filter, boolean isPrefix, boolean multipleSelection)
    {
        this.parent = parent;
        if(filter != null)
            filters = new StsFileFilter[]{new StsFileFilter(filter, isPrefix)};
        this.title = title;
        this.path = path;
        this.multipleSelection = multipleSelection;
        createChooser(path);
    }

    public StsFileChooser(Frame parent, String title, String path)
    {
        this(parent, title, path, ".*", false, false);
        fileSelectionMode = JFileChooser.DIRECTORIES_ONLY;
	}

    static public StsFileChooser createFileChooser(Frame parent, String title, String path)
    {
        return new StsFileChooser(parent, title, path, null, true, false);
    }

    static public StsFileChooser createFileChooser(Frame parent, String title, String path, String filter, boolean isPrefix, boolean multipleSelection)
    {
        return new StsFileChooser(parent, title, path, filter, isPrefix, multipleSelection);
    }

    static public StsFileChooser createFileChooserPrefix(Frame parent, String title, String path, String filter)
    {
        return new StsFileChooser(parent, title, path, filter, true, false);
    }

    static public StsFileChooser createFileChooserPostfix(Frame parent, String title, String path, String filter)
    {
        return new StsFileChooser(parent, title, path, filter, false, false);
    }

    static public StsFileChooser createMultiFileChooserPrefix(Frame parent, String title, String path, String filter)
    {
        return new StsFileChooser(parent, title, path, filter, true, true);
    }

    static public StsFileChooser createMultiFileChooserPostfix(Frame parent, String title, String path, String filter)
    {
        return new StsFileChooser(parent, title, path, filter, false, true);
    }

    static public StsFileChooser createDirectoryChooser(Frame parent, String title, String path)
    {
        StsFileChooser chooser = new StsFileChooser(parent, title, path, ".*", false, false);
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        return chooser;
    }

    static public StsFileChooser createMultiDirectoryChooser(Frame parent, String title, String path)
    {
        StsFileChooser chooser = createDirectoryChooser(parent, title, path);
        chooser.setMultiSelectionEnabled(true);
        return chooser;
    }

    protected void setMultiSelectionEnabled(boolean multiSelect)
    {
        chooser.setMultiSelectionEnabled(multiSelect);
    }

    // accessors

    private void createChooser()
    {
        String dirPath = "." + File.separator;
        createChooser(dirPath);
    }

    private void createChooser(String path)
    {
		if(chooser == null) chooser = new JFileChooser();
        chooser.setCurrentDirectory(new File(path));
//        chooser = new JFileChooser(path, new WindowsAltFileSystemView());
        chooser.setDialogTitle(title);

        // set the custom file type
        if (filters!=null)
        {
            for( int i=0; i<filters.length; i++ )
            	chooser.addChoosableFileFilter(filters[i]);
        	chooser.setFileFilter(filters[0]);
        }

        // set selection mode
       	chooser.setMultiSelectionEnabled(multipleSelection);
        chooser.setFileSelectionMode(fileSelectionMode);

        // add listener to get around bug in multiple selection
	/*
       	JList list = getFileList();
        if( list != null )
        {
        	list.addListSelectionListener(this);
            if( !multipleSelection )
	            list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        }
   */
    }

    public JFileChooser getFileChooser()
    {
        if(chooser == null) createChooser();
        return chooser;
    }

    public void setFileSelectionMode(int mode)
    {
        if(chooser == null) createChooser();
        chooser.setFileSelectionMode(mode);
    }

    public int getFileFilterIndex()
    {
        FileFilter currentFilter = chooser.getFileFilter();

        for( int i=0; i<filters.length; i++ )
        {
            if( currentFilter == filters[i] ) return i+1;
        }
        return 0;
    }

    /** pop up the selector and see if anything was selected */
    public boolean show()
    {
        return chooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION;
    }

    public int showReturnOption()
    {
        return chooser.showOpenDialog(parent);
    }

    public boolean showSave()
    {
        return chooser.showSaveDialog(parent) == JFileChooser.APPROVE_OPTION;
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	JList list = (JList) e.getSource();
        lastSelected = currentSelected;
		currentSelected = list.getSelectedIndices();
/*
		int numSelected = 0;
        if( lastSelected != null ) numSelected = lastSelected.length;
        System.out.println("last selected: " + numSelected);
		numSelected = 0;
        if( currentSelected != null ) numSelected = currentSelected.length;
        System.out.println("selected: " + numSelected);
*/
    }

    public String getFilePath()
    {
        waitForLoadFilesThreadToComplete();
        File selectedFile = chooser.getSelectedFile();
        if(selectedFile == null) return null;
        return selectedFile.getPath();
    }

    /** return an array of filenames (without paths) */
    public String[] getFilenames()
    {
        waitForLoadFilesThreadToComplete();
        File[] selectedFiles = chooser.getSelectedFiles();
        if(selectedFiles == null) return null;
        int nFilenames = selectedFiles.length;
        String[] filenames = new String[nFilenames];
        for(int n = 0; n < nFilenames; n++)
            filenames[n] = selectedFiles[n].getName();
        return filenames;
    }

    public String[] getPathnames()
    {
        waitForLoadFilesThreadToComplete();
        File[] selectedFiles = chooser.getSelectedFiles();
        if(selectedFiles == null) return null;
        int nFilenames = selectedFiles.length;
        String[] pathnames = new String[nFilenames];
        for(int n = 0; n < nFilenames; n++)
            pathnames[n] = selectedFiles[n].getAbsolutePath();
        return pathnames;
    }

    /** Workaround for deadlock in JFileChooser, possibly bug 6713352 */
    private static void waitForLoadFilesThreadToComplete()
    {
        Thread[] threads = new Thread[Thread.activeCount()];
        Thread.enumerate(threads);
        for (Thread t : threads)
        {
            if (t != null && t.getClass().getName().equals("javax.swing.plaf.basic.BasicDirectoryModel$LoadFilesThread"))
            {
                try
                {
                    t.join();
                }
                catch (InterruptedException e)
                {
                    // ignore
                }
            }
        }
    }


    public void cancelSelection()
    {
        ;
    }
    public String getFilename()
    {
        waitForLoadFilesThreadToComplete();
        File selectedFile = chooser.getSelectedFile();
        if(selectedFile == null) return null;
        return selectedFile.getName();
    /*
        String[] filenames = getFilenames();
        if(filenames == null || filenames.length != 1) return null;
        return filenames[0];
    */
    }

    /** return an array of filenames (without paths) */
/*
    public String[] getFilenames()
    {
        File[] files = chooser.getSelectedFiles();
        if( files == null ) files = new File[] { chooser.getSelectedFile() };

        if( files == null || files.length == 0 ) return null;

        String[] filenames = new String[files.length];
        for( int i=0; i<files.length; i++ ) filenames[i] = files[i].getName();
        return filenames;
    }
*/
    public File[] getFiles()
    {
        waitForLoadFilesThreadToComplete();
        if(chooser.isMultiSelectionEnabled())
        {
            return this.chooser.getSelectedFiles();
        }
        else
        {
            File file = chooser.getSelectedFile();
            return new File[] { file };
        }
    }

    public String getDirectoryPath()
    {
        return chooser.getCurrentDirectory().getPath();
    }


    public static void main(String[] args)
    {
    	try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
 //           StsFileChooser fc = StsFileChooser.createMultiDirectoryChooser(null, "test", ".");
            // System.out.println("you selected "+fc.getPathnames()[0]);
            String dirPath = "." + File.separator;
            StsFileChooser fc = StsFileChooser.createFileChooserPrefix(null, "Help", dirPath, "db.");
            // StsFileChooser fc = new StsFileChooser(null, "Help", dirPath, "db.", true, false);
            fc.show();
            System.out.println("you selected "+fc.getFilename());
        }
        catch(Exception e) { e.printStackTrace(); }
    }
}


















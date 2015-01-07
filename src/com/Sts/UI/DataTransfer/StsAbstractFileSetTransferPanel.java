package com.Sts.UI.DataTransfer;

import com.Sts.IO.*;
import com.Sts.DBTypes.*;

import javax.swing.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 14, 2010
 * Time: 8:26:08 AM
 * To change this template use File | Settings | File Templates.
 */

/** An abstractFileset is either a directory containing a set of files which satisfy a filenameFilter or a set of files in an StsJar file
 *  or in an StsWebStartJar file. This transfer panel allows the user to select the files s/he wishes to load.
 */
public class StsAbstractFileSetTransferPanel extends StsDirectoryObjectsTransferPanel
{
    StsAbstractFileSet fileset;

    public StsAbstractFileSetTransferPanel(String title, String currentDirectory, StsObjectTransferListener listener, int width, int height, StsAbstractFileSet fileset)
    {
        this.fileset = fileset;
        super.initialize(title, currentDirectory, listener, width, height);
    }

    public StsAbstractFile[] getAvailableObjects()
    {
        return fileset.getFiles();
    }

    public void setCurrentDirectory(String directory)
    {
        currentDirectory = directory;
        if(fileset == null) return;
        fileset.setDirectory(directory);
        initializeSetAvailableObjects();
    }

    public StsAbstractFile[] initializeAvailableObjects()
    {
        if(fileset == null) return null;
        return fileset.initializeAvailableFiles();
    }
    
    public static void main(String[] args)
    {
        TestAbstractFileSetTransferPanelListener listener = new TestAbstractFileSetTransferPanelListener();
        String currentDirectory = "H:\\FractureAnalysis\\Surfaces";
        StsFilenameFilter filenameFilter = new StsFilenameFilter(StsSurface.seismicGrp, "txt");
        StsFileSet fileset = StsFileSet.constructor(currentDirectory, filenameFilter);
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        }
        catch (Exception e) { }
        TestAbstractFileSetTransferPanel panel = new TestAbstractFileSetTransferPanel(currentDirectory, listener, 400, 100, fileset);
        com.Sts.Utilities.StsToolkit.createDialog(panel);
    }
}

class TestAbstractFileSetTransferPanel extends StsAbstractFileSetTransferPanel
{
    TestAbstractFileSetTransferPanel(String currentDirectory, TestAbstractFileSetTransferPanelListener listener, int width, int height, StsAbstractFileSet fileset)
    {
        super("Test Directory Objects Panel", currentDirectory, listener, width, height, fileset);
    }
}

class TestAbstractFileSetTransferPanelListener implements StsObjectTransferListener
{

    TestAbstractFileSetTransferPanelListener()
    {
    }

    public void addObjects(Object[] objects)
    {
        System.out.println("add Objects:");
        printObjects(objects);
    }

    private void printObjects(Object[] objects)
    {
        for (int n = 0; n < objects.length; n++)
            System.out.println("    " + objects[n].toString());
    }

    public void removeObjects(Object[] objects)
    {
        System.out.println("remove Objects:");
        printObjects(objects);
    }
    
    public void objectSelected(Object selectedObject)
    {
        System.out.println("selected Object:" + selectedObject.toString());
    }

    public boolean addAvailableObjectOk(Object object) { return true; }
}


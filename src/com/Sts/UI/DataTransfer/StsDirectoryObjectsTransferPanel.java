package com.Sts.UI.DataTransfer;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version c62e
 */


/**
 * This abstract panel provides a browse button for data directories.
 * The user uses the add, add-all to move files from the left list to the right,  and remove, remove-all buttons to move objects
 * from the right list to the left.
 * The StsObjectTransferListener (interface) can be provided which listens to events as objects are added and removed
 * from the right-hand side list.
 */
abstract public class StsDirectoryObjectsTransferPanel extends StsObjectsTransferPanel
{
    protected String currentDirectory;
    protected StsDirectoryBrowseGroupBox directoryGroupBox;
    /**
     * This method is called when the browseDatastore button is pushed.
     * The implementation should do whatever is necessary to find the datastore the user is looking for.
     * It could be a fileBrower, or databaseConnector and might involve bringing up a dialog box or whatever.
     *
     * @return datstore the user selects
     */

    public StsDirectoryObjectsTransferPanel()
    {
    }

    public void initialize(String title, String currentDirectory, StsObjectTransferListener listener, int width, int height)
    {
        this.currentDirectory = currentDirectory;
        super.initialize(title, listener, width, height);
    }

    public void constructTransferPanel()
    {
        directoryGroupBox = new StsDirectoryBrowseGroupBox(this, "currentDirectory");
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        addEndRow(directoryGroupBox);
        gbc.weighty = 1.0;
        super.constructTransferPanel();
    }

    public void initialize()
    {
        setCurrentDirectory(currentDirectory);
    }

    public int getSelectedIndex() { return selectedList.getSelectedIndex(); }

    public String getDirectoryName() { return getCurrentDirectory(); }

    public String getCurrentDirectory()
    {
        return currentDirectory;
    }

    public void setCurrentDirectory(String currentDirectory)
    {
        this.currentDirectory = currentDirectory;
        initializeSetAvailableObjects();
    }

    public void initializeSetAvailableObjects()
    {
        Object[] availableObjects = initializeAvailableObjects();
        setAvailableObjects(availableObjects);    
    }

    public static void main(String[] args)
    {
        TestDirectoryObjectsTransferPanelListener listener = new TestDirectoryObjectsTransferPanelListener();
        String currentDirectory = "c:\\data\\FractureAnalysis\\Surfaces";
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        }
        catch (Exception e) { }
        TestDirectoryObjectsTransferPanel panel = new TestDirectoryObjectsTransferPanel(currentDirectory, listener, 400, 100);
        com.Sts.Utilities.StsToolkit.createDialog(panel);
    }
}

class TestDirectoryObjectsTransferPanel extends StsDirectoryObjectsTransferPanel
{
    TestDirectoryObjectsTransferPanel(String currentDirectory, TestDirectoryObjectsTransferPanelListener listener, int width, int height)
    {
        super.initialize("Test Directory Objects Panel", currentDirectory, listener, width, height);
    }

    public Object[] getAvailableObjects()
    {
        return new String[] { "file1", "file2", "file3"};
    }

    public Object[] initializeAvailableObjects()
    {
       return new String[] { "file1", "file2", "file3"};
    }
}

class TestDirectoryObjectsTransferPanelListener implements StsObjectTransferListener
{

    TestDirectoryObjectsTransferPanelListener()
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

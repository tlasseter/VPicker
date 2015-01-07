package com.Sts.UI.DataTransfer;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.DBTypes.*;

import javax.swing.*;
import javax.swing.event.*;
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
 * This abstract panel provides a browse button for datastores; once selected a set of available projects are isVisible
 * in a comboBox;  once a project has been selected, a list of available objects appear in a left-hand list.  The user
 * uses the add, add-all, remove, and remove-all buttons to add or remove objects from a list on the right-hand side.
 * A datastore implements the StsDatastoreFace interface.  A listener must be provided which implements the StsObjectTransferFace
 * interface.  This listener responds to events as objects are added and removed from the right-hand side list.
 * The datastore implementation should be subclassed in such a way that it provides objects of a particular type (wells,
 * surfaces, prestack 3d datasets, etc) and properly implements the getProjectObjects() method for the particular type of
 * object desired.
 */
abstract public class StsDatastoreObjectTransferPanel extends StsObjectsTransferPanel
{
    /** datastore which contains projects which contain objects of desired type */
    protected StsDatastoreFace datastore;
    private String datastoreName = "none";
    private String projectName = "none";

    protected StsJPanel directoryGroupBox = new StsJPanel();
    private StsButton dataStoreBrowseButton = new StsButton("Get Data Store", "Browse for data store.", this, "browseDatastore");
    private StsStringFieldBean datastoreNameBean = new StsStringFieldBean(this, "datastoreName", null);
    private StsComboBoxFieldBean projectListBean = new StsComboBoxFieldBean(this, "projectName", "Project:", new String[]{"Please select project"});

    /**
     * This method is called when the browseDatastore button is pushed.
     * The implementation should do whatever is necessary to find the datastore the user is looking for.
     * It could be a fileBrower, or databaseConnector and might involve bringing up a dialog box or whatever.
     *
     * @return datstore the user selects
     */
    abstract public StsDatastoreFace getDatastore();

    public StsDatastoreObjectTransferPanel(/*StsDatastoreFace datastore*/ StsDatastoreFace datastore, StsObjectTransferListener listener, int width, int height)
    {
        initialize("Data Set Transfer Selection", listener, width, height);
        this.datastore = datastore;
        if (datastore != null) initializeDatastoreName();

        directoryGroupBox.gbc.fill = GridBagConstraints.NONE;
        directoryGroupBox.gbc.anchor = GridBagConstraints.WEST;
        directoryGroupBox.gbc.weightx = 0.0;
        //          directoryGroupBox.gbc.weighty = 0.0;
        directoryGroupBox.addToRow(dataStoreBrowseButton);
        directoryGroupBox.gbc.weightx = 1.0;
        directoryGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        directoryGroupBox.gbc.anchor = GridBagConstraints.EAST;
        datastoreNameBean.setColumns(30);
        directoryGroupBox.addEndRow(datastoreNameBean);
        directoryGroupBox.gbc.gridx = 1;
        directoryGroupBox.gbc.anchor = GridBagConstraints.WEST;
        directoryGroupBox.add(projectListBean);

        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weighty = 0.0;
        addEndRow(directoryGroupBox);
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty = 1.0;
        addEndRow(selectObjectsPanel);

        availableList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        availableList.addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent e)
                    {
                        availableObjectSelected(e);
                    }
                }
            );
        selectedList.addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent e)
                    {
                        selectedObjectSelected(e);
                    }
                }
            );
        setProjectObjects();
    }

    public void initializeDatastoreName()
    {
        if (datastore == null) return;
        String datastoreName = datastore.getDatastoreName();
        setDatastoreName(datastoreName);
        datastoreNameBean.setValue(datastoreName);
        this.revalidate();
        initializeProjectName();
    }

    public void initializeProjectName()
    {
        if (datastore == null) return;
        String[] projectNames = datastore.getProjectNames();
        if (projectNames == null || projectNames.length == 0) return;
        setProjectName(projectNames[0]);
        projectListBean.setListItems(projectNames);
        projectListBean.doSetValueObject(projectName);
    }

    public void browseDatastore()
    {
        datastore = getDatastore();
        if (datastore == null) return;
        initializeDatastoreName();
    }

    public void setDatastoreName(String name) { datastoreName = name; }

    public String getDatastoreName() { return datastoreName; }

    public void selectedObjectSelected(ListSelectionEvent e)
    {
        if (e.getValueIsAdjusting()) return;
        Object source = e.getSource();
        if (!(source instanceof JList)) return;
        Object selectedObject = selectedList.getSelectedValue();
        if (selectedObject == null) return;
        if (listener != null) listener.objectSelected(selectedObject);
    }

    public void availableObjectSelected(ListSelectionEvent e)
    {
        if (e.getValueIsAdjusting()) return;
        Object source = e.getSource();
        if (!(source instanceof JList)) return;
        Object selectedObject = availableList.getSelectedValue();
        if (selectedObject == null) return;
        if (listener != null) listener.objectSelected(selectedObject);
    }

    public void setProjectName(String projectName)
    {
        this.projectName = projectName;
        datastore.setProjectName(projectName);
        setProjectObjects();
    }

    public String getProjectName() { return projectName; }

    public void setProjectObjects()
    {
        availableListModel.clear();
        if (datastore == null)
            return;    // added by lkw to avoid a null pointer exception on the next line
        availableObjects = datastore.getObjects();
        if (availableObjects == null) return;
        for (int n = 0; n < availableObjects.length; n++)
            availableListModel.addElement(availableObjects[n]);
    }

    public static void main(String[] args)
    {
        TestDatastoreObjectTransferPanelListener listener = new TestDatastoreObjectTransferPanelListener();
        TestObjectTransferPanelDatastore datastore = new TestObjectTransferPanelDatastore();
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        }
        catch (Exception e) { }
        TestDatastoreObjectTransferPanel panel = new TestDatastoreObjectTransferPanel(datastore, listener, 400, 100);
        com.Sts.Utilities.StsToolkit.createDialog(panel);
    }
}

class TestDatastoreObjectTransferPanel extends StsDatastoreObjectTransferPanel
 {
     TestDatastoreObjectTransferPanel(TestObjectTransferPanelDatastore datastore, StsObjectTransferListener listener, int width, int height)
     {
         super(datastore, listener, 400, 200);
     }

     public Object[] getAvailableObjects()
     {
         return datastore.getObjects();
     }

     public Object[] initializeAvailableObjects()
     {
         datastore.initialize();
         return datastore.getObjects();
     }

     public StsDatastoreFace getDatastore()
     {
         return new TestObjectTransferPanelDatastore();
     }
 }

class TestDatastoreObjectTransferPanelListener implements StsObjectTransferListener
 {

     TestDatastoreObjectTransferPanelListener()
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

class TestObjectTransferPanelDatastore implements StsDatastoreFace //implements StsDatastoreFace
{
    String name = "Test datastore";
    String projectName;
    String[] projectNames;
    String[][] objectNames;

    TestObjectTransferPanelDatastore()
    {
        createProjects();
        createProjectObjects();
    }

    private void createProjects()
    {
        projectNames = new String[4];
        for (int n = 0; n < 4; n++)
            projectNames[n] = "project" + n;
        projectName = projectNames[0];
    }

    private void createProjectObjects()
    {
        objectNames = new String[4][3];
        for (int n = 0; n < 4; n++)
        {
            for (int i = 0; i < 3; i++)
                objectNames[n][i] = projectNames[n] + "-object" + i;
        }
    }

    public void setDatastoreName(String datastoreName)
    {
        name = datastoreName;
    }

    public void setDatastoreNames(String[] datastoreNames)
    {
        name = datastoreNames[0];
    }

    public String getDatastoreName() { return name; }
    public String[] getDatastoreNames() { return new String[] { name }; }

    public Object[] getObjects() { return objectNames; }

    public boolean initialize() { return true; }

    public String[] getProjectNames() { return projectNames; }

    public void setProjectName(String project)
    {
        this.projectName = project;
    }

    public String getProjectName() { return projectName; }

    public Object[] getProjectObjects()
    {
        for (int n = 0; n < projectNames.length; n++)
            if (projectNames[n] == projectName)
            {
                return objectNames[n];
            }
        return new Object[0];
    }

    public Object[] getAvailableObjects()
     {
         return objectNames;
     }

     public Object[] initializeAvailableObjects()
     {
         initialize();
         return objectNames;
     }
}

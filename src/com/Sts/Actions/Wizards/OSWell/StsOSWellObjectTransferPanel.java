package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Import.*;
import com.Sts.DBTypes.OpenSpirit.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 10, 2008
 * Time: 7:33:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsOSWellObjectTransferPanel extends StsDatastoreObjectTransferPanel
{
    JFrame parent = null;
    StsProgressPanel progressPanel = null;
    private String wbListName = "none";
    private String currentWbList = null;
    private JFileChooser chooseDirectory = null;

    private StsButton wbListBrowseButton = new StsButton("dir16x32", "Browse for Wellbore List.", this, "wbListBrowse");
    private StsStringFieldBean currentWbListBean = null;
    private StsButton refreshBtn = new StsButton("Refresh", "Press to refresh well selection list.", this, "refreshWellList");

    public StsOSWellObjectTransferPanel(StsOSWellDatastore datastore, StsObjectTransferListener listener, int width, int height)
    {
        this(datastore, listener, null, null, width, height);
    }

    public StsOSWellObjectTransferPanel(StsOSWellDatastore datastore, StsObjectTransferListener listener, JFrame parent, StsProgressPanel progressPanel, int width, int height)
    {
        super(datastore, listener, width, height);
        this.datastore = datastore;
        this.parent = parent;
        this.progressPanel = progressPanel;
        currentWbListBean = new StsStringFieldBean(this, "currentWbList", "Wellbore List: ");
        setCurrentWbList("None");
        directoryGroupBox.gbc.gridwidth = 1;
        directoryGroupBox.addToRow(wbListBrowseButton);
        currentWbListBean.setColumns(30);
        directoryGroupBox.addEndRow(currentWbListBean);
        directoryGroupBox.addEndRow(refreshBtn);
        currentWbListBean.setEditable(true);
    }

    //public StsDatastoreFace getDatastore()
    public void initializeDatastoreName()
    {
        if (progressPanel != null)
            progressPanel.setDescription("Retrieving projects from selected datastore. Please standby...");
        else
            StsMessageFiles.infoMessage("Retrieving projects from selected datastore. Please standby...");
        super.initializeDatastoreName();
    }

    public void initializeProjectName()
    {
        parent.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        if (progressPanel != null)
            progressPanel.setDescription("Getting wells from selected project. Please standby...");
        else
            StsMessageFiles.infoMessage("Getting wells from selected project. Please standby...");

        super.initializeProjectName();
        parent.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
    }

    public StsDatastoreFace getDatastore()
    {
        if (parent != null) parent.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        StsOpenSpiritImport ospImport = ((StsOSDatastore) datastore).getOpenSpiritImport();

        if (progressPanel != null)
            progressPanel.setDescription("Getting data servers from master installation. Please standby...");
        else
            StsMessageFiles.infoMessage("Getting data servers from master installation. Please standby...");
        String[] stores = ospImport.getOspDataServers();
        if (stores.length == 0)
        {
            if (parent != null) parent.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            return null;
        }
        StsComboBoxDialog dlg = new StsComboBoxDialog(parent, "Select OS Datastore", "OpenSpirit Datastores", stores, this, "datastoreName");
        String selStore = (String) dlg.getSelectedItem();
//		StsComboBoxFieldBean comboBoxFieldBean = new StsComboBoxFieldBean(StsComboBoxDialog.class, "Datastore", "OpenSpirit Datastores", stores);
//		StsComboBoxDialog dlg = new StsComboBoxDialog(parent, "Select an OpenSpirit Datastore", "Datastores:", stores, );
//		comboBoxFieldBean.setBeanObject(dlg);
//		String selStore = (String)comboBoxFieldBean.getValueObject();

        // parse the datastore name from the selected string of dsType dsVersion - name
        int iColon = selStore.indexOf("-");
        String name = selStore.substring(iColon + 2);

        super.setDatastoreName(name);

        if (parent != null) parent.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        return datastore;
    }

    public void wbListBrowse()
    {
        if (chooseDirectory == null) initializeChooseDirectory();

        chooseDirectory = new JFileChooser(currentWbList);
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_ONLY);
        chooseDirectory.setDialogTitle("Select or Enter Desired Wellbore List and Press Open");
        chooseDirectory.setApproveButtonText("Open Wellbore List");
        while (true)
        {
            int retVal = chooseDirectory.showOpenDialog(null);
            if (retVal != chooseDirectory.APPROVE_OPTION)
                break;
            File newFile = chooseDirectory.getSelectedFile();
            if (newFile == null) return;
            if (newFile.isFile())
            {
                setCurrentWbList(newFile.getAbsolutePath());
                break;
            }
            else
            {
                // Directory or nothing selected, strip off file name
                if (!StsYesNoDialog.questionValue(this, "Must select a wellbore list.\n\n Continue?"))
                    break;
            }
        }
    }

    private void initializeChooseDirectory()
    {
        chooseDirectory = new JFileChooser(currentWbList);
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_ONLY);
    }

    public void setCurrentWbList(String list)
    {
        if (list == null)
            return;
        currentWbList = list;
        if (currentWbListBean != null)
            currentWbListBean.setValue(currentWbList);

        if (list.equals("None")) return;

        // Open selected file and fill list.
        //setListFromWbList(currentWbList);
    }

    public String getCurrentWbList()
    {
        return currentWbList;
    }

    public String[] readWbListFromFile(File wbListFile)
    {
        return new String[]{"17470", "17471"};
    }

    public void setListFromWbList(String list)
    {
        File wbListFile = new File(list);
        String[] wbList = readWbListFromFile(wbListFile);
        if (wbList == null)
        {
            new StsMessage(this, StsMessage.WARNING, "No wellbores found in selected file.");
        }

        if (wbList.length > 1)
        {
            ((StsOSWellDatastore) datastore).runGetProjectObjects(this, progressPanel, new String[]{list}); // All wells in list
            //datastore.runGetProjectObjects(this, progressPanel, wbList); // All wells in list
            return;
        }
    }

    public void refreshWellList()
    {
        if (currentWbList.equals("None"))
            ((StsOSWellDatastore) datastore).runGetProjectObjects(this, progressPanel); // All wells in project
        else
        {
            // Open selected file and fill list.
            setListFromWbList(currentWbList);
        }
    }

    public void setAvailableFromWbList(String wbListName)
    {
        StsOpenSpiritImport ospImport = ((StsOSDatastore) datastore).getOpenSpiritImport();
        setAvailableObjects(ospImport.getWellsFromWellboreLists(new String[]{wbListName}));
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
}

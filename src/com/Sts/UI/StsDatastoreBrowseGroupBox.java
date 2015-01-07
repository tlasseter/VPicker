package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 1, 2008
 * Time: 3:12:17 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsDatastoreBrowseGroupBox extends StsGroupBox
{
    private JFileChooser chooseDirectory = null;
    private StsButton directoryBrowseButton = new StsButton("Add datastore", "Browse for new datastore.", this, "addDatastore");
    private StsComboBoxFieldBean datastoreComboBoxBean;
    private StsDatastoreFace[] datastores;

    static final String nullDatastoreName = "None";

    abstract public StsDatastoreFace[] getDatastores();

    public StsDatastoreBrowseGroupBox(String groupBoxName, Object transferPanel, String transferPanelDatastoreFieldName)
    {
        super(groupBoxName);
        
        datastores = getDatastores();
        if(datastores == null)
            datastoreComboBoxBean = new StsComboBoxFieldBean(transferPanel, transferPanelDatastoreFieldName, null, new String[] { nullDatastoreName });
        else
            datastoreComboBoxBean = new StsComboBoxFieldBean(transferPanel, transferPanelDatastoreFieldName, null, datastores);

        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.weightx = 0.0;
        addToRow(directoryBrowseButton);
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.EAST;
        addEndRow(datastoreComboBoxBean);
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weighty = 0.0;
    }

    public void addDatastore()
    {
        if (chooseDirectory == null) initializeChooseDirectory();
        StsModel model = StsModel.getCurrentModel();
        String currentDirectory = model.getProject().getRootDirString();
        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        chooseDirectory.setDialogTitle("Select or Enter new datastore and Press Open");
        chooseDirectory.setApproveButtonText("Open Datastore");
        while(true)
        {
            int retVal = chooseDirectory.showOpenDialog(null);
            if(retVal != chooseDirectory.APPROVE_OPTION)
                break;
            File newDirectory = chooseDirectory.getSelectedFile();
            if(newDirectory == null) continue;
            if(newDirectory.isDirectory())
            {
                setCurrentDirectory(newDirectory.getAbsolutePath());
                break;
            }
            else
            {
                // File or nothing selected, strip off file name
                String dirString = newDirectory.getPath();
                newDirectory = new File(dirString.substring(0,dirString.lastIndexOf(File.separator)));
                if(newDirectory.isDirectory())
                {
                    setCurrentDirectory(newDirectory.getAbsolutePath());
                    break;
                }
                if(!StsYesNoDialog.questionValue(this,"Must select the datastore.\n Continue?"))
                    break;
            }
        }
    }

    private void initializeChooseDirectory()
    {
        chooseDirectory = new JFileChooser(getCurrentDirectory());
        chooseDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    }

    private String getCurrentDirectory()
    {
        return (String)datastoreComboBoxBean.getValue();
    }

    private void setCurrentDirectory(String directory)
    {
        datastoreComboBoxBean.doSetValueObject(directory);
        datastoreComboBoxBean.setBeanObjectValue();
    }
}
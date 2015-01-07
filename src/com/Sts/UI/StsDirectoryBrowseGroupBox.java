package com.Sts.UI;

import com.Sts.UI.Beans.*;

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
public class StsDirectoryBrowseGroupBox extends StsGroupBox
{
    private JFileChooser chooseDirectory = null;
    private StsButton directoryBrowseButton = new StsButton("Get directory", "Browse for directory.", this, "browseDirectory");
    private StsStringFieldBean directoryPathnameBean;

    public StsDirectoryBrowseGroupBox(Object directoryObject, String directoryNameFieldName)
    {
        super("Directory selector");
        directoryPathnameBean = new StsStringFieldBean(directoryObject, directoryNameFieldName, null);
        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.weightx = 0.0;
        addToRow(directoryBrowseButton);
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.EAST;
        directoryPathnameBean.setColumns(30);
        addEndRow(directoryPathnameBean);
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weighty = 0.0;
    }

    public void browseDirectory()
    {
        if (chooseDirectory == null)
            initializeChooseDirectory();

        chooseDirectory = new JFileChooser(directoryPathnameBean.getValue());
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Open");
        chooseDirectory.setApproveButtonText("Open Directory");
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
                if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the Surface Files.\n\n Continue?"))
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
        return directoryPathnameBean.getValue();
    }

    private void setCurrentDirectory(String directory)
    {
        directoryPathnameBean.setValue(directory);
        directoryPathnameBean.setBeanObjectValue();
    }

    static public void main(String[] args)
    {
        String directory = "c:";

    }
}

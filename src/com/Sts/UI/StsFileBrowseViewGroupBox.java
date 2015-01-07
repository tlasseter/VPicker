package com.Sts.UI;

import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 4, 2008
 * Time: 11:40:13 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsFileBrowseViewGroupBox extends StsFileBrowseGroupBox
{
    StsButton viewButton;

    public StsFileBrowseViewGroupBox(String title, String directoryName, String fileFilter, boolean isPrefix)
    {
        super(title, directoryName, fileFilter, isPrefix);
    }

    protected void constructPanel()
    {
        gbc.fill = GridBagConstraints.NONE;
        gbc.weightx = 0.0;
        addToRow(fileBrowseButton);
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        filePathnameBean.setColumns(30);
        addToRow(filePathnameBean);
        gbc.fill = GridBagConstraints.NONE;
        gbc.weightx = 0.0;
        viewButton = new StsButton("View", "View selected file(s)", this, "view");
        addEndRow(viewButton);
    }

    public void view()
    {
        File[] selectedFiles = fileChooser.getFiles();
        StsTabbedTextAreaDialog dialog = new StsTabbedTextAreaDialog(null, "View files", true, 500, 500, selectedFiles);
        dialog.setVisible(true);
    }

    static public void main(String[] args)
    {
        StsFileBrowseViewGroupBox browseViewGroupBox = new StsFileBrowseViewGroupBox("test", ".", ".txt", false);
        StsToolkit.createDialog(browseViewGroupBox);
    }
}
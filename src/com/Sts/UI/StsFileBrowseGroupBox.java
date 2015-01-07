package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.filechooser.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 1, 2008
 * Time: 3:12:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFileBrowseGroupBox extends StsGroupBox
{
    private String filename;
    FileFilter fileFilter;
    protected StsFileChooser fileChooser = null;
    protected StsButton fileBrowseButton = new StsButton("Get file", "Browse for file.", this, "browseFile");
    protected StsStringFieldBean filePathnameBean;

    public StsFileBrowseGroupBox(String title, String directoryName, String fileFilter, boolean isPrefix)
    {
        super(title);
        fileChooser = StsFileChooser.createFileChooser(null, "Choose file", directoryName, fileFilter, isPrefix, false);
        filePathnameBean = new StsStringFieldBean(this, "filename", null);
        constructPanel();
    }

    protected void constructPanel()
    {
        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.weightx = 0.0;
        addToRow(fileBrowseButton);
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.EAST;
        filePathnameBean.setColumns(30);
        addToRow(filePathnameBean);
    }

    public void browseFile()
    {
        fileChooser.show();
        setFilename(fileChooser.getFilename());
        setCurrentFile(getFilename());
    }

    private void setCurrentFile(String filename)
    {
        filePathnameBean.setValue(filename);
    }

    public String getFilename()
    {
        return filename;
    }

    public void setFilename(String filename)
    {
        this.filename = filename;
    }

    public String getPathname() { return fileChooser.getFilePath(); }

    static public void main(String[] args)
    {
        String directoryName = ".";
        StsFileBrowseGroupBox browseGroupBox = new StsFileBrowseGroupBox("test", directoryName, "sensor.PVT", true);
        StsToolkit.createDialog(browseGroupBox);
    }
}
package com.Sts.Actions.Wizards.VspLoad;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2006
//Author:       TJLasseter
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

public class StsVspSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
    private boolean archiveIt = false;
    private StsFileTransferPanel selectionPanel;

    private StsFilenameFilter filenameFilter;

    private String currentDirectory = null;
    private JTextPane volumeInfoLbl = new JTextPane();
    private JPanel selectedVolInfoPanel = new JPanel();

    private StsWizard wizard;
    private StsWizardStep wizardStep;
    private StsAbstractFile[] selectedFiles;
    private StsModel model = null;

    private StsVsp seismicVolume = null;
    private StsSpectrum spectrum = null;

    StsGroupBox paletteBox = new StsGroupBox("Palette Selection");
	StsButton selectPaletteButton = new StsButton("Select Color Palette", "Select color palette.", this, "selectPalette");

    static private final String group = StsSeismicBoundingBox.groupVsp;
    static private final String format = StsSeismicBoundingBox.headerFormat;
    static private String headerFilePrefix = group + "." + format + ".";

    public StsVspSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void constructBeans()
    {
        model = wizard.getModel();
        if(model != null)
            currentDirectory = wizard.getModel().getProject().getRootDirString();
        else
            currentDirectory = System.getProperty("user.dirNo"); // standalone testing
        filenameFilter = new StsFilenameFilter(group, format);
        selectionPanel = new StsFileTransferPanel(currentDirectory, filenameFilter, this);
    }

    public StsAbstractFile[] getSelectedFiles()
    {
        if(selectedFiles != null) return selectedFiles;
        else return new StsFile[0];
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.BOTH;
        gbc.anchor = gbc.NORTH;
        gbc.weighty = 1.0;
        addEndRow(selectionPanel);

        paletteBox.addEndRow(selectPaletteButton);
        addEndRow(paletteBox);

        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        addEndRow(selectedVolInfoPanel);
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        String volName = null;

/*
        if(source == selectionPanel)
        {
            if(selectedFiles.length > 0)
                volName = selectionPanel.getSelectedFiles()[0].getFilenameStem();
        }
*/        

        // Get PostStack3d Information
        if(volName != null)
        {
            try
            {
                StsParameterFile.initialReadObjectFields(currentDirectory + File.separator + headerFilePrefix + volName, seismicVolume, StsSeismicBoundingBox.class, StsBoundingBox.class);
            }
            catch(Exception ex)
            {
                StsException.outputException("StsVolumeSelectPanel.valueChanged failed.", ex, StsException.WARNING);
                return;
            }
            String desc = new String("Name: " + volName + "\n" +
               "SegY File: " + seismicVolume.segyFilename + "\n" +
               "SegY Last Modified On: " + new Date(seismicVolume.segyLastModified).toString() + "\n" +
               "Data Range: " + String.valueOf(seismicVolume.dataMin) + " to " + String.valueOf(seismicVolume.dataMax) + "\n" +
               "X Origin: " + String.valueOf(seismicVolume.xOrigin) + " Y Origin: " + String.valueOf(seismicVolume.yOrigin));
            volumeInfoLbl.setText(desc);
        }
        else
            volumeInfoLbl.setText("No PostStack3d Selected");
        return;
    }

    public String[] getFilenameEndings(StsAbstractFile[] files)
    {
        int nFiles = files.length;
        String[] fileEndings = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            fileEndings[n] = filenameFilter.getFilenameEnding(files[n].getFilename());
        return fileEndings;
    }

	public void selectPalette()
	{
		StsSpectrumSelect ss = new StsSpectrumSelect(null);
		ss.setVisible(true);
		spectrum = ss.getSelectedSpectrum();
		selectPaletteButton.setText(spectrum.getName());
	}

    public StsSpectrum getSpectrum()
    {
        return spectrum;
    }
    
    public void availableFileSelected(StsFile selectedFile)
    {
    }
    
    public boolean hasDirectorySelection() { return true;  }
    public boolean hasReloadButton() { return false;  }
    public boolean hasOverrideButton() { return false; }
    public boolean hasArchiveItButton() { return false; }
    public void setArchiveIt(boolean archive) { archiveIt = archive; }
    public boolean getArchiveIt() { return archiveIt; }
    public void setOverrideFilter(boolean override) {  }
    public boolean getOverrideFilter() { return false; }
    public void setReload(boolean reload) { }
    public boolean getReload() { return false; }
    public void fileSelected(StsAbstractFile selectedFile) { }
    public void availableFileSelected(StsAbstractFile selectedFile) { }

    public void addFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
            selectedFiles = (StsFile[]) StsMath.arrayAddElement(selectedFiles, files[i]);
    }

    public void removeFiles(StsAbstractFile[] files)
    {
        for(int i=0; i<files.length; i++)
            selectedFiles = (StsFile[]) StsMath.arrayDeleteElement(selectedFiles, files[i]);
    }

    public void removeAllFiles()
    {
        if(selectedFiles == null) return;
       for(int i = 0; i<selectedFiles.length; i++)
           selectedFiles = (StsFile[]) StsMath.arrayDeleteElement(selectedFiles, selectedFiles[i]);

        selectedFiles = null;
    }

	public static void main(String[] args)
	{
		StsModel model = StsModel.constructDefaultMainWindow();
		StsActionManager actionManager = new StsActionManager(model);
		StsVspLoadWizard wizard = new StsVspLoadWizard(actionManager);
		StsVspSelectPanel panel = new StsVspSelectPanel(wizard, new StsVspSelect(wizard));
		StsToolkit.createDialog(panel, true, 200, 500);
	}
}

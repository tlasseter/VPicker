package com.Sts.Actions.Wizards.PostStack3dLoad;

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
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsVolumeSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew, ActionListener
{
	private StsFileTransferPanel selectionPanel;

	private StsFilenameFilter filenameFilter;

	private String currentDirectory = null;
	private StsGroupBox palettePanel = new StsGroupBox();
	private JTextPane volumeInfoLbl = new JTextPane();
	private StsGroupBox selectedVolInfoPanel = new StsGroupBox();
	private JTextField selectedPalette = new JTextField();
	private JButton paletteBtn = new JButton();

	private StsWizard wizard;
	private StsFile[] selectedFiles;
	private StsModel model = null;

	private StsSeismicVolume seismicVolume = null;
	private StsSpectrum spectrum = null;

	transient static String headerFilePrefix = "seis3d.txt.";

	public StsVolumeSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;

		try
		{
			constructBeans();
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private void constructBeans()
	{
		model = wizard.getModel();
		if (model == null)
			currentDirectory = System.getProperty("user.dirNo"); // standalone testing
		else
			currentDirectory = model.getProject().getRootDirString();

		filenameFilter = new StsFilenameFilter(StsSeismicBoundingBox.group3d, StsSeismicBoundingBox.headerFormat);
		seismicVolume = new StsSeismicVolume(false);
		selectionPanel = new StsFileTransferPanel(currentDirectory, filenameFilter, this, 300, 100, false);
	}

	public StsFile[] getSelectedFiles()
	{
		if (selectedFiles != null)
			return selectedFiles;
		else
			return new StsFile[0];
	}

	void jbInit() throws Exception
	{
		gbc.fill = gbc.BOTH;
		gbc.anchor = gbc.NORTH;
		gbc.weighty = 0.25;
		addEndRow(selectionPanel);

		selectedPalette.setBackground(UIManager.getColor("Menu.background"));
		selectedPalette.setBorder(BorderFactory.createEtchedBorder());
		selectedPalette.setText("Set Default Color Palette");
		paletteBtn.setText("Select Palette");
		paletteBtn.addActionListener(this);

		gbc.fill = gbc.BOTH;
		gbc.weighty = 0.0;
		palettePanel.gbc.fill = gbc.BOTH;
		palettePanel.gbc.weighty = 1.0;
		addEndRow(palettePanel);
		palettePanel.addToRow(selectedPalette, 2, 1.0);
		palettePanel.addEndRow(paletteBtn, 1, 0.0);

		gbc.fill = gbc.BOTH;
		gbc.weighty = 0.75;
		selectedVolInfoPanel.gbc.fill = gbc.BOTH;
		selectedVolInfoPanel.gbc.weighty = 1.0;
		selectedVolInfoPanel.addEndRow(volumeInfoLbl);
		addEndRow(selectedVolInfoPanel);
		volumeInfoLbl.setText("Selected PostStack3d Information");
		selectedVolInfoPanel.setMinimumSize(new Dimension(400, 100));
		volumeInfoLbl.setEditable(false);
		volumeInfoLbl.setBackground(Color.lightGray);
		volumeInfoLbl.setFont(new Font("Dialog", 0, 10));
	}

	public String[] getFilenameEndings(StsAbstractFile[] files)
	{
		int nFiles = files.length;
		String[] fileEndings = new String[nFiles];
		for (int n = 0; n < nFiles; n++)
			fileEndings[n] = filenameFilter.getFilenameEnding(files[n].getFilename());
		return fileEndings;
	}
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }
    
	public StsSpectrum getSpectrum()
	{
		return spectrum;
	}

	public void actionPerformed(ActionEvent e)
	{
		Object source = e.getSource();

		if (source == paletteBtn)
		{
			StsSpectrumSelect ss = new StsSpectrumSelect(null);
			ss.setVisible(true);
			spectrum = ss.getSelectedSpectrum();
			selectedPalette.setText(spectrum.getName());
		}
	}

	public void addFiles(StsAbstractFile[] files)
	{
		for (int i = 0; i < files.length; i++)
			selectedFiles = (StsFile[])StsMath.arrayAddElement(selectedFiles, files[i]);
	}

	public void removeFiles(StsAbstractFile[] files)
	{
		for (int i = 0; i < files.length; i++)
			selectedFiles = (StsFile[])StsMath.arrayDeleteElement(selectedFiles, files[i]);
	}

	public void removeAllFiles()
	{
		if (selectedFiles == null)
			return;
		for (int i = 0; i < selectedFiles.length; i++)
			selectedFiles = (StsFile[])StsMath.arrayDeleteElement(selectedFiles, selectedFiles[i]);

		selectedFiles = null;
	}

	public boolean hasDirectorySelection() { return true;  }
	public boolean hasReloadButton() { return false;  }
	public boolean hasOverrideButton() { return false; }
	public boolean hasArchiveItButton() { return false; }
	public void setArchiveIt(boolean archive) {}
	public boolean getArchiveIt() { return false; }
	public void setOverrideFilter(boolean override) {  }
	public boolean getOverrideFilter() { return false; }
	public void setReload(boolean reload) { }
	public boolean getReload() { return false; }

	public void fileSelected(StsAbstractFile selectedFile)
	{
		String volName = selectedFile.getFilename();

		// Get PostStack3d Information
		if (volName != null)
		{
			try
			{
				//StsParameterFile.readObjectFields(currentDirectory + File.separator + headerFilePrefix + volName, lineSet, StsSeismicBoundingBox.class, StsBoundingBox.class);
				StsParameterFile.initialReadObjectFields(currentDirectory + File.separator + volName, seismicVolume, StsSeismicBoundingBox.class, StsBoundingBox.class);
			}
			catch (Exception ex)
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
	}
}

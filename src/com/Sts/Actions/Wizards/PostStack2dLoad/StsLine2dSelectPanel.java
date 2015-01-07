package com.Sts.Actions.Wizards.PostStack2dLoad;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
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

public class StsLine2dSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew, ActionListener
{
	private StsFileTransferPanel selectionPanel;

	private StsFilenameFilter filenameFilter;
	private String group = StsSeismicBoundingBox.group2d;
	private String format = StsSeismicBoundingBox.headerFormat;

	private String currentDirectory = null;
	private StsGroupBox palettePanel = new StsGroupBox();
	private JTextPane volumeInfoLbl = new JTextPane();
	private StsGroupBox selectedVolInfoPanel = new StsGroupBox();
	private JTextField selectedPalette = new JTextField();
	private JButton paletteBtn = new JButton();

	private StsWizard wizard;
	private StsFile[] selectedFiles;
	private StsModel model = null;

	private StsSeismicLine2d seismicLine2d = null;
	private StsSpectrum spectrum = null;

	transient static String headerFilePrefix = "seis2d.txt.";

	public StsLine2dSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
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

		filenameFilter = new StsFilenameFilter(group, format);
		seismicLine2d = new StsSeismicLine2d(false);
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
		volumeInfoLbl.setText("Selected PostStack2d Information");
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

	public StsSpectrum getSpectrum()
	{
		return spectrum;
	}
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
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
			selectedFiles = (StsFile[])StsMath.arrayAddElement(selectedFiles, (StsFile)files[i]);
	}

	public void removeFiles(StsAbstractFile[] files)
	{
		for (int i = 0; i < files.length; i++)
			selectedFiles = (StsFile[])StsMath.arrayDeleteElement(selectedFiles, (StsFile)files[i]);
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

	public void fileSelected(StsAbstractFile abstractSelectedFile)
	{
        StsFile selectedFile = (StsFile)abstractSelectedFile;
        String volName = selectedFile.getFilenameStem();

		// Get PostStack3d Information
		if (volName != null)
		{
			try
			{
				StsParameterFile.initialReadObjectFields(currentDirectory + File.separator + headerFilePrefix + volName, seismicLine2d, StsPreStackLine.class, StsBoundingBox.class);
			}
			catch (Exception ex)
			{
				StsException.outputException("StsLine2dSelectPanel.valueChanged failed.", ex, StsException.WARNING);
				return;
			}
			String desc = new String("Name: " + volName + "\n" +
											 "SegY File: " + seismicLine2d.segyFilename + "\n" +
											 "SegY Last Modified On: " + new Date(seismicLine2d.segyLastModified).toString() + "\n" +
											 "Data Range: " + String.valueOf(seismicLine2d.dataMin) + " to " + String.valueOf(seismicLine2d.dataMax) + "\n" +
											 "X Origin: " + String.valueOf(seismicLine2d.xOrigin) + " Y Origin: " + String.valueOf(seismicLine2d.yOrigin));
			volumeInfoLbl.setText(desc);
		}
		else
			volumeInfoLbl.setText("No PostStack2d Selected");
	}
}


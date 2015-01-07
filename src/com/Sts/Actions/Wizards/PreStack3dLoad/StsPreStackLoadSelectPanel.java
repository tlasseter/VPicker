package com.Sts.Actions.Wizards.PreStack3dLoad;

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

public class StsPreStackLoadSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew, ActionListener
{
	private StsFileTransferPanel selectionPanel;

	private HeaderFilenameFilter filenameFilter;

	private String currentDirectory = null;
	private StsGroupBox palettePanel = new StsGroupBox();
	private JTextPane volumeInfoLbl = new JTextPane();
	private StsGroupBox selectedVolInfoPanel = new StsGroupBox();
	private StsGroupBox predecessorGroup = new StsGroupBox();
	private StsComboBoxFieldBean predecessorBean = new StsComboBoxFieldBean();
	private JTextField selectedPalette = new JTextField();
	private JButton paletteBtn = new JButton();

	private StsWizard wizard;
	private StsFile[] selectedFiles;
	private StsModel model = null;

	private StsPreStackLine3d seismicLine = null;
	private StsSpectrum spectrum = null;

	private StsStringFieldBean nameBean = new StsStringFieldBean();

	private String volName= "PreStack3d-0";

	static final String format = StsPreStackLine.headerFormat;

	public StsPreStackLoadSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
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
		{
		    currentDirectory = model.getProject().getDirectoryPath(StsProject.PreStackSeismicDirectory);  //use data directory if it has been set SWC 9/23/09
		    if (currentDirectory == null) currentDirectory = model.getProject().getRootDirString();
		}

		filenameFilter = new HeaderFilenameFilter(StsSeismicBoundingBox.group3dPrestack, format);
		seismicLine = new StsPreStackLine3d(false);
		selectionPanel = new StsFileTransferPanel(currentDirectory, filenameFilter, this, 300, 100, false);

		Object[] availableVolumes = model.getTrimmedList(StsPreStackLineSet3d.class);
		predecessorBean.initialize(wizard, "predecessor", "Previous PreStack3d:", availableVolumes);
		if (availableVolumes.length == 0)
			predecessorBean.setEditable(false);
		else
		{
			predecessorBean.setEditable(true);
			predecessorBean.setSelectedIndex(availableVolumes.length - 1);
		}

		// Derive next default name
		if (availableVolumes.length > 0)
		{
			StsPreStackLineSet3d lastVol = (StsPreStackLineSet3d)availableVolumes[availableVolumes.length - 1];
			setName("PreStack3d-" + (lastVol.getIndex() + 1));
			nameBean.setValue(volName);
		}
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

		nameBean.initialize(this, "name", volName, true, "PreStack3d Name:");

		gbc.fill = gbc.BOTH;
		gbc.weighty = 0.0;
		palettePanel.gbc.fill = gbc.BOTH;
		palettePanel.gbc.weighty = 1.0;
		addEndRow(palettePanel);
		palettePanel.addEndRow(nameBean, 3, 1.0);
		palettePanel.addToRow(selectedPalette, 2, 1.0);
		palettePanel.addEndRow(paletteBtn, 1, 0.0);
		predecessorGroup.add(predecessorBean);
		addEndRow(predecessorGroup);

		gbc.fill = gbc.BOTH;
		gbc.weighty = 0.75;
		selectedVolInfoPanel.gbc.fill = gbc.BOTH;
		selectedVolInfoPanel.gbc.weighty = 1.0;
		selectedVolInfoPanel.addEndRow(volumeInfoLbl);
		addEndRow(selectedVolInfoPanel);
		volumeInfoLbl.setText("Selected PreStack3d Information");
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

	public String getName()
	{
		return volName;
	}

	public void setName(String name)
	{
		volName = name;
	}

	public void addFiles(StsAbstractFile[] files)
	{
		for(int i=0; i<files.length; i++)
			 selectedFiles = (StsFile[]) StsMath.arrayAddElement(selectedFiles, files[i]);
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
				StsParameterFile.initialReadObjectFields(currentDirectory + File.separator + volName, seismicLine, StsPreStackLine.class, StsBoundingBox.class);
			}
			catch (Exception ex)
			{
				StsException.outputException("StsPreStackLoadSelectPanel.valueChanged failed.", ex, StsException.WARNING);
				return;
			}
			String desc = new String("Name: " + volName + "\n" +
											 "SegY File: " + seismicLine.segyFilename + "\n" +
											 "SegY Last Modified On: " + new Date(seismicLine.segyLastModified).toString() + "\n" +
											 "Data Range: " + String.valueOf(seismicLine.dataMin) + " to " + String.valueOf(seismicLine.dataMax) + "\n" +
											 "X Origin: " + String.valueOf(seismicLine.xOrigin) + " Y Origin: " + String.valueOf(seismicLine.yOrigin));
			volumeInfoLbl.setText(desc);
		}
		else
			volumeInfoLbl.setText("No PreStack3d Selected");
	}

	private boolean isObjectAlreadyLoaded(String name)
	{
		StsObject[] seismicVolumes = model.getObjectList(StsPreStackLineSet3d.class);
		for (int i = 0; i < seismicVolumes.length; i++)
		{
			StsPreStackLineSet3d vol = (StsPreStackLineSet3d)seismicVolumes[i];
			if (vol.isLineInVolume(name))
				return true;
		}
		return false;
	}

	final class HeaderFilenameFilter extends StsFilenameFilter
	{
		public HeaderFilenameFilter(String group, String format)
		{
            super(group, format);
		}

		public boolean accept(File dir, String filename)
		{
			if(!super.accept(dir, filename)) return false;
            String name = getFilenameEnding(filename);
            if(name == null) return false;
            return !isObjectAlreadyLoaded(name);
		}
	}
}

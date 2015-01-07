package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsPostStackSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
	// Going to hand build the global parameters panel in this class.
	private StsComboBoxFieldBean segyFormatComboBoxBean;
	private StsDoubleFieldBean scanPercentBean;
	private StsFloatFieldBean diskRequiredBean;

	private StsGroupBox segyFormatPanel = new StsGroupBox("Default Parameters");
	private StsFileTransferPanel selectionPanel;

	private String segyFormatFilename;

	// Z Domain
	public static final int TIME = 1;
	public static final int DEPTH = 2;

	private StsSeismicWizard wizard;
	private StsPostStackSelect wizardStep;
	private StsModel model;
	private String outputDirectory = null;

	private StsSeismicBoundingBox currentVolume = null;
	private StsFile currentSelectedFile = null;

    public static final String NO_SEGY_FORMAT_ASSIGNED = "Default";
    public static final String NO_SEGY_FILENAME = "none";

    public StsPostStackSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = (StsSeismicWizard)wizard;
		model = wizard.getModel();
		this.wizardStep = (StsPostStackSelect)wizardStep;
		try
		{
			String currentDirectory;
			if (model != null)
			{
				StsProject project = model.getProject();
				currentDirectory = project.getRootDirString();
				outputDirectory = project.getDataFullDirString();
			}
			else
			{
				currentDirectory = System.getProperty("user.dirNo"); // standalone testing
				outputDirectory = currentDirectory;
			}
			String[] filterStrings = new String[]
				{"sgy", "segy", "Segy", "SegY"};
			StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);

			selectionPanel = new StsFileTransferPanel(currentDirectory, filter, this); //, 300, 100);

			constructBeans();
			constructPanel();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private void constructBeans()
	{
		scanPercentBean = new StsDoubleFieldBean(wizard, "defaultScanPercent", 0.0, 100.0, "Scan %:", false);
		scanPercentBean.fixStep(0.01);
		segyFormatComboBoxBean = new StsComboBoxFieldBean(this, "segyFormatFilename", "SEGY format:");
		diskRequiredBean = new StsFloatFieldBean(this, "diskRequired", false, "Disk Required (MB):");
	}

	private void constructPanel() throws Exception
	{
		Color panelColor = Color.black;
		Font panelFont = new Font("Dialog", Font.PLAIN, 11);

		segyFormatPanel.setForeground(panelColor);
		segyFormatPanel.setFont(panelFont);

		segyFormatPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		segyFormatPanel.addToRow(segyFormatComboBoxBean);
		segyFormatPanel.addEndRow(scanPercentBean);
		segyFormatPanel.addToRow(diskRequiredBean);

		segyFormatComboBoxBean.setEditable(false);

		gbc.anchor = GridBagConstraints.NORTH;
		gbc.weighty = 1.0;
		gbc.fill = GridBagConstraints.BOTH;
		add(selectionPanel);
		gbc.weighty = 0.0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		add(segyFormatPanel, 1, 0.0);
	}

	public void initialize()
	{
		if (selectionPanel.getSelectedCount() > 0) return;
        byte segyFormatDataType = wizard.getSegyFormatDataType();
        initializeSegyFormats(segyFormatDataType);
		initializePanels();

		gbc.fill = GridBagConstraints.BOTH;
	}

	private void initializePanels()
	{
		if(segyFormatComboBoxBean != null) segyFormatComboBoxBean.getValueFromPanelObject();
		scanPercentBean.getValueFromPanelObject();
	}
    
    public void availableFileSelected(StsFile selectedFile)
    {
    }
    
	private void initializeSegyFormats(byte segyFormatDataType)
	{
        segyFormatComboBoxBean.removeAll();

        String[] segyFormatFilenames = StsSEGYFormat.getSegyFormatFilenames();
        int nFiles = segyFormatFilenames.length;
        segyFormatComboBoxBean.addItem(NO_SEGY_FORMAT_ASSIGNED);
        byte dataType = wizard.getSegyFormatDataType();
        for (int i = 0; i < nFiles; i++)
        {
            String stemname = segyFormatFilenames[i];
            if(StsSEGYFormat.getFilenameType(stemname) == dataType)
                segyFormatComboBoxBean.addItem(stemname);
        }
        segyFormatComboBoxBean.setEditable(nFiles > 0);
        segyFormatComboBoxBean.setSelectedIndex(0);
	}

    public String getSegyFormatFilename()
	{
		return segyFormatFilename;
	}

	public void setSegyFormatFilename(String filename)
	{
		StsSEGYFormat segyFormat;
        byte dataType = wizard.getSegyFormatDataType();
        wizard.setSkipReanalyzeTraces(true);
		try
		{
			if (filenameNotAssigned(filename))
                segyFormat = StsSEGYFormat.constructor(model, dataType);
			else
			{
				segyFormat = StsSEGYFormat.constructor(model, filename, dataType);
				if (segyFormat == null)
				{
					initializePanels();
					return;
				}
			}
			segyFormatFilename = filename;
			wizard.setSegyFormat(segyFormat);
			initializePanels();
		}
        catch(Exception e)
        {
            StsException.outputWarningException(this, "setSegyFormatFilename", e);
        }
		finally
		{
			wizard.setSkipReanalyzeTraces(false);
		}
	}

    private boolean filenameNotAssigned(String filename)
    {
        return filename == null || filename.equals(NO_SEGY_FORMAT_ASSIGNED) || filename.equals(NO_SEGY_FILENAME);
    }

    public void addFiles(StsAbstractFile[] files)
	{
		if (files == null || files.length == 0)
		{
			return;
		}
		wizard.setSkipReanalyzeTraces(true);
		StsFile firstAddedFile = null;
		for (int n = 0; n < files.length; n++)
		{
            StsFile file = (StsFile)files[n];
            if (wizard.addSegyVolume(file, outputDirectory) == null)
			{
				selectionPanel.selectSingleVolProgrammatically(file.getFilenameStem());
				selectionPanel.removeFiles();
			}
			else
			{
				if (firstAddedFile == null)
				{
					firstAddedFile = file;
				}
			}
		}
		fileSelected(firstAddedFile);
		initializePanels();
		updatePanel();
        wizard.checkForOutputFiles();
        wizard.setSkipReanalyzeTraces(false);
	}

	public boolean hasDirectorySelection()
	{
		return true;
	}

	public void removeFiles(StsAbstractFile[] files)
	{
		wizard.setSkipReanalyzeTraces(true);
        for (int n = 0; n < files.length; n++)
        {
            StsFile file = (StsFile)files[n];
            StsSeismicBoundingBox volume = wizard.getSegyVolume(file.getFilename());
            wizard.removeSegyVolume(volume);
        }

        resetVolumeList();
		initializePanels();
		wizard.setSkipReanalyzeTraces(false);
		updatePanel();
	}

    public void removeVolume(StsSeismicBoundingBox volume)
    {
        selectionPanel.selectSingleVolProgrammatically(volume.getName());
        selectionPanel.removeFiles();
        updatePanel();
    }

    public void moveToAvailableList(StsSeismicBoundingBox volume)
    {
        selectionPanel.selectSingleVolProgrammatically(volume.getName());
        selectionPanel.moveToAvailableList();
        updatePanel();
    }

	public void clearFiles()
	{
		selectionPanel.removeSelectedFiles();
        updatePanel();
    }

	public void addAvailableVolume(StsSeismicBoundingBox volume)
	{
		selectionPanel.addAvailableFile(volume.getSegyFilename());
        updatePanel();
    }

    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }

    public void removeAllFiles()
	{
		wizard.setSkipReanalyzeTraces(true);
		wizard.clearSegyVolumesList();
		setPanelObject(null);
		initializePanels();
		updatePanel();
		wizard.setSkipReanalyzeTraces(false);
//		wizard.disableNext();
	}

	public void fileSelected(StsAbstractFile selectedFile)
	{
		if (selectedFile == null)
		{
			currentSelectedFile = null;
			currentVolume = null;
		}
		if (selectedFile != currentSelectedFile)
		{
			currentSelectedFile = (StsFile)selectedFile;
			currentVolume = wizard.getSegyVolume(currentSelectedFile.getFilename());
		}
	}

    public float getDiskRequired()
    {
        return calcDiskRequired();
    }

    private float calcDiskRequired()
    {
        StsSeismicBoundingBox[] volumes = wizard.getSegyVolumes();
        return StsSeismicBoundingBox.calcDiskRequired(volumes);
	}

    public void resetVolumeList()
    {
        StsSeismicBoundingBox[] volumes = wizard.getSegyVolumes();
        if (volumes == null || volumes.length == 0)
        {
            diskRequiredBean.setValue(0);
        }
        else
        {
            diskRequiredBean.setValue(wizard.calcDiskRequired());
        }
    }

    public void changeHeaders()
	{
		wizardStep.changeHeaders();
	}

	// This method is ultimately called by the trace analysis callback
	public void updatePanel()
	{
        updateBeans();

        if ((wizard.getSegyVolumesList().size() > 0) ||
			 (wizard.getSegyVolumesToProcess() != null) && (wizard.getSegyVolumesToProcess().length > 0))
		{
			wizard.enableNext();
		}
		else
		{
			wizard.disableNext();
		}
	}

	public String getVolumeName()
	{
		if (currentVolume == null)
			return "none";
		else
			return currentVolume.getName();
	}

	public boolean hasReloadButton()
	{
		return false;
	}

	public void setReload(boolean reload)
	{
	}

	public boolean getReload()
	{
		return true;
	}

	public boolean hasArchiveItButton()
	{
		return false;
	}

	public void setArchiveIt(boolean reload)
	{
	}

	public boolean getArchiveIt()
	{
		return true;
	}

	public boolean hasOverrideButton()
	{
		return false;
	}

	public void setOverrideFilter(boolean override)
	{
	}

	public boolean getOverrideFilter()
	{
		return false;
	}
}

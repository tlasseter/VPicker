package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSegyVspFileFormatPanel extends StsJPanel implements StsSelectRowNotifyListener, StsTableModelListener 
{
	// Going to hand build the global parameters panel in this class.
	private StsComboBoxFieldBean segyFormatComboBoxBean;
	//private StsDoubleFieldBean scanPercentBean;

	//private StsGroupBox segyFormatPanel = new StsGroupBox("SEGY Format");
//	private StsViewHeadersPanel viewHeadersPanel = null;
	private StsPoststackUnitsPanel unitsPanel = null;
	private StsFileFormatPanel fileFormatPanel = null;
	private StsTablePanelNew volumeStatusTablePanel = null;

	public StsProgressPanel progressPanel = StsProgressPanel.constructor(5, 50);

	// Z Domain
	public static final int TIME = 1;
	public static final int DEPTH = 2;

	private StsSeismicWizard wizard;
    private StsSegyVspFileFormat wizardStep;

    private StsModel model;

	public static final String NO_SEGY_FORMAT_ASSIGNED = "Default";

	private StsSegyVsp currentVolume = null;

	public StsSegyVspFileFormatPanel(StsSeismicWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
        this.wizardStep = (StsSegyVspFileFormat)wizardStep;
        model = wizard.getModel();
//		viewHeadersPanel = new StsViewHeadersPanel(wizard.frame);
		try
		{
			unitsPanel = new StsPoststackUnitsPanel(this.wizard);
			fileFormatPanel = new StsFileFormatPanel(wizard, true);

            buildTablePanel();
			constructBeans();
			constructPanel();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "dataMin", "dataMax", "statusString"};
        String[] columnTitles = {"Name", "Data Min", "Data Max", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Volumes");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }

	private void constructBeans()
	{
		segyFormatComboBoxBean = new StsComboBoxFieldBean(this, "segyFormatFilename", "SEGY Format:");
	}

	private void constructPanel() throws Exception
	{
		Color panelColor = Color.black;
		Font panelFont = new Font("Dialog", Font.PLAIN, 11);

		unitsPanel.setForeground(panelColor);
		unitsPanel.setFont(panelFont);

		fileFormatPanel.setForeground(panelColor);
		fileFormatPanel.setFont(panelFont);

		segyFormatComboBoxBean.setEditable(false);

		gbc.anchor = GridBagConstraints.NORTH;
		gbc.weighty = 1.0;
		gbc.fill = GridBagConstraints.BOTH;
		addEndRow(volumeStatusTablePanel);
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.weighty = 0.0;
		add(unitsPanel);
		add(fileFormatPanel);
//        add(viewHeadersPanel);
        gbc.weighty = 1.0;
        add(progressPanel);
        gbc.anchor = GridBagConstraints.SOUTH;
        add(progressPanel);

        volumeStatusTablePanel.addSelectRowNotifyListener( this);
        volumeStatusTablePanel.addTableModelListener( this);
	}

    public void rowsSelected( int[] selectedIndices)
    {
        try
        {
            wizard.setSkipReanalyzeTraces(true);
            wizard.setSelectedVolumes( selectedIndices);
            initializePanels();
        }
        finally
        {
            wizard.setSkipReanalyzeTraces(false);
        }
    }

	private boolean firstTime = true;

	public void initialize()
	{
		wizard.setSkipReanalyzeTraces(true);
		if (firstTime)
		{
			initializeSegyFormats();
			firstTime = false;
		}
		initializePanels();
		updatePanel();
		volumeStatusTablePanel.setSelectAll(true);

//		viewHeadersPanel.setVolumes(wizard.getSegyVolumes(), 0);
		wizard.setSkipReanalyzeTraces(false);
    }

	private void initializePanels()
	{
//		viewHeadersPanel.initialize();
		unitsPanel.initialize();
		fileFormatPanel.initialize(wizard.getVolumeSegyDatasets());
	}

	private void initializeSegyFormats()
	{
		segyFormatComboBoxBean.removeAll();

		String[] segyFormatFilenames = StsSEGYFormat.getSegyFormatFilenames();
		int nFiles = segyFormatFilenames.length;
		segyFormatComboBoxBean.addItem(NO_SEGY_FORMAT_ASSIGNED);
		for (int i = 0; i < segyFormatFilenames.length; i++)
		{
			segyFormatComboBoxBean.addItem(segyFormatFilenames[i]);
		}
		segyFormatComboBoxBean.setEditable(nFiles > 0);
	}

	public String getSegyFormatFilename()
	{
		StsSEGYFormat f = wizard.getSegyFormat();
		if (f == null)
		{
			return NO_SEGY_FORMAT_ASSIGNED;
		}
		return f.getName();
	}

	public void setSegyFormatFilename(String filename)
	{
		StsSEGYFormat currentFormat = wizard.getSegyFormat();
		if (currentFormat.getName().equalsIgnoreCase(filename))
		{
			return;
		}

		StsSEGYFormat segyFormat;
		if (filename == null || filename.equals("") || filename.equals(NO_SEGY_FORMAT_ASSIGNED))
		{
			segyFormat = StsSEGYFormat.constructor(model, StsSEGYFormat.POSTSTACK);
		}
		else
		{
			segyFormat = StsSEGYFormat.constructor(model, filename, StsSEGYFormat.POSTSTACK);
			if (segyFormat == null)
			{
				progressPanel.setDescriptionAndLevel("Invalid SEGY Format", StsProgressBar.ERROR);
				progressPanel.appendLine(filename + " is an invalid SEGY Format");
				initializePanels();
				return;
			}
		}
		wizard.setSegyFormat(segyFormat);
		initializePanels();
	}

    public StsTablePanelNew getVolumeStatusTablePanel() { return volumeStatusTablePanel; }

    public void updateFileFormatPanel()
    {
 //       fileFormatPanel.initialize(wizard.getSelectedSegyDatasets());
    }

    public void removeRows( int firstRow, int lastRow)
    {
        for( int i = firstRow; i <= lastRow; i++)
            wizard.moveSegyVolumeToAvailableList( wizard.getSegyVolume(i));

        if (wizard.getSegyVolumesList().size() >  0)
            wizard.enableNext();
        else
            wizard.disableNext();
    }

    // This method is ultimately called by the trace analysis callback
	public void updatePanel()
	{
        volumeStatusTablePanel.replaceRows( wizard.getSegyVolumesList());
	}

	public String getVolumeName()
	{
		if (currentVolume == null)
			return "none";
		else
			return currentVolume.getName();
	}

	final class SegyFilenameFilter implements FilenameFilter
    {
		String[] filter = null;
		int length;

		public SegyFilenameFilter(String[] filter)
		{
			this.filter = new String[filter.length];
			for (int i = 0; i < filter.length; i++)
			{
				this.filter[i] = filter[i].toLowerCase();
			}
		}

		public boolean accept(File dir, String name)
		{
			for (int i = 0; i < filter.length; i++)
			{
				if (name.toLowerCase().endsWith(filter[i]))
					return true;
			}
			return false;
		}

		public String getFilenameStem(String filename)
		{
			int filterStart = 1;
			for (int i = 0; i < filter.length; i++)
			{
				if (filename.toLowerCase().indexOf(filter[i]) > 0)
				{
					filterStart = filename.toLowerCase().indexOf(filter[i]);
					break;
				}
			}
			return filename.substring(0, filterStart - 1);
		}
	}
}

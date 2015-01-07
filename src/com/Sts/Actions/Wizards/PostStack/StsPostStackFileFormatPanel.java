package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
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

public class StsPostStackFileFormatPanel extends StsJPanel implements StsSelectRowNotifyListener, StsTableModelListener
{
	protected StsPoststackUnitsPanel unitsPanel = null;
	protected StsFileFormatPanel fileFormatPanel = null;
	public StsTablePanelNew volumeStatusTablePanel = null;
    public StsProgressPanel progressPanel = StsProgressPanel.constructor(5, 50);

    // Z Domain
	public static final int TIME = 1;
	public static final int DEPTH = 2;

	protected StsSeismicWizard wizard;
    private StsPostStackFileFormat wizardStep;
    private StsModel model;

	public StsPostStackFileFormatPanel(StsSeismicWizard wizard, StsPostStackFileFormat wizardStep)
	{
		this.wizard = wizard;
        this.wizardStep = wizardStep;
        model = wizard.getModel();
//		viewHeadersPanel = new StsViewHeadersPanel(wizard.frame);
		try
		{
            buildSubPanels();
			constructPanel();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

    protected void buildSubPanels()
    {
        unitsPanel = new StsPoststackUnitsPanel(wizard);
        fileFormatPanel = new StsFileFormatPanel(wizard, true);
        buildTablePanel();
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

	private void constructPanel() throws Exception
	{
		Color panelColor = Color.black;
		Font panelFont = new Font("Dialog", Font.PLAIN, 11);

		unitsPanel.setForeground(panelColor);
		unitsPanel.setFont(panelFont);

		fileFormatPanel.setForeground(panelColor);
		fileFormatPanel.setFont(panelFont);

//		segyFormatComboBoxBean.setEditable(false);

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
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.SOUTH;
        add(progressPanel);

        volumeStatusTablePanel.addSelectRowNotifyListener(this);
        volumeStatusTablePanel.addTableModelListener(this);
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

	public boolean initialize()
	{
        // setup volume status table
		updatePanel();
        if(!volumeStatusTablePanel.setSelectAll(true)) return false;
        // initialize fileFormat fields
        initializePanels();
//        viewHeadersPanel.setVolumes(wizard.getSegyVolumes(), 0);
        wizard.analyzeHeaders();
        return true;
    }

	private void initializePanels()
	{
//        viewHeadersPanel.initialize();
		unitsPanel.initialize();
        wizard.setSkipReanalyzeTraces(true);
        fileFormatPanel.initialize(wizard.getSelectedSegyDatasets());
        wizard.setSkipReanalyzeTraces(false);
    }
/*
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
			segyFormat = StsSEGYFormat.constructor(model, filename, StsSEGYFormat.POSTSTACK, false);
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
*/
    public StsTablePanelNew getVolumeStatusTablePanel() { return volumeStatusTablePanel; }

    public void updateFileFormatPanel()
    {
        fileFormatPanel.initialize(wizard.getSelectedSegyDatasets());
    }
    // This method is ultimately called by the trace analysis callback
	public void updatePanel()
	{
       volumeStatusTablePanel.replaceRows( wizard.getSegyVolumesList());
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

	final class SegyFilenameFilter implements FilenameFilter {
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

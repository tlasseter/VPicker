package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

import java.awt.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSegyVspBatchPanel extends StsJPanel implements StsTableModelListener
{
	private StsTablePanelNew volumesToProcessStatusTablePanel = null;

	// Z Domain
	public static final int TIME = 1;
	public static final int DEPTH = 2;

	private StsSegyVspWizard wizard;

	public static final String NO_SEGY_FORMAT_ASSIGNED = "Default";

	public StsSegyVspBatchPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = (StsSegyVspWizard)wizard;
		try
		{
            buildTablePanel();
            constructPanel();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "statusString"};
        String[] columnTitles = {"Name", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumesToProcessStatusTablePanel = new StsTablePanelNew(tableModel);
        volumesToProcessStatusTablePanel.setLabel("Volumes Awaiting Processing");
        volumesToProcessStatusTablePanel.setSize(400, 100);
        volumesToProcessStatusTablePanel.initialize();
    }

	private void constructPanel() throws Exception
	{
		gbc.fill = GridBagConstraints.BOTH;
		gbc.weighty = 1.0;
		addEndRow(volumesToProcessStatusTablePanel);

        volumesToProcessStatusTablePanel.addTableModelListener( this);
	}

    public void initialize()
    {
        wizard.setSkipReanalyzeTraces(true);
        updatePanel();
        gbc.fill = GridBagConstraints.BOTH;
        wizard.setSkipReanalyzeTraces(false);
    }

    public void updatePanel()
    {
        if (wizard.getSegyVolumesToProcess() != null &&
            wizard.getSegyVolumesToProcess().length > 0)
        {
            wizard.enableNext();
        }
        else
        {
            wizard.disableNext();
            wizard.gotoFirstStep();
        }

        volumesToProcessStatusTablePanel.replaceRows( wizard.getSegyVolumesToProcessList());
        volumesToProcessStatusTablePanel.repaint();
    }

    public void removeRows( int firstRow, int lastRow)
    {
        StsSeismicBoundingBox[] segyVolumesToProcess = wizard.getSegyVolumesToProcess();
        if( segyVolumesToProcess == null) return;
        if( firstRow >= segyVolumesToProcess.length || lastRow >= segyVolumesToProcess.length) return;
        for( int i = firstRow; i <= lastRow; i++)
            wizard.removeSegyVolumeToProcess( segyVolumesToProcess[i]);
    }
}

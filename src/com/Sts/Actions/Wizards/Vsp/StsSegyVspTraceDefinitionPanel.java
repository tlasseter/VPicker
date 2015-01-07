package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.PostStack3d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSegyVspTraceDefinitionPanel extends StsPostStackTraceDefinitionPanel implements ChangeListener, ListSelectionListener, StsSelectRowNotifyListener, StsTableModelListener
{
    public StsSegyVspTraceDefinitionPanel(StsSeismicWizard wizard, StsWizardStep wizardStep)
	{
        super(wizard, wizardStep);
		try
		{
			overrideGeometryBean = new StsBooleanFieldBean(this, "overrideHeader", false, "Override Headers", true);
			showOnlyReqAttributesBean = new StsBooleanFieldBean(this, "selectedAttributes", false, "Only Show Required Attributes", false);
            buildTablePanel();
			constructPanel();
		}
		catch (Exception e)
		{
			StsException.systemError(this, "constructor", e.getMessage());
		}
	}

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "xOrigin", "yOrigin", "rowNumMin", "rowNumMax", "colNumMin", "colNumMax", "statusString"};
        String[] columnTitles = {"Name", "Start X", "Start Y", "Start Inline", "End Inline", "Start XLine", "End XLine", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Volumes");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }

	protected void constructPanel()
	{
		segyTable.setFont(new Font("Dialog", 3, 12));

		optionsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		optionsPanel.gbc.anchor = gbc.WEST;
		optionsPanel.addToRow(showOnlyReqAttributesBean);
		optionsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		if (wizardStep.getOverrideHeaderAllowed())
		{
			optionsPanel.addToRow(overrideGeometryBean);
		}
		attributePanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		attributePanel.gbc.anchor = gbc.WEST;
		attributePanel.addToRow(attNameLabel, 1, 0.0);
		attributePanel.addToRow(attNameTxt, 2, 1.0);
		attributePanel.addToRow(typeLabel, 1, 0.0);
		attributePanel.addToRow(headerFormatCombo, 3, 0.0);
		attributePanel.addToRow(byteLabel, 1, 0.0);
		attributePanel.addToRow(traceHdrPosTxt, 1, 0.0);
		attributePanel.gbc.anchor = gbc.EAST;
		attributePanel.addEndRow(addBtn, 1, 0.0);

		attributePanel.gbc.anchor = gbc.EAST;
		attributePanel.addToRow(descriptionLabel, 1, 0.0);
		attributePanel.gbc.anchor = gbc.WEST;
		attributePanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		attributePanel.addToRow(attDescTxt, 6, 1.0);
		attributePanel.gbc.anchor = gbc.EAST;
		attributePanel.addToRow(applyScalarLabel);
		attributePanel.gbc.anchor = gbc.WEST;
		attributePanel.addToRow(applyScalarBean, 1, 0.0);
		attributePanel.gbc.anchor = gbc.EAST;
		attributePanel.addEndRow(updateBtn, 1, 0.0);

		tableScrollPane1.setAutoscrolls(true);
		tableScrollPane1.getViewport().add(segyTable, null);
		tableScrollPane1.getViewport().setPreferredSize(new Dimension(550, 100));

		sliderPanel.setIncrementLabel("Step");
		sliderPanel.setTextColor(Color.black);
		sliderPanel.setValueLabel("Seq Trace:");
		sliderPanel.setSelected(true);
		sliderPanel.getCheckBoxSlider().setVisible(false);
		sliderPanel.addChangeListener(this);

		tracesPanel.gbc.anchor = gbc.NORTH;
		tracesPanel.gbc.fill = GridBagConstraints.BOTH;
		tracesPanel.add(tableScrollPane1, 1, 1.0);
		tracesPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		tracesPanel.add(sliderPanel, 1, 0.0);

		// now put it all together
		attributesPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		attributesPanel.gbc.anchor = gbc.NORTH;
		attributesPanel.add(optionsPanel, 1, 0.0);
		attributesPanel.add(attributePanel, 1, 0.0);
		attributesPanel.gbc.fill = GridBagConstraints.BOTH;
		attributesPanel.add(tracesPanel, 1, 1.0);

        gbc.fill = gbc.BOTH;
        add( volumeStatusTablePanel, 1, 1);
        attributesPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        add( attributesPanel, 1, 0);

        gbc.anchor = GridBagConstraints.SOUTH;
        add(progressPanel, 1, 1);

        segyTable.addListSelectionListener(this);
        volumeStatusTablePanel.addSelectRowNotifyListener(this);
        volumeStatusTablePanel.addTableModelListener(this);
	}

	public boolean initialize()
	{
		wizard.setSkipReanalyzeTraces(true);
        StsSEGYFormat segyFormat = wizard.getSegyFormat();
		currentColNum = 0;

        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if(selectedVolume != null)
        {
            traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(selectedVolume));
            setupTraceRecords(selectedVolume);
        }

		segyTable.setPreferredSize(new Dimension(25 * 75, 1500));

        headerFormatCombo.removeAllItems();
		String[] traceFormatStrings = StsSEGYFormat.headerFormatStrings;
		for (int i = 0; i < traceFormatStrings.length; i++)
			headerFormatCombo.addItem(traceFormatStrings[i]);

        applyScalarBean.removeAllItems();
		String[] scalarStrings = StsSEGYFormat.scalarAttributes;
		for (int i = 0; i < scalarStrings.length; i++)
			applyScalarBean.addItem(scalarStrings[i]);
		setPanelObject(segyFormat);

        volumeStatusTablePanel.replaceRows(wizard.getSegyVolumesList());
        volumeStatusTablePanel.setSelectionIndex(0);
        volumesSelected( wizard.getSelectedVolumes());

		wizard.setSkipReanalyzeTraces(false);
        wizard.rebuild();
        return true;
    }

	public void volumesSelected(StsSeismicBoundingBox[] volumes)
	{
		wizard.saveSelectedAttributes(wizard.getSelectedVolumes()[0]);
		wizard.setSelectedVolumes(volumes);
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if(selectedVolume != null)
        {
            traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(selectedVolume));
            setupTraceRecords(selectedVolume);
        }
		sliderPanel.setValue(0);
		setValues();
		initSlider(selectedVolume);
		valueChanged(null);
	}

	public void initSlider(StsSeismicBoundingBox selectedVolume)
	{
		int nTotalTraces = wizardStep.getNTraces(wizard.getSelectedVolumes()[0]);
		int[] startAndCount = wizardStep.getTraceStartAndCount(wizard.getSelectedVolumes()[0], 0);
		sliderPanel.initSliderValues(0.0f, nTotalTraces - 1, 1.0f, 0.0f);
		sliderPanel.setIncrementValue(startAndCount[1]);
		sliderPanel.setIncrement(startAndCount[1]);
	}

	public void setOverrideHeader(boolean override)
	{
		this.overrideGeometry = override;

		if (wizardStep.getOverrideHeaderAllowed() && override)
		{
			// Show Known Point Dialog
			JDialog dialog = new JDialog(((StsWizard)wizard).frame);
			StsSegySurveyDefinitionPanel panel = new StsSegySurveyDefinitionPanel(wizard, wizardStep, dialog);
			dialog.setSize(400, 500);
			dialog.setTitle("SegY Survey Definition");
			dialog.getContentPane().add(panel);
			dialog.setModal(true);
			dialog.setVisible(true);

			wizardStep.analyzeGrid();
		}
	}

    public void rowsSelected( int[] selectedIndices)
    {
        wizard.setSelectedVolumes( selectedIndices);
        volumesSelected( wizard.getSelectedVolumes());
    }

	public void valueChanged(ListSelectionEvent e)
	{
		if (segyTable.getNumberOfRows() > 0)
		{
			if (segyTable.getSelectedIndices().length == 0)
			{
				enableAttributeFields(false);
				return;
			}
			enableAttributeFields(true);
			attIndices = segyTable.getSelectedIndices();

			StsSEGYFormatRec record = null;
			if (showOnlyReqAttributesBean.isSelected())
			{
				record = requiredRecords[attIndices[0]];
			}
			else
			{
				record = allRecords[attIndices[0]];
			}
			attNameTxt.setText(record.getUserName());
			if(record.getRequired())
			   attNameTxt.setEnabled(false);
			attDescTxt.setText(record.getDescription());
			headerFormatCombo.setSelectedItem(StsSEGYFormat.sampleFormatStrings[record.format]);
			traceHdrPosTxt.setValue(record.getLoc() + 1);
			applyScalarBean.setSelectedItem(record.getApplyScalar());
		}
	}

	/**
	* Called when the "update" button is clicked
	*/
	public void update()
	{
		StsSEGYFormatRec record = null;
		if (showOnlyReqAttributesBean.isSelected())
		{
			record = requiredRecords[attIndices[0]];
		}
		else
		{
			record = allRecords[attIndices[0]];
		}

		record.setUserName(attNameTxt.getText());
		record.setDescription(attDescTxt.getText());
		int format = StsSEGYFormat.getSampleFormatFromString((String)headerFormatCombo.getSelectedItem());
		record.setFormat(format);
		record.setLoc(traceHdrPosTxt.getIntValue() - 1);
		int selected = applyScalarBean.getSelectedIndex();
		record.setApplyScalar(StsSEGYFormat.scalarAttributes[selected]);
		StsSEGYFormat segyFormat = wizard.getSegyFormat();
		if(segyFormat != null)
			segyFormat.setFormatChanged(true);

		wizardStep.analyzeGrid();
	}

	public void updatePanel()
	{
        volumeStatusTablePanel.replaceRows(wizard.getSegyVolumesList());
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if(selectedVolume != null)
            setupTraceRecords(selectedVolume);
        setValues();
        if (wizard.getVolumes().length > 0)
            wizardStep.getWizard().enableNext();
        else
            wizardStep.getWizard().disableNext();
	}

	private boolean ignoreChangeEvent = false;

	public boolean getSelectedAttributes()
	{
		return selectedAttributes;
	}

	public boolean getOverrideHeader()
	{
		return this.overrideGeometry;
	}

    public void stateChanged(ChangeEvent e)
	{
		StsSliderBean source = (StsSliderBean)e.getSource();
		if (source == sliderPanel)
		{
			if (ignoreChangeEvent)
			{
				return;
			}
			if (!source.isDraggingSlider() && !source.isInitializing())
			{
				currentColNum = (int)sliderPanel.getValue();
				setValues();
			}
		}
		return;
	}

    static public void main(String[] args)
	{
		JFrame frame = new JFrame("Test Panel");
		frame.setSize(300, 200);
//        frame.addKeyListener();

		Container contentPane = frame.getContentPane();
		StsSegyVspTraceDefinitionPanel panel = new StsSegyVspTraceDefinitionPanel(null, null);
		contentPane.add(panel);
		StsToolkit.centerComponentOnScreen(frame);
		frame.pack();
		frame.setVisible(true);
	}
}

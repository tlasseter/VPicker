package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.TraceDefinition.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

abstract public class StsPostStackTraceDefinitionPanel extends StsFieldBeanPanel implements ChangeListener, ListSelectionListener,
                                                                                StsSelectRowNotifyListener, StsTableModelListener
{
    abstract public void buildTablePanel();

    protected StsSeismicWizard wizard;
    protected StsPostStackTraceDefinition wizardStep;
	protected int currentColNum = -1;
    protected StsSEGYFormatRec[] allRecords = null;

	protected StsSEGYFormatRec[] requiredRecords;
	protected int[] attIndices = null;
	protected StsJPanel attributesPanel = StsJPanel.addInsets();

	// options panel and it's components
	protected StsGroupBox optionsPanel = new StsGroupBox("Options");
	protected StsBooleanFieldBean showOnlyReqAttributesBean = null;
    protected StsBooleanFieldBean overrideGeometryBean = null;

    // attribute panel and it's components
	protected StsGroupBox attributePanel = new StsGroupBox("Attribute");
	protected JLabel attNameLabel = new JLabel("Name:");
	protected JTextField attNameTxt = new JTextField();
	protected JLabel typeLabel = new JLabel("Type:");
	protected JComboBox headerFormatCombo = new JComboBox();
	protected JLabel byteLabel = new JLabel("Byte:");
	protected StsIntFieldBean traceHdrPosTxt = new StsIntFieldBean(true);
	protected StsToggleButton addBtn = new StsToggleButton("add/remove", "Add/remove record from required collection.", this, "add", "remove");
	protected JLabel descriptionLabel = new JLabel("Description:");
	protected JTextField attDescTxt = new JTextField();
	protected JLabel applyScalarLabel = new JLabel("Apply Scalar:");
    protected StsFloatFieldBean scalarValue = new StsFloatFieldBean(false);
	protected JComboBox applyScalarBean = new JComboBox();
	protected StsButton updateBtn = new StsButton("Update", "Updates fields for this record.", this, "update");
	public StsTablePanelNew volumeStatusTablePanel = null;

	// traces panel
	protected StsGroupBox tracesPanel = new StsGroupBox("Traces");
	protected StsTablePanel segyTable = new StsTablePanel(false);
	protected JScrollPane tableScrollPane1 = new JScrollPane();

	// slider panel
	protected StsSliderBean sliderPanel = new StsSliderBean(true, false);
	protected StsProgressPanel progressPanel = StsProgressPanel.constructor(5, 25);

	protected int headerSize = StsSEGYFormat.defaultTraceHeaderSize;
	protected boolean overrideGeometry = false;
	protected boolean applyScalar = true;
	protected boolean selectedAttributes = false;

    static public final boolean debug = false;

    public StsPostStackTraceDefinitionPanel(StsSeismicWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = (StsPostStackTraceDefinition)wizardStep;

		try
		{
			overrideGeometryBean = new StsBooleanFieldBean(wizard, "overrideGeometry", false, "Override geometry (use survey definition step)", false);
			showOnlyReqAttributesBean = new StsBooleanFieldBean(this, "selectedAttributes", false, "Only Show Required Attributes", false);
			buildTablePanel();
            constructPanel();
		}
		catch (Exception e)
		{
			StsException.systemError(this, "constructor", e.getMessage());
		}
	}

    protected void constructPanel()
    {
        segyTable.setFont(new Font("Dialog", 3, 12));

        optionsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        optionsPanel.gbc.anchor = gbc.WEST;
        optionsPanel.add(showOnlyReqAttributesBean);
        optionsPanel.add(overrideGeometryBean);

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
        attributePanel.addToRow(scalarValue);
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
        // StsSEGYFormat segyFormat = wizard.getSegyFormat();
		currentColNum = 0;
        //TODO getTraceHeaderSize directory from segyFormat which is available directly from wizard rather than thru segyVolume
        ArrayList segyVolumesList = wizard.getSegyVolumesList();
         if(segyVolumesList.size() == 0) return false;
        StsSeismicBoundingBox selectedVolume = (StsSeismicBoundingBox)segyVolumesList.get(0);
        traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(selectedVolume));
        setupTraceRecords(selectedVolume);

        segyTable.setPreferredSize(new Dimension(25 * 75, 1500));

        headerFormatCombo.removeAllItems();
		String[] traceFormatStrings = StsSEGYFormat.headerFormatStrings;
		for (int i = 0; i < traceFormatStrings.length; i++)
			headerFormatCombo.addItem(traceFormatStrings[i]);

        applyScalarBean.removeAllItems();
		String[] scalarStrings = StsSEGYFormat.scalarAttributes;
		for (int i = 0; i < scalarStrings.length; i++)
			applyScalarBean.addItem(scalarStrings[i]);

        volumeStatusTablePanel.replaceRows(wizard.getSegyVolumesList());
        //volumeStatusTablePanel.setSelectionIndex(0);
        if(!volumeStatusTablePanel.setSelectAll(true)) return false;
        volumesSelected(wizard.getSelectedVolumes());

		wizard.setSkipReanalyzeTraces(false);
        wizard.rebuild();
        return true;
    }

	public void volumesSelected(StsSeismicBoundingBox[] volumes)
	{
        if( volumes == null || volumes.length == 0) return;
        StsSeismicBoundingBox selectedVolume = volumes[0];
        traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(selectedVolume));
		setupTraceRecords(selectedVolume);
		sliderPanel.setValue(0);
		setValues();
		initSlider(selectedVolume);
		valueChanged(null);
	}

	public void setupTraceRecords(StsSeismicBoundingBox selectedVolume)
	{
        if( selectedVolume == null) return;
        StsSEGYFormat segyFormat = selectedVolume.getSegyFormat();
		allRecords = segyFormat.getAllTraceRecords();
		requiredRecords = segyFormat.getRequiredTraceRecords();
        setPanelObject(segyFormat);
    }

	public void updateRowStatus()
	{
		if (!showOnlyReqAttributesBean.isSelected())
		{
			StsSEGYFormatRec record = null;
			for (int i = 0; i < allRecords.length; i++)
			{
				segyTable.setRowType(i, StsTablePanelNew.NOT_HIGHLIGHTED);
				record = allRecords[i];
				if (record.isRequired())
					segyTable.setRowType(i, StsTablePanelNew.NOT_EDITABLE);
			}
		}
		else
		{
			for (int i = 0; i < requiredRecords.length; i++)
			    segyTable.setRowType(i, StsTablePanelNew.NOT_EDITABLE);
		}
	}

	public void initSlider(StsSeismicBoundingBox selectedVolume)
	{
        if (selectedVolume == null) return;
        int nTotalTraces = wizardStep.getNTraces(selectedVolume);
		int[] startAndCount = wizardStep.getTraceStartAndCount(selectedVolume, 0);
		sliderPanel.initSliderValues(0.0f, nTotalTraces - 1, 1.0f, 0.0f);
		sliderPanel.setIncrementValue(startAndCount[1]);
		sliderPanel.setIncrement(startAndCount[1]);
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
            setScalarValue( record.getApplyScalar());
            if(record.userRequired)
            {
                addBtn.setSelected(true);
                addBtn.setEnabled(true);
            }
            else if(!record.required)
            {
                addBtn.setSelected(false);
                addBtn.setEnabled(true);
            }
            else
                addBtn.setEnabled(false);
        }
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

	protected void enableAttributeFields(boolean enabled)
	{
		attNameTxt.setEnabled(enabled);
		headerFormatCombo.setEnabled(enabled);
		traceHdrPosTxt.setEnabled(enabled);
		addBtn.setEnabled(enabled);
		applyScalarBean.setEnabled(enabled);
		applyScalarBean.setEditable(enabled);
		updateBtn.setEnabled(enabled);
		attDescTxt.setEnabled(enabled);
		traceHdrPosTxt.setEditable(enabled);
	}

    public void resetSelectedVolumeStatus()
    {
        wizard.resetSelectedTraceStatus();
    }

    public void add()
	{
		StsSEGYFormatRec record = null;
		if (showOnlyReqAttributesBean.isSelected()) return;

        // Adding Attribute
        for (int i = 0; i < attIndices.length; i++)
        {
            record = allRecords[attIndices[i]];
            if(!record.required && !record.userRequired)
            {
                requiredRecords = (StsSEGYFormatRec[])StsMath.arrayAddElement(requiredRecords, record);
                record.userRequired = true;
            }
        }

        record.setUserName(attNameTxt.getText());
        record.setDescription(attDescTxt.getText());
        int format = StsSEGYFormat.getSampleFormatFromString((String)headerFormatCombo.getSelectedItem());
        record.setFormat(format);
        record.setLoc(traceHdrPosTxt.getIntValue() - 1);
        int selected = applyScalarBean.getSelectedIndex();
        record.setApplyScalar(StsSEGYFormat.scalarAttributes[selected]);
	}

    public void remove()
	{
		StsSEGYFormatRec record = null;
		if (showOnlyReqAttributesBean.isSelected())
		{
			// Removing attribute

			// doing in reverse order so we don't need to correct the index as we delete
			for (int i = attIndices.length - 1; i >= 0 ; i--)
			{
                int index = attIndices[i];
                record = requiredRecords[index];
				if(!record.required && record.userRequired)
				{
					requiredRecords = (StsSEGYFormatRec[])StsMath.arrayDeleteElement(requiredRecords, index);
					record.userRequired = false;
				}
			}
		}
		else
		{
			// Adding Attribute
			for (int i = 0; i < attIndices.length; i++)
			{
				record = allRecords[attIndices[i]];
                if(!record.required && record.userRequired)
				{
					record.userRequired = false;
				}

			}
        }

        record.setUserName(attNameTxt.getText());
        record.setDescription(attDescTxt.getText());
        int format = StsSEGYFormat.getSampleFormatFromString((String)headerFormatCombo.getSelectedItem());
        record.setFormat(format);
        record.setLoc(traceHdrPosTxt.getIntValue() - 1);
        int selected = applyScalarBean.getSelectedIndex();
        record.setApplyScalar(StsSEGYFormat.scalarAttributes[selected]);
	}
    /**
	* Called when the "update" button is clicked
	*/
	public void update()
	{
		StsSEGYFormatRec record;
		if (showOnlyReqAttributesBean.isSelected())
		{
			record = requiredRecords[attIndices[0]];
		}
		else
		{
			record = allRecords[attIndices[0]];
		}

        boolean changed = false;
        record.setUserName(attNameTxt.getText());
		record.setDescription(attDescTxt.getText());
		int format = StsSEGYFormat.getSampleFormatFromString((String)headerFormatCombo.getSelectedItem());
		changed = changed | record.setFormat(format);
		changed = changed | record.setLoc(traceHdrPosTxt.getIntValue() - 1);
        int selected = applyScalarBean.getSelectedIndex();
        changed = changed | record.setApplyScalarSelectedIndex(selected);
        if(!changed) return;
		wizard.getSegyFormat().setFormatChanged(true);
        wizard.resetSelectedTraceStatus();
        reanalyze();
    }

	public void updatePanel()
	{
        if(debug) StsException.systemDebug(this, "updatePanel", " time:" + System.currentTimeMillis());
        volumeStatusTablePanel.replaceRows(wizard.getSegyVolumesList());
        setupTraceRecords(wizard.getFirstSelectedVolume());
        setValues();
        segyTable.repaint();
    }

    public void repaint()
    {
        if(debug) StsException.systemDebug(this, "repaint", " time:" + System.currentTimeMillis());
        super.repaint();
    }

    private boolean ignoreChangeEvent = false;

	public boolean setValues()
	{
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if(selectedVolume == null) return false;
        if(allRecords == null) return false;

		wizardStep.initSetValues(selectedVolume);

		segyTable.removeListSelectionListener(this);

		int[] startAndCount = wizardStep.getTraceStartAndCount(wizard.getSelectedVolumes()[0], currentColNum);
		int nCols  = startAndCount[1];

		if (nCols < 0)
		{
			return false;
		}

		ignoreChangeEvent = true;
		if (! wizardStep.getOverrideStep())
		{
			sliderPanel.setIncrementValue(nCols);
			sliderPanel.setIncrement(nCols);
		}
		sliderPanel.setValue(startAndCount[0]);
		ignoreChangeEvent = false;

		boolean littleEndian = selectedVolume.getIsLittleEndian();

		Object[] data = null;
		if (showOnlyReqAttributesBean.isSelected())
		{
			data = requiredRecords;
		}
		else
		{
			data = allRecords;
		}

		int xyScalePos = 0;
		int xyScaleFmt = 0;
		int edScalePos = 0;
		int edScaleFmt = 0;
		StsSEGYFormat segyFormat = selectedVolume.getSegyFormat();
		if (segyFormat != null)
		{
			StsSEGYFormatRec rec = segyFormat.getTraceRec("CO-SCAL");
			xyScalePos = rec.getLoc();
			xyScaleFmt = rec.getFormat();
			rec = segyFormat.getTraceRec("ED-SCAL");
			edScalePos = rec.getLoc();
			edScaleFmt = rec.getFormat();
		}

		byte[][] traceBinaryHeader = new byte[nCols][];
		int trace = startAndCount[0];
		int startTrace = trace;
		for (int col = 0; col < nCols; col++, trace++)
		{
			traceBinaryHeader[col] = wizardStep.getTraceHeaderBinary(selectedVolume, trace);
			if (traceBinaryHeader[col] == null)
			{
				nCols = col;
				break;
			}
		}

		StsSegyTraceTableModel model = new StsSegyTraceTableModel(data,
																  startTrace,
																  nCols,
																  xyScalePos,
																  xyScaleFmt,
																  edScalePos,
																  edScaleFmt,
																  traceBinaryHeader,
																  littleEndian);
		segyTable.setPreferredSize(new Dimension(nCols*75,1500));
		segyTable.setTableModel(model);
		segyTable.revalidate();	
		model.fireTableDataChanged();
		updateRowStatus();
		segyTable.addListSelectionListener(this);
		return true;
	}

	public boolean getSelectedAttributes()
	{
		return selectedAttributes;
	}

	public void setSelectedAttributes(boolean val)
	{
		selectedAttributes = val;
		if (selectedAttributes)
		{
			addBtn.setText("Remove");
		}
		else
		{
			addBtn.setText("Add");
		}
		setValues();
		valueChanged(null);
		this.tableScrollPane1.getViewport().setViewPosition(new Point(0, 0));
	}

	public boolean getOverrideGeometry()
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

	public int getRequiredRowCount()
	{
		return requiredRecords.length;
	}

	public StsSEGYFormatRec getRequiredRecord(int i)
	{
		return requiredRecords[i];
	}

    public void setScalarValue( String applyScalar)
    {
        float scalar = 0.0f;
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if(selectedVolume == null) return;
        StsSEGYFormat segyFormat = selectedVolume.getSegyFormat();
        boolean littleEndian = selectedVolume.getIsLittleEndian();
        byte[] traceBinaryHeader = wizardStep.getTraceHeaderBinary(selectedVolume, 0);
        if( applyScalar.equals("CO-SCAL"))
            scalar = segyFormat.getFloatCoordinateScale( traceBinaryHeader, littleEndian);
        else if( applyScalar.equals("ED-SCAL"))
            scalar = segyFormat.getFloatElevationScale( traceBinaryHeader, littleEndian);
        scalarValue.setValue( scalar);
    }

    public void reanalyze()
    {
        wizardStep.analyzeGrid();
    }
}

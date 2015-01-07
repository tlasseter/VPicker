package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack3d.*;
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
import java.util.List;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSegyTraceDefinitionPanel extends StsFieldBeanPanel implements ChangeListener, ListSelectionListener,
                                                                    StsSelectRowNotifyListener, StsTableModelListener
{
	private StsPreStackWizard wizard;
	private StsSegyPreStackTraceDefinition wizardStep;
	private int currentColNum = -1;
	public StsSEGYFormatRec[] fmtRecs = null;

	public List requiredRecords = new LinkedList();
	private int[] attIndices = null;

	private StsJPanel attributesPanel = StsJPanel.addInsets();

	// options panel and it's components
	private StsGroupBox optionsPanel = new StsGroupBox("Options");
	private StsBooleanFieldBean overrideCoorBean = null;
	private StsBooleanFieldBean selectedTableChk = null;

	// attribute panel and it's components
	private StsGroupBox attributePanel = new StsGroupBox("Attribute");
	private JLabel attNameLabel = new JLabel("Name:");
	private JTextField attNameTxt = new JTextField();
	private JLabel typeLabel = new JLabel("Type:");
	private JComboBox headerFormatCombo = new JComboBox();
	private JLabel byteLabel = new JLabel("Byte:");
	private StsIntFieldBean traceHdrPosTxt = new StsIntFieldBean(true);
	private StsButton addBtn = new StsButton("Add", "Adds field to this record.", this, "add");
	private JLabel descriptionLabel = new JLabel("Description:");
	private JTextField attDescTxt = new JTextField();
	private JLabel applyScalarLabel = new JLabel("Apply Scalar:");
	private JComboBox applyScalarBean = new JComboBox();
	private StsButton updateBtn = new StsButton("Update", "Updates fields for this record.", this, "update");
    protected StsTablePanelNew volumeStatusTablePanel = null;

	// traces panel
	private StsGroupBox tracesPanel = new StsGroupBox("Traces");
	private StsTablePanel segyTable = new StsTablePanel(false);
	private JScrollPane tableScrollPane1 = new JScrollPane();

	// slider panel
	private StsSliderBean sliderPanel = new StsSliderBean(true, false);

	public StsProgressPanel progressPanel = StsProgressPanel.constructor(5, 50);

	private JSplitPane jSplitPane1 = new JSplitPane();
	private StsJPanel bottomPanel = StsJPanel.addInsets();

	protected int headerSize = StsSEGYFormat.defaultTraceHeaderSize;
	protected boolean overrideCoor = false;
	protected boolean applyScalar = true;
	private boolean selectedAttributes = false;

	public StsSegyTraceDefinitionPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = (StsPreStackWizard)wizard;
		this.wizardStep = (StsSegyPreStackTraceDefinition)wizardStep;

		try
		{
			overrideCoorBean = new StsBooleanFieldBean(this, "overrideHeader", false, "Override Headers", true);
			selectedTableChk = new StsBooleanFieldBean(this, "selectedAttributes", false, "Only Show Required Attributes", false);
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

    private void constructPanel()
    {
        segyTable.setFont(new Font("Dialog", 3, 12));

        optionsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        optionsPanel.gbc.anchor = gbc.WEST;
        optionsPanel.addToRow(selectedTableChk);
        optionsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        if (wizardStep.getOverrideHeaderAllowed())
        {
            optionsPanel.addToRow(overrideCoorBean);
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

    public void initialize()
    {
        wizard.setSkipReanalyzeTraces(true);
        StsSEGYFormat segyFormat = wizard.getSegyFormat();
        currentColNum = 0;

        ArrayList segyVolumesList = wizard.getSegyVolumesList();
        if( segyVolumesList != null && segyVolumesList.size() > 0)
        {
            StsSeismicBoundingBox segyVolume = (StsSeismicBoundingBox)wizard.getSegyVolumesList().get(0);
            traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(segyVolume));
            setupTraceRecords();
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
    }

    public void volumesSelected(StsSeismicBoundingBox[] volumes)
	{
        if( volumes == null) return;
        StsSeismicBoundingBox segyVolume = wizard.getFirstSegyVolume();
        if( segyVolume == null) return;

        wizard.saveSelectedAttributes( segyVolume);
		wizard.setSelectedVolumes(volumes);

		traceHdrPosTxt.setValueAndRange(0, 0, wizardStep.getTraceHeaderSize(volumes[0]));

		setupTraceRecords();
		sliderPanel.setValue(0);
		setValues();
		initSlider();
		valueChanged(null);
	}

	public void setupTraceRecords()
	{
        StsSeismicBoundingBox segyVolume = wizard.getFirstSegyVolume();
        if( segyVolume == null) return;

        fmtRecs = wizardStep.getAllTraceRecords( segyVolume);
		StsSEGYFormatRec[] records = wizardStep.getRequiredTraceRecords( segyVolume);
		requiredRecords.clear();
		StsSEGYFormat segyFormat = wizard.getSegyFormat();
		if((records != null) && (segyFormat != null))
		{
			for (int i = 0; i < records.length; i++)
			{
				requiredRecords.add(fmtRecs[segyFormat.getTraceRecIndex(records[i].getName())]);
			}
		}
		else
		{
			for (int i = 0; i < fmtRecs.length; i++)
			{
                            if (fmtRecs[i].required || fmtRecs[i].userRequired)
					requiredRecords.add(fmtRecs[i]);
			}
		}
	}

	public void updateRowStatus()
	{
		if (!selectedTableChk.isSelected())
		{
			StsSEGYFormatRec record = null;
			for (int i = 0; i < fmtRecs.length; i++)
			{
				segyTable.setRowType(i, StsTablePanelNew.NOT_HIGHLIGHTED);
				record = fmtRecs[i];
				if (requiredRecords.contains(record))
				{
					segyTable.setRowType(i, StsTablePanelNew.NOT_EDITABLE);
				}
			}
		}
		else
		{
			for (int i = 0; i < requiredRecords.size(); i++)
			{
                                StsSEGYFormatRec rec = ((StsSEGYFormatRec)requiredRecords.get(i));
                                if (rec.required || rec.userRequired)
				{
					segyTable.setRowType(i, StsTablePanelNew.NOT_EDITABLE);
				}
				else
				{
					segyTable.setRowType(i, StsTablePanelNew.NOT_HIGHLIGHTED);
				}
			}
		}
	}

	public void initSlider()
	{
        StsSeismicBoundingBox segyVolume = wizard.getFirstSegyVolume();
        if( segyVolume == null) return;

        int nTotalTraces = wizardStep.getNTraces( segyVolume);
		int[] startAndCount = wizardStep.getTraceStartAndCount( segyVolume, 0);
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
			if (selectedTableChk.isSelected())
			{
				record = (StsSEGYFormatRec)requiredRecords.get(attIndices[0]);
			}
			else
			{
				record = fmtRecs[attIndices[0]];
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

    public void removeRows( int firstRow, int lastRow)
    {
        for( int i = firstRow; i <= lastRow; i++)
            wizard.moveSegyVolumeToAvailableList( wizard.getSegyVolume(i));

        if (wizard.getSegyVolumesList().size() >  0)
            wizard.enableNext();
        else
            wizard.disableNext();
    }

    public void setOverrideHeader(boolean override)
    {
        this.overrideCoor = override;

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

            wizardStep.analyzeGeometry();
        }
    }

	private void enableAttributeFields(boolean enabled)
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

	public void add()
	{
		StsSEGYFormatRec record = null;
		if (selectedTableChk.isSelected())
		{
			// Removing attribute

			// doing in reverse order so we don't need to correct the index as we delete
			for (int i = attIndices.length - 1; i >= 0 ; i--)
			{
				record = (StsSEGYFormatRec)requiredRecords.get(attIndices[i]);
				if (record.required)
				{
					new StsMessage(((StsWizard)wizard).frame, StsMessage.WARNING, "Cannot remove "
								   + record.getUserName() + " - it is a required attribute.");
				}
				else
				{
					requiredRecords.remove(record);
					record.userRequired = false;
				}
			}

		}
		else
		{
			// Adding Attribute
			boolean found = false;
			for (int i = 0; i < attIndices.length; i++)
			{
				record = fmtRecs[attIndices[i]];
				record.userRequired = true;
				found = requiredRecords.contains(record);
				if (!found)
				{
					requiredRecords.add(record);
				}
			}

			fmtRecs[attIndices[0]].setUserName(attNameTxt.getText());
			fmtRecs[attIndices[0]].setDescription(attDescTxt.getText());
			int format = StsSEGYFormat.getSampleFormatFromString((String)headerFormatCombo.getSelectedItem());
			fmtRecs[attIndices[0]].setFormat(format);
			fmtRecs[attIndices[0]].setLoc(traceHdrPosTxt.getIntValue() - 1);
			int selected = applyScalarBean.getSelectedIndex();
			fmtRecs[attIndices[0]].setApplyScalar(StsSEGYFormat.scalarAttributes[selected]);
		}
		wizard.getSegyFormat().setFormatChanged(true);

        wizardStep.analyzeGeometry();

		setValues();
	}

	/**
	* Called when the "update" button is clicked
	*/
	public void update()
	{
		StsSEGYFormatRec record = null;
		if (selectedTableChk.isSelected())
		{
			record = (StsSEGYFormatRec)requiredRecords.get(attIndices[0]);
		}
		else
		{
			record = fmtRecs[attIndices[0]];
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


		wizardStep.analyzeGeometry();
	}

    public void updatePanel()
    {
        volumeStatusTablePanel.replaceRows(wizard.getSegyVolumesList());
        setupTraceRecords();
        setValues();
        if (wizard.getVolumes().length > 0)
            wizardStep.getWizard().enableNext();
        else
            wizardStep.getWizard().disableNext();
    }

	private boolean ignoreChangeEvent = false;

	public boolean setValues()
	{
		if ((wizard.getSelectedVolumes() == null) || wizard.getSelectedVolumes().length < 1 || (fmtRecs == null))
		{
			return false;
		}

        StsSeismicBoundingBox volume = wizard.getFirstSegyVolume();
        wizardStep.initSetValues(volume);

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

		boolean littleEndian = volume.getIsLittleEndian();

		Object[] data = null;
		if (selectedTableChk.isSelected())
		{
			data = requiredRecords.toArray();
		}
		else
		{
			data = fmtRecs;
		}

		int xyScalePos = 0;
		int xyScaleFmt = 0;
		int edScalePos = 0;
		int edScaleFmt = 0;
		StsSEGYFormat segyFormat = wizard.getSegyFormat();
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
			traceBinaryHeader[col] = wizardStep.getTraceHeaderBinary(wizard.getSelectedVolumes()[0], trace);
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

	public boolean getOverrideHeader()
	{
		return this.overrideCoor;
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
		return requiredRecords.size();
	}

	public StsSEGYFormatRec getRequiredRecord(int i)
	{
		return (StsSEGYFormatRec)requiredRecords.get(i);
	}

	static public void main(String[] args)
	{
		JFrame frame = new JFrame("Test Panel");
		frame.setSize(300, 200);
//        frame.addKeyListener();

		Container contentPane = frame.getContentPane();
		StsSegyTraceDefinitionPanel panel = new StsSegyTraceDefinitionPanel(null, null);
		contentPane.add(panel);
		StsToolkit.centerComponentOnScreen(frame);
		frame.pack();
		frame.setVisible(true);
	}
}

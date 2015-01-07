package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsPreStackSurveyDefinitionPanel extends StsJPanel implements StsSelectRowNotifyListener
{
	private StsPreStackWizard wizard;
    private StsPreStackSurveyDefinition wizardStep;
    private StsSeismicBoundingBox[] volumes = null;
	private boolean override = false;
	private boolean useDecimation = false;
	private int inlineDecimation = 1;
	private int xlineDecimation = 1;
	private boolean isSurveyEditable = false;

    private StsGroupBox surveyDefinitionBox;
	private StsJPanel overridePanel;
	private StsJPanel parameterPanel;
    private StsJPanel bottomPanel;
    private StsGroupBox buttonBox;
    private Sts3PointSurveyDefinitionPanel surveyDefinition3PointPanel;

    private StsBooleanFieldBean overrideChk;

    private StsComboBoxFieldBean surveyDefinitionsComboBox;
	private StsIntFieldBean inlineDecimationBean;
	private StsIntFieldBean xlineDecimationBean;
	private StsDoubleFieldBean xOriginBean;
	private StsDoubleFieldBean yOriginBean;
	private StsFloatFieldBean zStartBean;
	private StsFloatFieldBean inlineOriginBean;
	private StsFloatFieldBean xlineOriginBean;
	private StsFloatFieldBean xIncBean;
	private StsFloatFieldBean yInclBean;
	private StsFloatFieldBean zIncBean;
	private StsFloatFieldBean angleBean;
    private StsBooleanFieldBean isCCWBean;

    //private SurveyDefinition[] surveyDefinitions;
	private ArrayList surveyDefinitions = new ArrayList();
	private SurveyDefinition currentSurveyDefinition = null;
	private SurveyDefinition defaultSurveyDefinition = null;
	private SurveyDefinition userSurveyDefinition = null;
    private SurveyDefinition user3PointSurveyDefinition = null;

    private StsTablePanelNew volumeStatusTablePanel = null;

    public StsPreStackSurveyDefinitionPanel(StsPreStackWizard wizard, StsPreStackSurveyDefinition wizardStep)
	{
		super(true); // true adds insets
		this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
		{
            buildTablePanel();
			constructBeans();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "xOrigin", "yOrigin", "rowNumMin", "rowNumMax", "colNumMin", "colNumMax", "statusString"};
        String[] columnTitles = {"Name", "Start X", "Start Y", "Start Inline", "End Inline", "Start XLine", "End XLine", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Volumes");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }

	private void constructBeans()
	{
		overrideChk = new StsBooleanFieldBean(this, "override", false, "Override");
		surveyDefinitionsComboBox = new StsComboBoxFieldBean(this, "surveyDefinition", "Coor Source");
		surveyDefinitionsComboBox.setEditable(false);
		inlineDecimationBean = new StsIntFieldBean(this, "inlineDecimation", 0, 1000, "Inline Decimation:", true);
		xlineDecimationBean = new StsIntFieldBean(this, "xlineDecimation", 0, 1000, "Crossline Decimation:", true);

		xOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "xOrigin", false, "X Origin:", true);
		yOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "yOrigin", false, "Y Origin:", true);
		zStartBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "zMin", false, "Z minimum:", true);
		inlineOriginBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "rowNumMin", false, "Origin Inline:");
		xlineOriginBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "colNumMin", false, "Origin Crossline:");
		xIncBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "xInc", false, "X Interval:", true);
		yInclBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "yInc", false, "Y Interval:", true);
		zIncBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "zInc", false, "Z Interval:", true);
		angleBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "angle", 0.0f, 360.0f, "Angle:");
		angleBean.setEditable(false);
        isCCWBean = new StsBooleanFieldBean(StsSeismicBoundingBox.class, "isXLineCCW", "Line->XLine CCW");
    }

    public void rowsSelected( int[] selectedIndices)
    {
        try
        {
            wizard.setSkipReanalyzeTraces(true);
            wizard.setSelectedVolumes( selectedIndices);
        }
        finally
        {
            wizard.setSkipReanalyzeTraces(false);
        }
    }

	public void initialize()
	{
        wizard.setSkipReanalyzeTraces(true);
        volumes = wizard.getSegyVolumes();
        volumeStatusTablePanel.replaceRows( wizard.getSegyVolumesList());
        initialiseAvailableSurveyDefinitions();
        buildPanel();
        wizard.setSkipReanalyzeTraces(false);
	}

	// Survey definitions can come from already loaded volumes, the project, or the first segyLine being processed
	// in that order of priority
	private void initialiseAvailableSurveyDefinitions()
	{
		surveyDefinitions.clear();
		// add volumes being processed to surveys
		addSurveyDefinitions(volumes, "Pre-stack segy volume: ");
		defaultSurveyDefinition = (SurveyDefinition)surveyDefinitions.get(0);
		currentSurveyDefinition = defaultSurveyDefinition;

		// user-definable survey
		StsSeismicBoundingBox userDefinitionVolume = createUserSurveyVolume(volumes[0]);
		userSurveyDefinition = addSurveyDefinition(userDefinitionVolume, "User-defined volume ");

        StsSeismicBoundingBox userDefinition3PointVolume = create3PointUserSurveyVolume(volumes[0]);
		user3PointSurveyDefinition = addSurveyDefinition(userDefinition3PointVolume, "User-defined 3 point volume ");

        if (wizard.getModel().getProject().getGridDefinition() != null)
		{
			// add project to surveys
			StsRotatedGridBoundingBox boundingBox = wizard.getModel().getProject().getRotatedBoundingBox();
			StsSeismicBoundingBox seismicBoundingBox = new StsSeismicBoundingBox(boundingBox, false);
			addSurveyDefinition(seismicBoundingBox, "Project volume");
		}
        // add loaded prestack volumes to surveys
		StsObject[] preStackSeismicVolumes = wizard.getModel().getObjectList(StsPreStackLineSet3d.class);
		addSurveyDefinitions(preStackSeismicVolumes, "Pre-stack volume: ");

		// add loaded poststack volumes to surveys
		StsObject[] seismicVolumes = wizard.getModel().getObjectList(StsSeismicVolume.class);
		addSurveyDefinitions(seismicVolumes, "Post-stack volume: ");

	    surveyDefinitionsComboBox.setListItems(surveyDefinitions.toArray());
		surveyDefinitionsComboBox.setEditable(override);
		checkSetSurveyIsEditable();
	}

	private void addSurveyDefinitions(StsObject[] surveyObjects, String description)
	{
		for (int n = 0; n < surveyObjects.length; n++)
		{
			addSurveyDefinition(surveyObjects[n], description);
		}
	}

	private SurveyDefinition addSurveyDefinition(StsObject surveyObject, String description)
	{
		SurveyDefinition surveyDefinition = new SurveyDefinition((StsSeismicBoundingBox)surveyObject, description);
		surveyDefinitions.add(surveyDefinition);
		return surveyDefinition;
	}

    private void checkSetSurveyIsEditable()
	{
		isSurveyEditable = currentSurveyDefinition == userSurveyDefinition;
	}

	private StsSeismicBoundingBox createUserSurveyVolume(StsSeismicBoundingBox referenceVolume)
	{
		StsSeismicBoundingBox userDefinitionVolume = new StsSeismicBoundingBox(referenceVolume, false);
        userDefinitionVolume.xMin = 0.0f;
        userDefinitionVolume.yMin = 0.0f;
		return userDefinitionVolume;
	}

	private StsSeismicBoundingBox create3PointUserSurveyVolume(StsSeismicBoundingBox referenceVolume)
	{
		StsSeismicBoundingBox userDefinitionVolume = new StsSeismicBoundingBox(referenceVolume, false);
        userDefinitionVolume.xMin = 0.0f;
        userDefinitionVolume.yMin = 0.0f;
		return userDefinitionVolume;
	}

    private void buildPanel()
	{
        removeAll();

        gbc.fill = GridBagConstraints.BOTH;
		volumeStatusTablePanel.addSelectRowNotifyListener(this);
        add(volumeStatusTablePanel);

        surveyDefinitionBox = new StsGroupBox("Survey Definition");
	    overridePanel = StsJPanel.addInsets();
	    parameterPanel = StsJPanel.addInsets();
        bottomPanel = StsJPanel.addInsets();
        buttonBox = new StsGroupBox();

        gbc.fill = GridBagConstraints.HORIZONTAL;

        surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		surveyDefinitionBox.gbc.anchor = gbc.NORTH;
		overridePanel.addToRow(overrideChk);
		overridePanel.addEndRow(surveyDefinitionsComboBox);
		surveyDefinitionBox.add(overridePanel);

		parameterPanel.addToRow(xOriginBean);
		parameterPanel.addEndRow(yOriginBean);
		parameterPanel.addToRow(inlineOriginBean);
		parameterPanel.addEndRow(xlineOriginBean);
		parameterPanel.addToRow(xIncBean);
		parameterPanel.addEndRow(yInclBean);
		parameterPanel.addToRow(zStartBean);
		parameterPanel.addEndRow(zIncBean);
		parameterPanel.addToRow(angleBean);
        parameterPanel.addEndRow(isCCWBean);
        parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox);

        surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		surveyDefinitionBox.add(parameterPanel);

        bottomPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		bottomPanel.gbc.anchor = gbc.NORTH;
		bottomPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		bottomPanel.add(surveyDefinitionBox, 1, 0.0);
        if(override && currentSurveyDefinition == this.user3PointSurveyDefinition)
        {
            surveyDefinition3PointPanel = new Sts3PointSurveyDefinitionPanel(wizard, wizardStep);
            surveyDefinition3PointPanel.initialize(volumes[0]);
            bottomPanel.add(surveyDefinition3PointPanel);
            StsButton acceptBtn = new StsButton("Accept", "Recompute geometry based on this trace info.", this, "accept");
            StsButton cancelBtn = new StsButton("Cancel", "Ignore this trace info.", this, "cancel");
            buttonBox.addToRow(acceptBtn);
            buttonBox.addEndRow(cancelBtn);
            bottomPanel.add(buttonBox);
        }
        bottomPanel.gbc.weighty = 1.0;
		bottomPanel.gbc.fill = GridBagConstraints.BOTH;
		bottomPanel.addEndRow(Box.createVerticalGlue(), 1, 1.0);

        gbc.fill = GridBagConstraints.BOTH;
		add(bottomPanel);

/*
		jSplitPane1.setContinuousLayout(false);
		jSplitPane1.setOrientation(JSplitPane.VERTICAL_SPLIT);
		jSplitPane1.setOneTouchExpandable(true);

		jSplitPane1.add(volumeStatusTablePanel, JSplitPane.TOP);
		jSplitPane1.add(bottomPanel, JSplitPane.BOTTOM);

		this.gbc.fill = this.gbc.BOTH;

		topPanel.gbc.anchor = topPanel.gbc.NORTH;
		topPanel.gbc.fill = topPanel.gbc.BOTH;
		topPanel.add(jSplitPane1);
		this.gbc.fill = this.gbc.BOTH;
		this.add(topPanel);
*/

        wizard.rebuild();
    }

	public boolean getOverride() { return override; }

	public void setOverride(boolean value)
	{
		override = value;
		if(override)
		{
			parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox);
			surveyDefinitionsComboBox.setEditable(true);
			parameterPanel.setEditable(isSurveyEditable);
		}
		else
		{
			currentSurveyDefinition = defaultSurveyDefinition;
            surveyDefinitionsComboBox.setValueObject(currentSurveyDefinition);
            parameterPanel.setPanelObject(defaultSurveyDefinition.boundingBox);
			surveyDefinitionsComboBox.setEditable(false);
			parameterPanel.setEditable(false);
		}
        buildPanel();
    }

    public int getInlineDecimation() { return inlineDecimation; }
	public void setInlineDecimation(int value) { inlineDecimation = value; }
	public int getXlineDecimation() { return xlineDecimation; }
	public void setXlineDecimation(int value) { xlineDecimation = value; }
	public boolean getUseDecimation() { return useDecimation; }

	public SurveyDefinition getSurveyDefinition()
	{
		return this.currentSurveyDefinition;
	}

	private void setSurveyDefinition(SurveyDefinition surveyDefinition)
	{
        if(currentSurveyDefinition == surveyDefinition) return;

        currentSurveyDefinition = surveyDefinition;
		parameterPanel.setPanelObject(surveyDefinition.boundingBox);
		checkSetSurveyIsEditable();
		setBeansEditable(isSurveyEditable);
		if(currentSurveyDefinition == defaultSurveyDefinition)
		{
			inlineDecimation = Math.round(volumes[0].rowNumInc/surveyDefinition.boundingBox.rowNumInc);
			xlineDecimation = Math.round(volumes[0].colNumInc/surveyDefinition.boundingBox.colNumInc);
			useDecimation = false;
		}
        else if(currentSurveyDefinition == user3PointSurveyDefinition)
        {
			inlineDecimation = 1;
			xlineDecimation = 1;
			useDecimation = true;
        }
        else
		{
			inlineDecimation = 1;
			xlineDecimation = 1;
			useDecimation = true;
		}
		inlineDecimationBean.setValue(inlineDecimation);
		xlineDecimationBean.setValue(xlineDecimation);

        buildPanel();
    }

	private void setBeansEditable(boolean editable)
	{
		parameterPanel.setEditable(editable);
		angleBean.setEditable(editable);
	}

	public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
	{
		return this.currentSurveyDefinition.boundingBox;
	}

	// This method is ultimately called by the trace analysis callback
	public void updatePanel()
	{
        volumeStatusTablePanel.replaceRows( wizard.getSegyVolumesList());
        if (wizard.getSegyVolumesList().size() >  1)
		{
			wizardStep.enableNext();
		}
		else
		{
			wizard.disableNext();
		}
	}

	class SurveyDefinition
	{
		String surveyNameString;
		StsSeismicBoundingBox boundingBox;

        SurveyDefinition(StsSeismicBoundingBox boundingBox, String description)
		{
			String name = boundingBox.getName();
			if(name == null)
				surveyNameString = description;
			else
				surveyNameString = new String(description + name);
			this.boundingBox = boundingBox;
		}

	    public String toString() { return surveyNameString; }
	}

    class SurveyDefinition3Point extends SurveyDefinition
    {
        SurveyDefinition3Point(StsSeismicBoundingBox boundingBox, String description)
        {
            super(boundingBox, description);   
        }
    }

    public void accept()
    {
        if(volumes[0].checkKnownPointsOK()) return;

        new StsMessage(null, StsMessage.INFO,
            "All fields must be specified. Fill in missing fields or cancel.");
        return;
    }


    public void cancel()
    {
        return;
    }
}

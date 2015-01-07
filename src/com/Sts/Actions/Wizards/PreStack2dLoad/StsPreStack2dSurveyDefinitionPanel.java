package com.Sts.Actions.Wizards.PreStack2dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsPreStack2dSurveyDefinitionPanel extends StsJPanel
{
	private StsPreStackLoad2dWizard wizard;
	private StsPreStackLineSet3d[] volumes = null;
	private boolean override = false;
//	private boolean migrated = false;
	private boolean useDecimation = false;
//	private int inlineDecimation = 1;
//	private int xlineDecimation = 1;
	private boolean isSurveyEditable = false;


	StsGroupBox surveyDefinitionBox = new StsGroupBox("Survey Definition");
	StsJPanel overridePanel = StsJPanel.addInsets();
	StsJPanel parameterPanel = StsJPanel.addInsets();
	StsJPanel infoBox = StsJPanel.addInsets();

	StsBooleanFieldBean overrideChk = new StsBooleanFieldBean();
	//StsBooleanFieldBean migratedChk = new StsBooleanFieldBean();

	StsComboBoxFieldBean surveyDefinitionsComboBox;
	StsDoubleFieldBean xOriginBean;
	StsDoubleFieldBean yOriginBean;
	StsFloatFieldBean zStartBean;
	StsFloatFieldBean zIncBean;

	StsFloatFieldBean xMinBean;
	StsFloatFieldBean yMinBean;
	StsFloatFieldBean xMaxBean;
	StsFloatFieldBean yMaxBean;

	private transient SurveyDefinition[] surveyDefinitions;
	private transient SurveyDefinition currentSurveyDefinition = null;
	private transient SurveyDefinition defaultSurveyDefinition = null;
	private transient SurveyDefinition userSurveyDefinition = null;

	public StsPreStack2dSurveyDefinitionPanel(StsPreStackLoad2dWizard wizard, StsWizardStep wizardStep)
	{
		super(true); // true adds insets
		this.wizard = wizard;
		setSize(700,500);
		try
		{
			constructBeans();
			jbInit();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public void constructBeans()
	{
		overrideChk = new StsBooleanFieldBean(this, "override", false, "Override");
		surveyDefinitionsComboBox = new StsComboBoxFieldBean(this, "surveyDefinition", "Coor Source");
		surveyDefinitionsComboBox.setEditable(override);
		//migratedChk = new StsBooleanFieldBean(this, "migrated", false, "Migrated");
		xOriginBean = new StsDoubleFieldBean(StsRotatedGridBoundingBox2.class, "xOrigin", false, "X Origin:", true);
		yOriginBean = new StsDoubleFieldBean(StsRotatedGridBoundingBox2.class, "yOrigin", false, "Y Origin:", true);

		xMinBean = new StsFloatFieldBean(StsRotatedGridBoundingBox2.class, "xMin", false, "X Min:", true);
		yMinBean = new StsFloatFieldBean(StsRotatedGridBoundingBox2.class, "yMin", false, "Y Min:", true);
		xMaxBean = new StsFloatFieldBean(StsRotatedGridBoundingBox2.class, "xMax", false, "X Max:", true);
		yMaxBean = new StsFloatFieldBean(StsRotatedGridBoundingBox2.class, "yMax", false, "Y Max:", true);

		zStartBean = new StsFloatFieldBean(StsRotatedGridBoundingBox.class, "zMin", false, "Z minimum:", true);
		zIncBean = new StsFloatFieldBean(StsRotatedGridBoundingBox2.class, "zInc", false, "Z Interval:", true);
	}

	public void initialize()
	{
		initialiseAvailableSurveyDefinitions();
		parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox2);
	}

	// Survey definitions can come from already loaded volumes, the project
	private void initialiseAvailableSurveyDefinitions()
	{
        surveyDefinitionsComboBox.removeAll();
        surveyDefinitions = null;

		// add lines being processed to surveys
		StsRotatedGridBoundingBox lineMerge = wizard.getLinesBoundingBoxes();
		SurveyDefinition prestackSegyVolumeSurveys = addSurveyDefinition(lineMerge, "Prestack segy Line(s): ");
		defaultSurveyDefinition = prestackSegyVolumeSurveys;
		currentSurveyDefinition = defaultSurveyDefinition;
		setSurveyDefinition(currentSurveyDefinition);

		// user-definable survey
		StsRotatedGridBoundingBox userDefinitionVolume = createUserSurveyVolume(lineMerge);
		userSurveyDefinition = addSurveyDefinition(userDefinitionVolume, "User-defined volume ");

		// add project to surveys
		StsRotatedGridBoundingBox boundingBox = wizard.getModel().getProject().getRotatedBoundingBox();
		addSurveyDefinition(boundingBox, "Project volume");

		// add loaded prestack volumes to surveys
		StsObject[] preStackSeismicVolumes = wizard.getModel().getObjectList(StsPreStackLineSet3d.class);
		addSurveyDefinitions(preStackSeismicVolumes, "Prestack volume: ");

		// add loaded poststack volumes to surveys
		StsObject[] seismicVolumes = wizard.getModel().getObjectList(StsSeismicVolume.class);
		addSurveyDefinitions(seismicVolumes, "Poststack volume: ");

		surveyDefinitionsComboBox.setListItems(surveyDefinitions);
		surveyDefinitionsComboBox.setEditable(override);
		checkSetSurveyIsEditable();
        configureUI();
	}

	private void checkSetSurveyIsEditable()
	{
		isSurveyEditable = (currentSurveyDefinition == defaultSurveyDefinition || currentSurveyDefinition == userSurveyDefinition);
	}

	private StsRotatedGridBoundingBox createUserSurveyVolume(StsRotatedGridBoundingBox referenceVolume)
	{
		StsRotatedGridBoundingBox userDefinitionVolume = new StsRotatedGridBoundingBox(referenceVolume);
		if(userDefinitionVolume.xMin == StsParameters.largeFloat)  userDefinitionVolume.xMin = 0.0f;
		if(userDefinitionVolume.yMin == StsParameters.largeFloat)  userDefinitionVolume.yMin = 0.0f;
		return userDefinitionVolume;
	}

	private void jbInit() throws Exception
	{
		gbc.fill = GridBagConstraints.HORIZONTAL;

		//add(migratedChk);

		surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		overridePanel.addToRow(overrideChk);
		overridePanel.addEndRow(surveyDefinitionsComboBox);
		surveyDefinitionBox.add(overridePanel);

//		decimationPanel.addToRow(inlineDecimationBean);
//		decimationPanel.addEndRow(xlineDecimationBean);
//		surveyDefinitionBox.add(decimationPanel);

		parameterPanel.addToRow(xOriginBean);
		parameterPanel.addEndRow(yOriginBean);
		parameterPanel.addToRow(xMinBean);
		parameterPanel.addEndRow(yMinBean);
		parameterPanel.addToRow(xMaxBean);
		parameterPanel.addEndRow(yMaxBean);

//		parameterPanel.addToRow(inlineOriginBean);
//		parameterPanel.addEndRow(xlineOriginBean);

//		parameterPanel.addToRow(xIncBean);
//		parameterPanel.addEndRow(yIncBean);

		parameterPanel.addToRow(zStartBean);
		parameterPanel.addEndRow(zIncBean);

//		parameterPanel.addToRow(nXBean);
//		parameterPanel.addEndRow(nYBean);

//		parameterPanel.add(angleBean);

		surveyDefinitionBox.add(parameterPanel);

		add(surveyDefinitionBox);

		gbc.fill = GridBagConstraints.BOTH;
		add(infoBox);
		//setOverride(true);
	}

	public boolean getOverride() { return override; }
	public void setOverride(boolean value)
	{
		override = value;
		configureUI();
	}

    public void configureUI()
    {
        if(override)
        {
            parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox2);
            surveyDefinitionsComboBox.setEditable(true);
            parameterPanel.setEditable(isSurveyEditable);
        }
        else
        {
            parameterPanel.setPanelObject(defaultSurveyDefinition.boundingBox2);
            surveyDefinitionsComboBox.setSelectedIndex(0);
            surveyDefinitionsComboBox.setEditable(false);
            parameterPanel.setEditable(false);
        }
    }

//	public boolean getMigrated() { return migrated; }
//	public void setMigrated(boolean value) { migrated = value; }

//	public int getInlineDecimation() { return inlineDecimation; }
//	public void setInlineDecimation(int value) { inlineDecimation = value; }
//	public int getXlineDecimation() { return xlineDecimation; }
//	public void setXlineDecimation(int value) { xlineDecimation = value; }
	public boolean getUseDecimation() { return useDecimation; }

	public SurveyDefinition getSurveyDefinition() { return this.currentSurveyDefinition; }


	private void setSurveyDefinition(SurveyDefinition surveyDefinition)
	{
		currentSurveyDefinition = surveyDefinition;
		parameterPanel.setPanelObject(surveyDefinition.boundingBox2);
		checkSetSurveyIsEditable();
		setBeansEditable(isSurveyEditable);
     	{
//			inlineDecimation = 1;
//			xlineDecimation = 1;
//			useDecimation = true;
		}
//		inlineDecimationBean.setValue(inlineDecimation);
//		xlineDecimationBean.setValue(xlineDecimation);

//		surveyDefinition.boundingBox2.computeNxNy();
	}

	private void setBeansEditable(boolean editable)
	{
//		surveyDefinitionsComboBox.setEditable(editable);
//		decimationPanel.setEditable(editable);
		parameterPanel.setEditable(editable);
	}

	public StsRotatedGridBoundingBox2 getSurveyDefinitionBoundingBox() { return this.currentSurveyDefinition == null ? null: this.currentSurveyDefinition.boundingBox2; }

	private SurveyDefinition[] addSurveyDefinitions(StsObject[] surveyObjects, String description)
	{
		SurveyDefinition[] surveyDefinitions = new SurveyDefinition[surveyObjects.length];
		for(int n = 0; n < surveyObjects.length; n++)
			surveyDefinitions[n] = addSurveyDefinition(surveyObjects[n], description);
		return surveyDefinitions;
	}

	private SurveyDefinition addSurveyDefinition(StsObject surveyObject, String description)
	{
		SurveyDefinition surveyDefinition = new SurveyDefinition((StsRotatedGridBoundingBox)surveyObject, description);
		surveyDefinitions = (SurveyDefinition[])StsMath.arrayAddElement(surveyDefinitions, surveyDefinition, SurveyDefinition.class);
		return surveyDefinition;
	}


	class SurveyDefinition
	{
		String surveyNameString;
		transient StsRotatedGridBoundingBox2 boundingBox2;

		SurveyDefinition(StsRotatedGridBoundingBox boundingBox, String description)
		{
			boundingBox2 = new StsRotatedGridBoundingBox2(false);
			boundingBox2.initializeToBoundingBox(boundingBox);
			String name = boundingBox2.getName();
			if(name == null)
				surveyNameString = description;
			else
				surveyNameString = new String(description + name);
			this.boundingBox2 = boundingBox2;
		}

		public String toString() { return surveyNameString; }
	}

	class StsRotatedGridBoundingBox2 extends StsRotatedGridBoundingBox
	{
//		public long nX, nY;
		/** default constructor */
		public StsRotatedGridBoundingBox2()
		{
			super(false);
		}
		public StsRotatedGridBoundingBox2(boolean persistent)
        {
	        super(persistent);
        }
		public StsRotatedGridBoundingBox2(StsRotatedGridBoundingBox box)
		{
			super(box);
		}
    /*
        public void computeNxNy()
		{
			nX = (long)((xMax - xMin)/xInc);
			nY = (long)((yMax - yMin)/yInc);

			nXBean.setValue(nX);
			nYBean.setValue(nY);
		}

		public void computeXincYinc()
		{
			xInc = ((xMax - xMin)/((float) nX));
			yInc = ((yMax - yMin)/((float) nY));
			xIncBean.setValue(xInc);
			yIncBean.setValue(yInc);
		}

		public void setNX(long nx) { this.nX = nx; computeXincYinc();}
		public long getNX() { return nX; }
		public void setNY(long ny) { this.nY = ny; computeXincYinc();}
		public long getNY() { return nY; }

		public void setXInc(float xInc)
		{
			this.xInc = xInc;
			computeNxNy();
		}
		public void setYInc(float yInc)
		{
			this.yInc = yInc;
			computeNxNy();
		}
		public float getYInc()
		{
			return yInc;
		}
		public float getXInc()
		{
			return xInc;
		}
    */
		public void setXMin(float xMin)
		{
			super.setXMin(xMin);
//			computeNxNy();
		}
		public void setYMin(float yMin)
		{
			super.setYMin(yMin);
//			computeNxNy();
		}
		public void setXMax(float xMax)
		{
			super.setXMax(xMax);
//			computeNxNy();
		}
		public void setYMax(float yMax)
		{
			super.setYMax(yMax);
//			computeNxNy();
		}

	}
}

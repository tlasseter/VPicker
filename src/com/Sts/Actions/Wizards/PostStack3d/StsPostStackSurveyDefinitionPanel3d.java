package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

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

public class StsPostStackSurveyDefinitionPanel3d extends StsJPanel implements StsSelectRowNotifyListener
{
	private StsSeismicWizard wizard;
    private StsPostStackSurveyDefinition3d wizardStep;
    private StsGroupBox surveyDefinitionBox;
	private StsJPanel overridePanel;
	private StsJPanel parameterPanel;
    private StsJPanel surveyDefinition3PointPanel;

	private StsComboBoxFieldBean surveyDefinitionsComboBox;
	private StsDoubleFieldBean xOriginBean;
	private StsDoubleFieldBean yOriginBean;
	private StsFloatFieldBean inlineOriginBean;
	private StsFloatFieldBean xlineOriginBean;
	private StsFloatFieldBean xIncBean;
	private StsFloatFieldBean yIncBean;
	private StsFloatFieldBean angleBean;
    private StsBooleanFieldBean isCCWBean;
    private StsFloatFieldBean zMinBean;
    private StsFloatFieldBean zIncBean;

	private ArrayList surveyDefinitions = new ArrayList();
	private StsPostStackSurveyDefinitionPanel3d.SurveyDefinition currentSurveyDefinition = null;
	private StsPostStackSurveyDefinitionPanel3d.SurveyDefinition userSurveyDefinition = null;
    private StsPostStackSurveyDefinitionPanel3d.SurveyDefinition3Point user3PointSurveyDefinition = null;

    private StsTablePanelNew statusTablePanel = null;
    private StsButton analyzeButton;
    public StsProgressPanel progressPanel;

    static final public float nullValue = StsParameters.nullValue;
    static final public double doubleNullValue = StsParameters.doubleNullValue;

    public StsPostStackSurveyDefinitionPanel3d(StsSeismicWizard wizard, StsPostStackSurveyDefinition3d wizardStep)
	{
		super(true); // true adds insets
		this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
		{
            buildTablePanel();
            constructBeans();
            buildPanel();
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
        statusTablePanel = new StsTablePanelNew(tableModel);
        statusTablePanel.setLabel("Volumes");
        statusTablePanel.setSize(400, 100);
        statusTablePanel.initialize();
    }

	private void constructBeans()
	{
        currentSurveyDefinition = new SurveyDefinition("none");
        surveyDefinitionsComboBox = new StsComboBoxFieldBean(this, "surveyDefinition", "Coor Source", new Object[] { currentSurveyDefinition });
		xOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "xOrigin", "X Origin:");
		yOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "yOrigin", "Y Origin:");
		inlineOriginBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "rowNumMin", "Origin Inline:");
		xlineOriginBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "colNumMin", "Origin Crossline:");
		xIncBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "xInc", "X Interval:");
		yIncBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "yInc", "Y Interval:");
		angleBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "angle", 0.0f, 360.0f, "Angle:");
        isCCWBean = new StsBooleanFieldBean(StsSeismicBoundingBox.class, "isXLineCCW", "Line->XLine CCW");
        zMinBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "zMin", "Z min:");
        zIncBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "zInc", "Z Inc:");
        analyzeButton = new StsButton("Analyze", "Analyze geometry.", this, "analyzeAngle");
        progressPanel = StsProgressPanel.constructor(5, 50);
    }

	public void initialize()
	{
		wizard.setSkipReanalyzeTraces(true);
        statusTablePanel.replaceRows( wizard.getSegyVolumesList());
        statusTablePanel.setSelectAll(true);
        initialiseAvailableSurveyDefinitions();
	    buildPanel();
		wizard.setSkipReanalyzeTraces(false);
        wizard.disableNext();
    }

    public void analyzeAngle()
    {
	    currentSurveyDefinition.analyzeGeometry();
    }


    // Survey definitions can come from already loaded volumes, the project, or the first segyLine being processed
	// in that order of priority
	private void initialiseAvailableSurveyDefinitions()
	{
		surveyDefinitions.clear();
        if(wizard.getSegyVolumesList().size() == 0) return;
        StsSegyVolume segyVolume = (StsSegyVolume)wizard.getSegyVolumesList().get(0);

        // User-definable survey
		StsSeismicBoundingBox userDefinitionVolume = createUserSurveyVolume( segyVolume);
		userSurveyDefinition = addSurveyDefinition(userDefinitionVolume, "User-defined volume ");
        currentSurveyDefinition = userSurveyDefinition;
        StsSeismicBoundingBox userDefinition3PointVolume = create3PointUserSurveyVolume( segyVolume);
		user3PointSurveyDefinition = add3PointSurveyDefinition(userDefinition3PointVolume, "User-defined 3 point volume ");

        if (wizard.getModel().getProject().getGridDefinition() != null)
		{
			// Add project to surveys
			StsRotatedGridBoundingBox boundingBox = wizard.getModel().getProject().getRotatedBoundingBox();
			StsRotatedGridBoundingBox seismicBoundingBox = new StsRotatedGridBoundingBox(boundingBox, false);
			addSurveyDefinition(seismicBoundingBox, "Project volume");
		}

		// Add loaded poststack volumes to surveys
		StsObject[] seismicVolumes = wizard.getModel().getObjectList(StsSeismicVolume.class);
		addSurveyDefinitions(seismicVolumes, "Post-stack volume: ");

        StsObject[] prestackVolumes = wizard.getModel().getObjectList(StsPreStackLineSet3d.class);
		addSurveyDefinitions(prestackVolumes, "Pre-stack volume: ");

        surveyDefinitionsComboBox.setListItems(surveyDefinitions.toArray());
	}

	private void addSurveyDefinitions(StsObject[] surveyObjects, String description)
	{
		for (int n = 0; n < surveyObjects.length; n++)
			addSurveyDefinition(surveyObjects[n], description);
	}

	private StsPostStackSurveyDefinitionPanel3d.SurveyDefinition addSurveyDefinition(StsObject surveyObject, String description)
	{
		SurveyDefinition surveyDefinition = new SurveyDefinition((StsSeismicBoundingBox)surveyObject, description);
		surveyDefinitions.add(surveyDefinition);
		return surveyDefinition;
	}


	private SurveyDefinition3Point add3PointSurveyDefinition(StsObject surveyObject, String description)
	{
		SurveyDefinition3Point surveyDefinition = new SurveyDefinition3Point((StsSeismicBoundingBox)surveyObject, description);
		surveyDefinitions.add(surveyDefinition);
		return surveyDefinition;
	}

	private StsSeismicBoundingBox createUserSurveyVolume(StsSeismicBoundingBox referenceVolume)
	{
		StsSeismicBoundingBox userDefinitionVolume = new StsSeismicBoundingBox(referenceVolume, false);
        initializeGeometry(userDefinitionVolume);
        return userDefinitionVolume;
	}

	private StsSeismicBoundingBox create3PointUserSurveyVolume(StsSeismicBoundingBox referenceVolume)
	{
		StsSeismicBoundingBox userDefinitionVolume = new StsSeismicBoundingBox(referenceVolume, false);
        initializeGeometry(userDefinitionVolume);
        userDefinitionVolume.setHorizontalScalar(referenceVolume.getHorizontalScalar());
		return userDefinitionVolume;
	}

    private void initializeGeometry(StsSeismicBoundingBox userDefinitionVolume)
    {
        userDefinitionVolume.xMin = 0.0f;
        userDefinitionVolume.yMin = 0.0f;
        userDefinitionVolume.xInc = 1.0f;
        userDefinitionVolume.yInc = 1.0f;
        userDefinitionVolume.xMax = userDefinitionVolume.nCols-1;
        userDefinitionVolume.yMin = userDefinitionVolume.nRows-1;
    }

    private void buildPanel()
	{
        removeAll();

        gbc.fill = GridBagConstraints.BOTH;
		statusTablePanel.addSelectRowNotifyListener(this);
        add(statusTablePanel);

        surveyDefinitionBox = new StsGroupBox("Survey Definition");
	    overridePanel = StsJPanel.addInsets();
	    parameterPanel = StsJPanel.addInsets();

	    surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		surveyDefinitionBox.gbc.anchor = gbc.NORTH;
		overridePanel.add(surveyDefinitionsComboBox);
        overridePanel.gbc.anchor = gbc.CENTER;
        overridePanel.gbc.gridwidth = 2;
        overridePanel.add(analyzeButton);
        surveyDefinitionBox.add(overridePanel);

		parameterPanel.addToRow(xOriginBean);
		parameterPanel.addEndRow(yOriginBean);
		parameterPanel.addToRow(inlineOriginBean);
		parameterPanel.addEndRow(xlineOriginBean);
		parameterPanel.addToRow(xIncBean);
		parameterPanel.addEndRow(yIncBean);
		parameterPanel.addToRow(zMinBean);
        parameterPanel.addEndRow(zIncBean);
		parameterPanel.addToRow(angleBean);
        parameterPanel.addEndRow(isCCWBean);
        parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox);

        surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		surveyDefinitionBox.add(parameterPanel);
        add(surveyDefinitionBox);
        if(currentSurveyDefinition == this.user3PointSurveyDefinition)
        {
            surveyDefinition3PointPanel = user3PointSurveyDefinition.get3PointSurveyDefinitionPanel();
            add(surveyDefinition3PointPanel);
        }
        gbc.anchor = GridBagConstraints.SOUTH;
        add(progressPanel);

        wizard.rebuild();
    }

    public float getZInc()
	{
		return currentSurveyDefinition.boundingBox.getZInc();
	}

    public float getZMin()
	{
		return currentSurveyDefinition.boundingBox.getZMin();
	}

    public void setZInc(float zInc)
	{
		currentSurveyDefinition.boundingBox.setZInc(zInc);
	}

    public void setZMin(float zMin)
	{
		currentSurveyDefinition.boundingBox.setZMin(zMin);
	}

    public StsPostStackSurveyDefinitionPanel3d.SurveyDefinition getSurveyDefinition()
	{
		return this.currentSurveyDefinition;
	}

    private void setSurveyDefinition(StsPostStackSurveyDefinitionPanel3d.SurveyDefinition surveyDefinition)
    {
        if(currentSurveyDefinition == surveyDefinition) return;
        currentSurveyDefinition = surveyDefinition;
        buildPanel();
        parameterPanel.setPanelObject(surveyDefinition.boundingBox);
        boolean editable = surveyDefinition instanceof SurveyDefinition3Point || surveyDefinition instanceof SurveyDefinitionUserDefined;
        parameterPanel.setEnabled(editable);
    }

     public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
     {
         return this.currentSurveyDefinition.boundingBox;
     }

	public void rowsSelected( int[] indices) {}

    class SurveyDefinition
	{
		String surveyNameString;
		StsSeismicBoundingBox boundingBox;

        SurveyDefinition(String name)
		{
            boundingBox = new StsSeismicBoundingBox(false);
            surveyNameString = name;
		}

        SurveyDefinition(StsSeismicBoundingBox boundingBox, String description)
		{
			String name = boundingBox.getName();
			if(name == null)
				surveyNameString = description;
			else
				surveyNameString = new String(description + name);
			this.boundingBox = boundingBox;
		}

        /** this is implementation for a user-defined survey definition volume.
         *  We are directly setting the angle, so we return true indicating angle analysis has been completed.
         * @return true if angle analysis has been completed (running angleAnalyzer is unnecessary).
         */
        boolean analyzeGeometry()
        {
            progressPanel.initialize(1);
            boolean ok = true;
            if(boundingBox.xInc <= 0)
            {
                ok = false;
                progressPanel.appendLine("xInc cannot be <= zero.");
            }
            if(boundingBox.yInc <= 0)
            {
                ok = false;
                progressPanel.appendLine("yInc cannot be <= zero.");
            }
            progressPanel.finished();
            if(ok)
            {
                double xOrigin = boundingBox.getXOrigin();
                double yOrigin = boundingBox.getYOrigin();
                float ilineOrigin = boundingBox.getRowNumMin();
                float xlineOrigin = boundingBox.getColNumMin();
                StsSeismicBoundingBox[] volumes = wizard.getSelectedVolumes();
                for(int n = 0; n < volumes.length; n++)
                {
                    StsSeismicBoundingBox volume = volumes[n];
                    if( volume.status < StsSeismicBoundingBox.STATUS_TRACES_OK && volume.status > StsSeismicBoundingBox.STATUS_GRID_BAD) continue;
                    volume.xInc = boundingBox.xInc;
                    volume.yInc = boundingBox.yInc;
                    volume.setAngle(boundingBox.angle);
                    volume.xMin = boundingBox.getXFromColNum(volume.colNumMin);
                    volume.yMin = boundingBox.getYFromRowNum(volume.rowNumMin);
                    volume.xMax = boundingBox.getXFromColNum(volume.colNumMax);
                    volume.yMax = boundingBox.getYFromRowNum(volume.rowNumMax);
                    volume.xOrigin = xOrigin;
                    volume.yOrigin = yOrigin;
                    setFileStatus(StsSeismicBoundingBox.STATUS_GEOMETRY_OK, volume);
//                    volume.setOverrideGeometry(boundingBox);
                }
                progressPanel.setDescriptionAndLevel("Survey override OK", StsProgressBar.INFO);
                wizard.enableNext();
            }
            else
            {
                progressPanel.setDescriptionAndLevel("Override parameters bad", StsProgressBar.ERROR);
                return false;
            }
            return false;
        }


        private void setFileStatus( int status, StsSeismicBoundingBox volume)
        {
            if( volume != null)
            {
                volume.status = status;
                statusTablePanel.setValueAt( StsSeismicBoundingBox.statusText[volume.status], volume, "statusString");
            }
        }

        public String toString() { return surveyNameString; }
	}

    class SurveyDefinitionUserDefined extends SurveyDefinition
    {
       SurveyDefinitionUserDefined(StsSeismicBoundingBox boundingBox, String description)
        {
            super(boundingBox, description);
        }
    }

    class SurveyDefinition3Point extends SurveyDefinition
    {
        KnownPoint[] knownPoints = new KnownPoint[3];

        SurveyDefinition3Point(StsSeismicBoundingBox boundingBox, String description)
        {
            super(boundingBox, description);
            for(int n = 0; n < 3; n++)
                knownPoints[n] = new KnownPoint();
        }


        StsJPanel get3PointSurveyDefinitionPanel()
        {
            StsGroupBox knownBox = new StsGroupBox("Known points");
            knownBox.gbc.fill = gbc.HORIZONTAL;
            for(int n = 0; n < 3; n++)
            {
                StsFieldBean[] beans = knownPoints[n].getFieldBeans();
                knownBox.addRowOfBeans(beans);
            }
            return knownBox;
        }

        /** this is implementation for a a 3-point user-defined survey definition volume.
         *  We need to run the angleAnalyzer for this case, so we return false.
         * @return false as angle analysis needs to be run.
         */
        boolean analyzeGeometry()
        {
            if(!knownPointsOK())
            {
                new StsMessage(null, StsMessage.INFO, "All fields must be specified. Fill in missing fields or cancel.");
                return false;
            }

            StsSegyVolume segyVolume = (StsSegyVolume)wizard.getSegyVolumesList().get(0);
            analyzeGeometryWithKnownPoints( segyVolume.getSegyFormat());
            super.analyzeGeometry();
            return false;
        }

        private boolean knownPointsOK()
        {
            for(int i = 0; i < 3; i++)
                if(!knownPoints[i].ok()) return false;
            return true;
        }

        boolean analyzeGeometryWithKnownPoints(StsSEGYFormat segyFormat)
        {
            StsSEGYFormat.TraceHeader[] surveyTraces = new StsSEGYFormat.TraceHeader[3];
            float horizontalScalar = boundingBox.getHorizontalScalar();
            for(int n = 0; n < 3; n++)
            {
                surveyTraces[n] = segyFormat.constructTraceHeader();
                surveyTraces[n].x = (float)knownPoints[n].x * horizontalScalar;
                surveyTraces[n].y = (float)knownPoints[n].y * horizontalScalar;
                surveyTraces[n].xLine = knownPoints[n].crossline;
                surveyTraces[n].iLine = knownPoints[n].inline;
            }

            // compute rotation angle and bin spacings in line and xline directions
            if(!boundingBox.analyzeAngle(surveyTraces[0], surveyTraces[1], surveyTraces[2]))
                return false;
            boundingBox.initializeRange(surveyTraces[0]);

            boundingBox.xOrigin = knownPoints[0].x * horizontalScalar;
            boundingBox.yOrigin = knownPoints[0].y * horizontalScalar;

            return true;
        }

        class KnownPoint
        {
            float inline = nullValue;
            float crossline = nullValue;
            double x = doubleNullValue;
            double y = doubleNullValue;
            StsFieldBean[] beans = new StsFieldBean[]
            {
                new StsDoubleFieldBean(this, "x", "Point 1: X"),
                new StsDoubleFieldBean(this, "y", "Y"),
                new StsFloatFieldBean(this, "inline", "Inline"),
                new StsFloatFieldBean(this, "xline", "Xline")
            };

            KnownPoint()
            {
            }

            StsFieldBean[] getFieldBeans() { return beans; }

            boolean ok()
            {
                if(x == doubleNullValue) return false;
                if(y == doubleNullValue) return false;
                if(inline == nullValue) return false;
                if(crossline == nullValue) return false;
                return true;
            }

            double getX() { return x; }
            void setX(double x) { this.x = x; }

            double getY() { return y; }
            void setY(double y) { this.y = y; }

            float getInline() { return inline; }
            void setInline(float inline) { this.inline = inline; }

            float getXline() { return crossline; }
            void setXline(float crossline) { this.crossline = crossline; }
        }
    }
}
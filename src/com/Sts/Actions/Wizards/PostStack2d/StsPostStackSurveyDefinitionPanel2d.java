package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

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

public class StsPostStackSurveyDefinitionPanel2d extends StsJPanel implements StsSelectRowNotifyListener, StsTableModelListener
{
	private StsPostStack2dWizard wizard;
    private StsPostStackSurveyDefinition2d wizardStep;
    private StsSeismicBoundingBox[] volumes = null;
	private boolean isSurveyEditable = false;

    private StsGroupBox surveyDefinitionBox;
	private StsJPanel surveyPanel;
	private StsJPanel parameterPanel;

	private StsComboBoxFieldBean surveyDefinitionsComboBox;
	private StsDoubleFieldBean xOriginBean;
	private StsDoubleFieldBean yOriginBean;
    private StsFloatFieldBean xMinBean;
    private StsFloatFieldBean yMinBean;
    private StsFloatFieldBean xMaxBean;
    private StsFloatFieldBean yMaxBean;
    private StsFloatFieldBean cdpIntervalBean;

	private ArrayList surveyDefinitions = new ArrayList();
	private SurveyDefinition currentSurveyDefinition = new SurveyDefinition();
	private SurveyDefinition defaultSurveyDefinition = new SurveyDefinition();
	private SurveyDefinition userSurveyDefinition = new SurveyDefinition();

    private StsTablePanelNew volumeStatusTablePanel = null;

	public StsPostStackSurveyDefinitionPanel2d(StsPostStack2dWizard wizard, StsPostStackSurveyDefinition2d wizardStep)
	{
		super(true); 
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
        String[] columnNames = {"stemname", "colNumMin", "colNumMax", "assocLineName","statusString"};
        String[] columnTitles = {"Name", "Start Cdp", "End Cdp", "Assoc Name", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Lines");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }

	private void constructBeans()
	{
		surveyDefinitionsComboBox = new StsComboBoxFieldBean(this, "surveyDefinition", "Coor Source");
		surveyDefinitionsComboBox.setEditable(false);
        xOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "xOrigin", false, "X Origin:", true);
        yOriginBean = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "yOrigin", false, "Y Origin:", true);
        xMinBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "xMin", false, "X Min:", true);
        yMinBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "yMin", false, "Y Min:", true);
        xMaxBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "xMax", false, "X Max:", true);
        yMaxBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "yMax", false, "Y Max:", true);
        cdpIntervalBean = new StsFloatFieldBean(StsSeismicBoundingBox.class, "cdpInterval", false, "Cdp Interval:", true);
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

	// Survey definitions come from already loaded lines or user defined
	private void initialiseAvailableSurveyDefinitions()
	{
		surveyDefinitions.clear();

        StsSeismicBoundingBox unknownVolume = createUserSurveyVolume(volumes[0]);
        addSurveyDefinition(unknownVolume, "None");

        StsObject[] availableLines = wizard.getModel().getObjectList( StsPreStackLine2d.class);
		addSurveyDefinitions( availableLines, "2d prestack: ");

        availableLines = wizard.getModel().getObjectList( StsSeismicLine2d.class);
		addSurveyDefinitions( availableLines, "2d poststack: ");

//****wrw implemented 10/18/07 but feature turned off for now
//        StsSeismicBoundingBox userDefinitionVolume = createUserSurveyVolume(volumes[0]);
//        userSurveyDefinition = addSurveyDefinition(userDefinitionVolume, "User-defined line");

        defaultSurveyDefinition = (StsPostStackSurveyDefinitionPanel2d.SurveyDefinition)surveyDefinitions.get(0);
        currentSurveyDefinition = defaultSurveyDefinition;

	    surveyDefinitionsComboBox.setListItems(surveyDefinitions.toArray());
		surveyDefinitionsComboBox.setEditable(false);
		checkSetSurveyIsEditable();
	}

	private void addSurveyDefinitions(StsObject[] surveyObjects, String description)
	{
        for (int n = 0; n < surveyObjects.length; n++)
		{
            if( surveyObjects[n] != null)
                addSurveyDefinition(surveyObjects[n], description);
		}
	}

	private StsPostStackSurveyDefinitionPanel2d.SurveyDefinition addSurveyDefinition(StsObject surveyObject, String description)
	{
		StsPostStackSurveyDefinitionPanel2d.SurveyDefinition surveyDefinition =
                new StsPostStackSurveyDefinitionPanel2d.SurveyDefinition( surveyObject, description);
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

    private void buildPanel()
	{
        removeAll();

        surveyDefinitionBox = new StsGroupBox("Survey Definition");
	    surveyPanel = StsJPanel.addInsets();
	    parameterPanel = StsJPanel.addInsets();

        surveyPanel.gbc.fill = GridBagConstraints.NONE;
	    surveyDefinitionBox.gbc.fill = GridBagConstraints.NONE;
        surveyDefinitionBox.gbc.anchor = gbc.WEST;

		surveyPanel.add(surveyDefinitionsComboBox);
		surveyDefinitionBox.add(surveyPanel);

//****wrw implemented 10/18/07 but feature turned off for now
/*
        parameterPanel.addToRow(xOriginBean);
        parameterPanel.addEndRow(yOriginBean);
        parameterPanel.addToRow(xMinBean);
        parameterPanel.addEndRow(yMinBean);
        parameterPanel.addToRow(xMaxBean);
        parameterPanel.addEndRow(yMaxBean);
        parameterPanel.addEndRow(cdpIntervalBean);

        parameterPanel.setPanelObject(currentSurveyDefinition.boundingBox);
		surveyDefinitionBox.add(parameterPanel);
*/

        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.weighty = 1.0;
        add(volumeStatusTablePanel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        add(surveyDefinitionBox);

        volumeStatusTablePanel.addSelectRowNotifyListener(this);
        volumeStatusTablePanel.addTableModelListener(this);

        wizard.rebuild();
    }

    public StsPostStackSurveyDefinitionPanel2d.SurveyDefinition getSurveyDefinition()
    {
        return this.currentSurveyDefinition;
    }

    private boolean checkSurveyCdpRange(StsPostStackSurveyDefinitionPanel2d.SurveyDefinition surveyDefinition)
    {
        if( surveyDefinition.line instanceof StsSeismicLine)
        {
            StsSeismicBoundingBox selectedLine = wizard.getFirstSelectedVolume();
            StsSeismicLine assocLine = (StsSeismicLine)surveyDefinition.line;
            int nAssocCdp = assocLine.getNCols();
            int nSelectedCdp = selectedLine.getNCols();
            if( nAssocCdp < nSelectedCdp)
            {
                new StsMessage(wizard.frame, StsMessage.WARNING, "Number of association Cdps " + nAssocCdp +
                                                                 "\n is less than the number of input line Cdps " + nSelectedCdp);
                return false;
            }
        }
        return true;
    }

    private void setAssociation(StsPostStackSurveyDefinitionPanel2d.SurveyDefinition surveyDefinition)
    {
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if( selectedVolume == null) return;
        StsSegyLine2d selectedLine = (StsSegyLine2d)wizard.getSegyVolume(selectedVolume.getSegyFilename());
        if( surveyDefinition.line instanceof StsSeismicLine)
        {
            StsSeismicLine assocLine = (StsSeismicLine)surveyDefinition.line;
            selectedLine.assocLineName = assocLine.getName();
            selectedLine.setAssocLine(assocLine);
            selectedLine.setAssocLineName(assocLine.getName());
        }
        else if( surveyDefinition.toString() == "User-defined line")
        {
            selectedLine.assocLineName = "User-defined";
            StsSeismicBoundingBox seismicBoundingBox = (StsSeismicBoundingBox)surveyDefinition.line;
            seismicBoundingBox.setName( "User-defined line");
            selectedLine.setAssocLine( seismicBoundingBox);
            selectedLine.setAssocLineName( "User-defined line");
        }
        else if( surveyDefinition.toString().compareTo("None") == 0)
        {
            selectedLine.assocLineName = "Unknown";
        }
    }

	private void setSurveyDefinition(StsPostStackSurveyDefinitionPanel2d.SurveyDefinition surveyDefinition)
	{
        if(currentSurveyDefinition == surveyDefinition) return;
        if( !checkSurveyCdpRange ( surveyDefinition)) return;

        currentSurveyDefinition = surveyDefinition;
        setAssociation( surveyDefinition);
		checkSetSurveyIsEditable();

        parameterPanel.setPanelObject(surveyDefinition.boundingBox);
        parameterPanel.setEditable(isSurveyEditable);

        buildPanel();
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

	public void volumesSelected(StsSeismicBoundingBox[] volume)
	{
        if( volume == null) return;
        StsSeismicBoundingBox segyVolume = wizard.getFirstSegyVolume();
        if( segyVolume == null) return;
        
        wizard.saveSelectedAttributes( segyVolume);
        wizard.setSelectedVolumes(volume);

        setSurveyDefinition( defaultSurveyDefinition);
        surveyDefinitionsComboBox.setValueObject(currentSurveyDefinition);
        surveyDefinitionsComboBox.setEditable(true);
        StsSegyLine2d selectedLine = (StsSegyLine2d)segyVolume;

        buildPanel();
    }

	// This method is ultimately called by the trace analysis callback
	public void updatePanel()
	{
        if ((wizard.getSegyVolumesList().size() > 0) ||
			 (wizard.getSegyVolumesToProcess() != null) && (wizard.getSegyVolumesToProcess().length > 0))
		{
			wizard.enableNext();
		}
		else
		{
			wizard.disableNext();
		}
		volumeStatusTablePanel.replaceRows( wizard.getSegyVolumesList());
	}

    public void rowsSelected( int[] selectedIndices)
    {
        wizard.setSelectedVolumes( selectedIndices);
        volumesSelected( wizard.getSelectedVolumes());
    }

    public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
    {
        return currentSurveyDefinition.boundingBox;
    }

    class SurveyDefinition
	{
        String surveyNameString;
        StsSeismicBoundingBox boundingBox;
        StsObject line;

        SurveyDefinition(StsObject line, String description)
        {
            boundingBox = (StsSeismicBoundingBox)line;
            String name = boundingBox.getName();
            if(name == null)
                surveyNameString = description;
            else
                surveyNameString = new String(description + name);

            this.line = line;
        }

        SurveyDefinition()
        {
            boundingBox = new StsSeismicBoundingBox(false);
            surveyNameString = "None";
        }

	    public String toString() { return surveyNameString; }
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

package com.Sts.Actions.Wizards.HandVelocity;

import com.Sts.DB.DBCommand.*;
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

public class StsHandVelocitySurveyDefinitionPanel extends StsJPanel implements StsSelectRowNotifyListener, StsTableModelListener
{
	private StsHandVelocityWizard wizard;
    private StsHandVelocitySurveyDefinition wizardStep;
    private StsSeismicBoundingBox[] handVelVolumes = null;

    private StsGroupBox surveyDefinitionBox;
	private StsJPanel surveyPanel;

	private StsComboBoxFieldBean surveyDefinitionsComboBox;
    private JTextPane hintText = new JTextPane();

	public  ArrayList surveyDefinitions = new ArrayList();
	private StsHandVelocitySurveyDefinitionPanel.SurveyDefinition currentSurveyDefinition = null;
	private StsHandVelocitySurveyDefinitionPanel.SurveyDefinition defaultSurveyDefinition = null;
    private StsHandVelocitySurveyDefinitionPanel.SurveyDefinition userSurveyDefinition = null;

    private StsTablePanelNew volumeStatusTablePanel = null;

	public StsHandVelocitySurveyDefinitionPanel( StsHandVelocityWizard wizard, StsHandVelocitySurveyDefinition wizardStep)
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
        String[] columnNames = {"stemname", "rowNumMin", "rowNumMax","colNumMin", "colNumMax", "assocLineName","statusString"};
        String[] columnTitles = {"Name", "Start Inline", "End Inline", "Start Xline", "End Xline", "Assoc", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.setLabel("Handvels");
        volumeStatusTablePanel.initialize();
    }

	private void constructBeans()
	{
		surveyDefinitionsComboBox = new StsComboBoxFieldBean(this, "surveyDefinition", "Coor Source");
		surveyDefinitionsComboBox.setEditable(false);
    }

	public void initialize()
	{
        handVelVolumes = wizard.getHandVelVolumes();
        volumeStatusTablePanel.replaceRows( wizard.getHandVelVolumesList());
		initialiseAvailableSurveyDefinitions();
		rowsSelected(new int[] {0});
		setSurveyDefinition((StsHandVelocitySurveyDefinitionPanel.SurveyDefinition)surveyDefinitions.get(1));
		surveyDefinitionsComboBox.setSelectedItem(surveyDefinitions.get(1));
	    buildPanel();
	}

	// Survey definitions come from already loaded lines or volumes
	private void initialiseAvailableSurveyDefinitions()
	{
		surveyDefinitions.clear();

        StsSeismicBoundingBox unknownVolume = createUserSurveyVolume( handVelVolumes[0]);
        addSurveyDefinition(unknownVolume, "None");

        StsObject[] availableVolumes = wizard.getModel().getObjectList( StsPreStackLine2d.class);
		addSurveyDefinitions( availableVolumes, "2d prestack: ");

        availableVolumes = wizard.getModel().getObjectList( StsPreStackLineSet3d.class);
		addSurveyDefinitions( availableVolumes, "3d prestack: ");

        defaultSurveyDefinition = (StsHandVelocitySurveyDefinitionPanel.SurveyDefinition)surveyDefinitions.get(0);
        currentSurveyDefinition = defaultSurveyDefinition;

	    surveyDefinitionsComboBox.setListItems(surveyDefinitions.toArray());
		//surveyDefinitionsComboBox.setEditable(false);
		//surveyDefinitionsComboBox.setSelectedItem(surveyDefinitions.get(1));
	}

	private void addSurveyDefinitions( StsObject[] surveyObjects, String description)
	{
        for (int n = 0; n < surveyObjects.length; n++)
		{
            if( surveyObjects[n] != null)
                addSurveyDefinition(surveyObjects[n], description);
		}
	}

	private StsHandVelocitySurveyDefinitionPanel.SurveyDefinition addSurveyDefinition( StsObject surveyObject, String description)
	{
		StsHandVelocitySurveyDefinitionPanel.SurveyDefinition surveyDefinition =
                new StsHandVelocitySurveyDefinitionPanel.SurveyDefinition( surveyObject, description);
		surveyDefinitions.add(surveyDefinition);
		return surveyDefinition;
	}

	private StsSeismicBoundingBox createUserSurveyVolume( StsSeismicBoundingBox referenceVolume)
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

        surveyPanel.gbc.fill = GridBagConstraints.NONE;
	    surveyDefinitionBox.gbc.fill = GridBagConstraints.NONE;
        surveyDefinitionBox.gbc.anchor = gbc.WEST;

		surveyPanel.add( surveyDefinitionsComboBox);
		surveyDefinitionBox.add( surveyPanel);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.weighty = 1.0;
        add( volumeStatusTablePanel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        add( surveyDefinitionBox);

        hintText.setText( "Use the Velocity Analysis Wizard to extract handVels for the associated lines.");
        hintText.setEditable(false);
        hintText.setBackground(Color.lightGray);
        hintText.setFont(new Font("Dialog", 0, 10));
        add( hintText);

        volumeStatusTablePanel.addSelectRowNotifyListener(this);
        volumeStatusTablePanel.addTableModelListener(this);

        wizard.rebuild();
    }

    private boolean checkSurveyCdpRange( StsHandVelocitySurveyDefinitionPanel.SurveyDefinition surveyDefinition)
    {
        if( surveyDefinition.surveyNameString.compareTo("None") == 0) return true;

        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        int nSelectedRows = selectedVolume.getNRows();
        int nSelectedCols = selectedVolume.getNCols();

        if( surveyDefinition.volume instanceof StsPreStackLineSet3d)
        {
            int nAssocRows = 0;
            StsPreStackLineSet3d assocVolume = (StsPreStackLineSet3d)surveyDefinition.volume;
            nAssocRows = assocVolume.getNRows();
            if(nAssocRows > 0 &&  nAssocRows < nSelectedRows)
            {
                new StsMessage(wizard.frame, StsMessage.WARNING, "Cdp range of the associated volume, " + nAssocRows + ", " +
                                                 "\n is less than the cdp range of the handvel file, " + nSelectedRows + ".");
                return false;
            }
        }
        int nAssocCols = 0;
        if( surveyDefinition.volume instanceof StsPreStackLine2d)
        {
            StsPreStackLine2d assocVolume = (StsPreStackLine2d)surveyDefinition.volume;
            nAssocCols = assocVolume.getNCols();
        }
        else if( surveyDefinition.volume instanceof StsPreStackLineSet3d)
        {
            StsPreStackLineSet3d assocVolume = (StsPreStackLineSet3d)surveyDefinition.volume;
            nAssocCols = assocVolume.getNCols();
        }

        if( nAssocCols < nSelectedCols)
        {
            new StsMessage(wizard.frame, StsMessage.WARNING, "Cdp range of the associated volume, " + nAssocCols + ", " +
                                             "\n is less than the cdp range of the handvel file, " + nSelectedCols + ".");
            return false;
        }
        return true;
    }

    private void setAssociation( StsHandVelocitySurveyDefinitionPanel.SurveyDefinition surveyDefinition)
    {
        StsSeismicBoundingBox selectedVolume = wizard.getFirstSelectedVolume();
        if( selectedVolume == null) return;

        if( surveyDefinition.volume instanceof StsPreStackLine2d)
        {
            StsPreStackLine2d assocLine = (StsPreStackLine2d)surveyDefinition.volume;
            selectedVolume.assocLineName = assocLine.getName();
            selectedVolume.statusString = StsSeismicBoundingBox.STATUS_OK_STR;
            assocLine.setHandVelName( selectedVolume.getName());
            assocLine.setHandVelFiles( wizard.getHandVelFilesList() );

            StsChangeCmd cmd = new StsChangeCmd( assocLine, selectedVolume.getName(), "handVelName", false);
            wizard.getModel().addTransactionCmd("update StsPreStackLine2d handVelName", cmd);
           // cmd = new StsChangeCmd( assocLine, wizard.handVelSelect.panel.getCurrentDirectory(), "handVelDir", false);
           // wizard.getModel().addTransactionCmd("update StsPreStackLine2d handVelDir", cmd);
        }
        else if( surveyDefinition.volume instanceof StsPreStackLineSet3d)
        {
            StsPreStackLineSet3d assocVolume = (StsPreStackLineSet3d)surveyDefinition.volume;
            selectedVolume.assocLineName = assocVolume.getName();
            assocVolume.setHandVelName( selectedVolume.getName());
            assocVolume.setHandVelFiles( wizard.getHandVelFilesList() );

            StsChangeCmd cmd = new StsChangeCmd( assocVolume, selectedVolume.getName(), "handVelName", false);
            wizard.getModel().addTransactionCmd("update StsPreStackLineSet3d handVelName", cmd);
           // cmd = new StsChangeCmd( assocVolume, wizard.handVelSelect.panel.getCurrentDirectory(), "handVelDir", false);
           // wizard.getModel().addTransactionCmd("update StsPreStackLineSet3d handVelDir", cmd);
        }
        else if( surveyDefinition.toString().compareTo("None") == 0)
        {
            selectedVolume.assocLineName = "Unknown";
        }
    }

	private void setSurveyDefinition( StsHandVelocitySurveyDefinitionPanel.SurveyDefinition surveyDefinition)
	{
        if(currentSurveyDefinition == surveyDefinition) return;
        //TODO we need to smarten up checkSurveyCdpRange to deal with different spacings between the data and velocity volumes 
//        if( !checkSurveyCdpRange ( surveyDefinition)) return;

        currentSurveyDefinition = surveyDefinition;
        setAssociation( surveyDefinition);

        volumeStatusTablePanel.repaint();
    }

    public StsHandVelocitySurveyDefinitionPanel.SurveyDefinition getSurveyDefinition()
    {
        return this.currentSurveyDefinition;
    }

    public void removeRows( int firstRow, int lastRow)
    {
        if (wizard.getHandVelFilesList().size() > 0)
            wizard.enableNext();
        else
            wizard.disableNext();
    }

	public void removeVolume( StsSeismicBoundingBox volume){}

    public void volumesSelected( StsSeismicBoundingBox[] volume)
	{
        if( volume == null) return;

        setSurveyDefinition( defaultSurveyDefinition);
        surveyDefinitionsComboBox.setValueObject(currentSurveyDefinition);
        surveyDefinitionsComboBox.setEditable(true);
    }

    public void rowsSelected( int[] selectedIndices)
    {
        wizard.setSelectedVolumes( selectedIndices);
        volumesSelected( wizard.getSelectedVolumes());
    }

	public void updatePanel()
	{
        volumeStatusTablePanel.replaceRows( wizard.getHandVelVolumesList());
        if ((wizard.getHandVelVolumesList().size() > 0))
			wizard.enableNext();
		else
			wizard.disableNext();
	}

	class SurveyDefinition
	{
        String surveyNameString;
        StsSeismicBoundingBox boundingBox;
        StsObject volume;

        SurveyDefinition(StsObject volume, String description)
        {
            boundingBox = (StsSeismicBoundingBox)volume;
            String name = boundingBox.getName();
            if(name == null)
                surveyNameString = description;
            else
                surveyNameString = new String( description + name);

            this.volume = volume;
        }

	    public String toString() { return surveyNameString; }
	}

    public void cancel() { return;}
}
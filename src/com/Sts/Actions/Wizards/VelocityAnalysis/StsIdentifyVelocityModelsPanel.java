package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsIdentifyVelocityModelsPanel extends StsJPanel
{
    private StsVelocityAnalysisWizard wizard;
    private StsIdentifyVelocityModels wizardStep;

    private String velocityUnits = StsParameters.VEL_UNITS_NONE;
    double scaleMultiplier = 1.0f;
    public boolean createProfiles = false;
    public int numVelPts = 10;
    private ArrayList modelDefinitions = new ArrayList();
    private StsIdentifyVelocityModelsPanel.ModelDefinition currentModelDefinition;

    private StsGroupBox velocityGroupBox;
    private StsComboBoxFieldBean velocityModelBean;
    StsComboBoxFieldBean velocityUnitsBean;
    private StsBooleanFieldBean createProfilesBean;
    private StsIntFieldBean numVelPtsBean;

    public StsIdentifyVelocityModelsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsVelocityAnalysisWizard)wizard;
        this.wizardStep = (StsIdentifyVelocityModels)wizardStep;
        try
        {
             constructBeans();
             constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void constructBeans()
    {
        velocityGroupBox = new StsGroupBox("Input Velocity Model");
        velocityModelBean = new StsComboBoxFieldBean(this, "inputVelocityModelOrVolume", "Previous Velocity:");
        velocityModelBean.setToolTipText("Input Velocity Data or Model");
        velocityUnitsBean = new StsComboBoxFieldBean(this, "velocityUnits", "Velocity units:", StsParameters.VEL_UNITS );
        createProfilesBean = new StsBooleanFieldBean(this, "createProfiles", "Create Velocity Profile");
        createProfilesBean.setEditable(false);
        createProfilesBean.setEnabled(false);
        numVelPtsBean = new StsIntFieldBean(this, "numVelPts", 1,100, "Number of Velocity Points", true);
        numVelPtsBean.setEditable(false);
        numVelPtsBean.setEnabled(false);
    }

    public void constructPanel()
    {
        velocityGroupBox.add(velocityModelBean);
        velocityGroupBox.add(velocityUnitsBean);
        velocityGroupBox.add(createProfilesBean);
        velocityGroupBox.add(numVelPtsBean);
        add(velocityGroupBox);
    }

    public void initialize()
    {
        setModelDefinitions();
        setVolumes();
    }

    private StsIdentifyVelocityModelsPanel.ModelDefinition getModelDefinition( Object modelType)
    {
        StsIdentifyVelocityModelsPanel.ModelDefinition modelDefinition = null;
        if( modelType instanceof StsSeismicLineSet)
            modelDefinition = new StsIdentifyVelocityModelsPanel.ModelDefinition( modelType, "lineset:");
        else if( modelType instanceof StsSeismicVolume)
            modelDefinition = new StsIdentifyVelocityModelsPanel.ModelDefinition( modelType, "velocity volume:");
        else if( modelType instanceof StsPreStackVelocityModel)
        {
            StsPreStackVelocityModel velocityModel = (StsPreStackVelocityModel)modelType;
            if( velocityModel.getNProfiles() > 0)
                modelDefinition = new StsIdentifyVelocityModelsPanel.ModelDefinition( modelType, "velocity model:");
        }
        else
            modelDefinition = new StsIdentifyVelocityModelsPanel.ModelDefinition( modelType, "");

        return modelDefinition;
    }

    private void setModelDefinitions()
    {
        modelDefinitions.clear();
        
        StsPreStackLineSet volume = wizard.getPreStackVolume();
        Object[] modelTypes = volume.getAvailableModelsAndVolumesList();
        for( int i=0; i<modelTypes.length; i++)
        {
            StsIdentifyVelocityModelsPanel.ModelDefinition modelDefinition;
            modelDefinition = getModelDefinition( modelTypes[i]);
            if( modelDefinition != null)
                modelDefinitions.add( modelDefinition);
        }
    }

    public StsIdentifyVelocityModelsPanel.ModelDefinition getInputVelocityModelOrVolume()
    {
        Object modelType = wizard.getInputVelocityModelOrVolume();
        currentModelDefinition = getModelDefinition( modelType);
        return currentModelDefinition;
    }

    public void setInputVelocityModelOrVolume( StsIdentifyVelocityModelsPanel.ModelDefinition modelDefinition)
    {
        if( modelDefinition.model instanceof StsSeismicBoundingBox)
        {
            createProfilesBean.setEditable(true);
            numVelPtsBean.setEditable(true);
            estimateVelocityUnits((StsSeismicBoundingBox)modelDefinition.model);
        }
        velocityModelBean.setSelectedItem( modelDefinition);
        wizard.setInputVelocityModelOrVolume( modelDefinition.model);
        currentModelDefinition = modelDefinition;
    }

    private void setVolumes()
    {
        StsPreStackLineSet volume = wizard.getPreStackVolume();
        Object modelType = volume.getInputVelocityModel();
        StsIdentifyVelocityModelsPanel.ModelDefinition modelDefinition = getModelDefinition( modelType);
        
        velocityModelBean.setListItems( modelDefinitions.toArray());
        velocityModelBean.setSelectedItem( modelDefinition);
    }

    private void estimateVelocityUnits( StsSeismicBoundingBox seismicVolume)
    {
        velocityUnits = StsParameters.estimateVelocityUnits(seismicVolume.dataMin);
        velocityUnitsBean.setValueObject(velocityUnits);
    }

    public void setVelocityUnits(String unitsString)
    {
        velocityUnits = unitsString;
        setVelScaleMultiplier();
    }

    public String getVelocityUnits() { return velocityUnits; }

    private void setVelScaleMultiplier()
    {
        StsProject project = wizard.getModel().getProject();
        scaleMultiplier = project.calculateVelScaleMultiplier(velocityUnits);
    }

    public double getVelScaleMultiplier() { return scaleMultiplier; }

    public void setCreateProfiles(boolean isSelected) {createProfiles = isSelected;}
    public boolean getCreateProfiles() {return createProfiles;}
    public void setNumVelPts(int numVelPts){ this.numVelPts = numVelPts;}
    public int getNumVelPts() { return numVelPts; }

    class ModelDefinition
    {
        String modelNameString;
        Object model;

        ModelDefinition( Object model, String description)
        {
            if( model.toString() != null)
                modelNameString = new String( description + model.toString());
            else
                modelNameString = new String( description);
            this.model = model;
        }

        public String toString() { return modelNameString; }
    }
}
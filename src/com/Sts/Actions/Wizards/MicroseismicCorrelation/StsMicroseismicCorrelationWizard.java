package com.Sts.Actions.Wizards.MicroseismicCorrelation;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard;
import com.Sts.MVC.*;
import com.Sts.DBTypes.*;

public class StsMicroseismicCorrelationWizard extends StsWizard
{
    private StsComputeSetup computeSetup;
    private StsSelectSensors selectSensors;
    private StsProcessAttribute processAttribute;

    String name = "Attribute";
    boolean floatVolume = false;
    public Object[] selectedSensors = new Object[0];

    private StsWizardStep[] mySteps =
    {
        computeSetup = new StsComputeSetup(this),
        selectSensors = new StsSelectSensors(this),
        processAttribute = new StsProcessAttribute(this)
    };

    public StsMicroseismicCorrelationWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 400);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Compute Correlation Attributes");
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public Object getSelectedSensors()
    {
        return selectedSensors;
    }

    public void setSelectedSensors(Object sensor)
    {
        if(sensor == null) return;
        selectedSensors = selectSensors.getSelectedSensors();
    }

    public StsSeismicVolume  getSelectedVolume()
    {
        return computeSetup.panel.getSelectedVolume();
    }

    public String getAttributeName()
    {
        return StsComputeSetupPanel.OPERATORS[getOperator()] + "-" + getSelectedVolume().getName().substring(0,10);
    }
    
    public int getOperator()
    {
        return computeSetup.panel.getOperator();
    }

	public boolean checkPrerequisites()
	{
		StsDynamicSensorClass sensorClass = (StsDynamicSensorClass) model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
		StsObject[] sensors = sensorClass.getSensors();
		if (sensors.length == 0) {
			reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
			return false;
		}
        int nSeismicVolumes = getSeismicVolumes().length;
        nSeismicVolumes = nSeismicVolumes + getVirtualVolumes().length;
		if(nSeismicVolumes < 1)
        {
			reasonForFailure = "No virtual or seismic volumes exist. Must have at least one seismic volume.";
			return false;
		}
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
        return true;
	}

    public StsSeismicVolume[] getSeismicVolumes()
    {
        return (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
    }

    public StsVirtualVolume[] getVirtualVolumes()
    {
        StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class);
        return vvClass.getVirtualVolumes();
    }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsMicroseismicCorrelationWizard vvWizard = new StsMicroseismicCorrelationWizard(actionManager);
        vvWizard.start();
    }
}
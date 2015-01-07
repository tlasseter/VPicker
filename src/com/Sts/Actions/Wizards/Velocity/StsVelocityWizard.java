package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsVelocityWizard extends StsWizard
{
	private StsEditVelocity editVelocity = new StsEditVelocity(this);
    private StsVolumeDefine volumeDefine = new StsVolumeDefine(this);
    //private StsFromSonics fromSonics = new StsFromSonics(this);
    private StsAssignVelocity assignVelocity = new StsAssignVelocity(this);
	private StsCreateVelocityModel createVelocity = new StsCreateVelocityModel(this);
	private StsWizardStep[] wizardSteps = new StsWizardStep[]{ editVelocity, volumeDefine, assignVelocity, createVelocity };

	private StsModelSurface[] selectedSurfaces = null;
    private StsSeismicVolume inputVelocityVolume = null;
    private String velocityTypeString = V_AVG_STRING;

    private float topDepthDatum = 0.0f;
    private float topTimeDatum = 0.0f;
    private float minVelocity;
    private float maxVelocity;
    static final public String V_AVG_STRING = StsParameters.V_AVG_STRING;
    static final public String V_INSTANT_STRING = StsParameters.V_INSTANT_STRING;
    static final public String[] VELOCITY_TYPE_STRINGS = new String[] { V_AVG_STRING, V_INSTANT_STRING };

	public StsVelocityWizard(StsActionManager actionManager)
	{
		super(actionManager);
		addSteps(wizardSteps);
        disableFinish();
    }

	public void previous()
	{
        if(currentStep == assignVelocity)
        {
            gotoStep(editVelocity);
        }
        else
            gotoPreviousStep();
	}

	public void next()
	{
        if(currentStep == editVelocity)
        {
 //           if(editVelocity.panel.getFromSonics() == true)
 //               gotoStep(fromSonics);
 //           else
            {
                if(editVelocity.panel.getVolume() == null)
                {
                    if(!editVelocity.panel.useWellControl)
                        gotoStep(assignVelocity);
                    else
                        gotoStep(createVelocity);
                }
                else
                {
                    gotoStep(volumeDefine);
                }
            }
        }
        else if(currentStep == volumeDefine)
        {
            gotoStep(createVelocity);
        }
//        else if(currentStep == fromSonics)
//            gotoStep(createVelocity);
        else
            gotoNextStep();
    }

	public boolean start()
	 {
		 dialog.setTitle("Create/Edit Velocity Model");
		 return super.start();
	 }

	 public boolean end()
	 {
         model.setActionStatus(getClass().getName(), StsModel.STARTED);
		 return super.end();
    }

	public void setSurfaces(StsModelSurface[] surfaces)
	{
		this.selectedSurfaces = surfaces;
	}

	public StsModelSurface[] getSelectedSurfaces() { return selectedSurfaces; }

    public void setVelocityVolume(StsSeismicVolume volume)
    {
        this.inputVelocityVolume = volume;
		byte type = StsParameters.getVelocityTypeFromString(velocityTypeString);
        if(volume != null) inputVelocityVolume.setType(type);
    }

    public StsSeismicVolume getVelocityVolume() { return inputVelocityVolume; }

    public void setVelocityTypeString(String typeString)
    {
		velocityTypeString = typeString;
		byte type = StsParameters.getVelocityTypeFromString(typeString);
		if(inputVelocityVolume != null) inputVelocityVolume.setType(type);
    }

    public String getVelocityTypeString() { return velocityTypeString; }

    public void setTopDepthDatum(float depth) { this.topDepthDatum = depth;  }
    public float getTopDepthDatum() { return topDepthDatum; }
    public void setTopTimeDatum(float time) { this.topTimeDatum = time; }
    public float getTopTimeDatum() { return topTimeDatum; }

    public boolean constructVelocityModel(StsProgressPanel panel)
    {
		boolean useWellTD = false;
        double scaleMultiplier = volumeDefine.panel.getVelScaleMultiplier();
		float newTimeInc = assignVelocity.panel.getZInc();
		float [] oneWayConstantIntervalVelocities =  assignVelocity.panel.getIntervalVelocities();
		float oneWayConstantBottomIntervalVelocity = assignVelocity.panel.btmVelocity;

        boolean useWellControl = editVelocity.panel.getWellControl();

		float markerFactor = editVelocity.panel.getMarkerFactor();
        int gridType =          editVelocity.panel.getGridType();
        // Convert to m/msec
	 /*
        if(intervalVelocities != null)
        {
            for (int i = 0; i < intervalVelocities.length; i++)
                intervalVelocities[i] = intervalVelocities[i] / 1000.0f;
        }
    */
        StsSeismicVelocityModel velocityModel = model.getProject().constructVelocityModel(selectedSurfaces, inputVelocityVolume,
            topTimeDatum, topDepthDatum, minVelocity, maxVelocity, scaleMultiplier, newTimeInc,
				oneWayConstantIntervalVelocities, oneWayConstantBottomIntervalVelocity, useWellControl, markerFactor, gridType, panel);

        if(velocityModel != null) return true;

        new StsMessage(model.win3d, StsMessage.WARNING, "Velocity wizard failed to build velocity model.\n" +
                       "Change parameters and try again or cancel.");
        return false;
    }

    public void setMinVelocity(float minVelocity) { this.minVelocity = minVelocity;  }
    public float getMinVelocity() { return minVelocity; }
    public void setMaxVelocity(float maxVelocity) { this.maxVelocity = maxVelocity;  }
    public float getMaxVelocity() { return maxVelocity; }
}


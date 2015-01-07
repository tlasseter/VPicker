package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.MVC.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsParameters;
import com.Sts.DBTypes.StsSeismicVelocityModel;

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
public class StsSeisVelWizard extends StsWizard
{
	//private StsEditVelocity editVelocity = new StsEditVelocity(this);
    private StsVolumeDefineSV volumeDefine = new StsVolumeDefineSV(this);
    private StsFromSonics fromSonics = new StsFromSonics(this);
    //private StsAssignVelocity assignVelocity = new StsAssignVelocity(this);
	private StsCreateVelocityVolume createVelocity = new StsCreateVelocityVolume(this);
	private StsWizardStep[] wizardSteps = new StsWizardStep[]{ volumeDefine, fromSonics , createVelocity};

	//private StsModelSurface[] selectedSurfaces = null;
    private StsSeismicVolume inputVelocityVolume = null;
    private String velocityTypeString = V_AVG_STRING;

    private float topDepthDatum = 0.0f;
    private float topTimeDatum = 0.0f;
    private float minVelocity = 1000;
    private float maxVelocity = 15000;
	private boolean useSonic = true;
	private boolean useVelf = false;
	private String[] velfList = null;
    static final public String V_AVG_STRING = StsParameters.V_AVG_STRING;
    static final public String V_INSTANT_STRING = StsParameters.V_INSTANT_STRING;
    static final public String[] VELOCITY_TYPE_STRINGS = new String[] { V_AVG_STRING, V_INSTANT_STRING };

	public StsSeisVelWizard(StsActionManager actionManager)
	{
		super(actionManager);
		addSteps(wizardSteps);
        disableFinish();
    }

	public void previous()
	{
        if(currentStep == fromSonics)
        {
            gotoStep(volumeDefine);
        }
		if(currentStep == createVelocity)
        {
	       gotoStep(fromSonics);
        }

        else
            gotoPreviousStep();
	}

	public void next()
	{
        if(currentStep == volumeDefine)
        {
                gotoStep(fromSonics);
        }
        else
            gotoNextStep();
    }

	public boolean start()
	 {
		 dialog.setTitle("Create Initial Seismic Velocity Volume");
		 return super.start();
	 }

	 public boolean end()
	 {
         model.setActionStatus(getClass().getName(), StsModel.STARTED);
		 return super.end();
    }


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
		//float newTimeInc = assignVelocity.panel.getZInc();
		float [] intervalVelocities = null;
//		if(editVelocity.panel.getFromSonics() == true)
//		   //useWellTD=true;
//	   inputVelocityVolume =
//	    else
//			intervalVelocities = assignVelocity.panel.getIntervalVelocities();

//		boolean useWellControl = editVelocity.panel.getWellControl();
		// Convert to m/msec
//		if(intervalVelocities != null)
//		{
//			for (int i = 0; i < intervalVelocities.length; i++)
//				intervalVelocities[i] = intervalVelocities[i] / 1000.0f;
//		}
		double newTimeInc=0.04;
		StsSeismicVelocityModel velocityModel = model.getProject().constructVelocityModelSV(topTimeDatum, topDepthDatum, minVelocity, maxVelocity, scaleMultiplier, newTimeInc,panel);

		if(velocityModel != null) return true;

		new StsMessage(model.win3d, StsMessage.WARNING, "Velocity wizard failed to build velocity model.\n" +
					   "Change parameters and try again or cancel.");
		return false;
	}

	public boolean constructVelocityVolume(StsProgressPanel panel)
	{
		double scaleMultiplier = volumeDefine.panel.getVelScaleMultiplier();
	    try
        {
            if(constructVelocityModel(panel))
            {
                inputVelocityVolume = model.getProject().constructVelocityVolume(inputVelocityVolume, topTimeDatum, topDepthDatum,
                        minVelocity, maxVelocity, scaleMultiplier, useSonic, useVelf, velfList, panel);
            }
            else
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Velocity wizard failed to build velocity model.\n" +
                    "Change parameters and try again or cancel.");
                return false;
            }
	    }
        catch (Exception e)
        {
		   StsException.outputWarningException(this, "constructVelocityVolume", e);
	    }

	    if(inputVelocityVolume != null) return true;
		new StsMessage(model.win3d, StsMessage.WARNING, "Velocity wizard failed to build velocity volume.\n" +
				"Change parameters and try again or cancel.");
        return false;
    }

    public void setMinVelocity(float minVelocity) { this.minVelocity = minVelocity;  }
    public float getMinVelocity() { return minVelocity; }
	public void setMaxVelocity(float maxVelocity) { this.maxVelocity = maxVelocity;  }
    public float getMaxVelocity() { return maxVelocity; }
	public void setUseSonic(boolean val) { this.useSonic = val;  }
    public boolean getUseSonic() { return useSonic; }
	public void setUseVelf(boolean val) { this.useVelf = val;  }
    public boolean getUseVelf() { return useVelf; }
	public void setVelfList(String[] val) { this.velfList = val; }
	public String[] getVelfList() { return velfList; }
}


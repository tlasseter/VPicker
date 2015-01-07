package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsSemblanceRangeProperties extends StsSeismicPanelProperties implements StsSerializable
{
    /** percent range for residual semblance panel: multiplier ranges from 1-range/100 to 1+range/100 */
    public float percentRange = 30.0f;
    /** percent range icnrement for residual semblance panel; number of columns is 2*range/rangeInc */
    public float percentRangeInc = 0.5f;
	public float zMin;
	public float zMax;
	public float zInc;
    public float zMinLimit;
    public float zMaxLimit;
    /** percent increment between panels on a VVS display */
    public float vvsPercentRangeInc = 3.0f;  //1% change is too small SWC 6/3/09
    public float velPickThreshold = 50.0f; //minimum time between velocity picks (ms)
    public boolean autoSaturateColors = true;

//	public StsColorscale semblanceColorscale = null;

    transient public float velocityMinLimit = 0.0f;
	transient public float velocityMaxLimit = 30.0f; //sometimes velocities get faster than 25 SWC 6/3/09
	transient public float velocityIncMinLimit = 0.001f;
	transient public float velocityIncMaxLimit = 5.0f;
	transient public float percentRangeIncMinLimit;
	transient public float percentRangeIncMaxLimit;

	transient StsBooleanFieldBean autoSaturateBean;
    transient StsFloatFieldBean velocityMinBean;
	transient StsFloatFieldBean velocityMaxBean;
	transient StsFloatFieldBean velocityIncBean;
	transient StsFloatFieldBean zMinBean;
	transient StsFloatFieldBean zMaxBean;
	// transient StsFloatFieldBean zIncBean;  //don't need this - zinc is always = sample rate of data
	transient StsFloatFieldBean velPickThresholdBean;
	transient StsFloatFieldBean percentRangeBean;
	transient StsFloatFieldBean percentRangeIncBean;
    transient StsFloatFieldBean vvsPercentRangeIncBean;
    transient StsComboBoxFieldBean velocityTypeBean;

    static final String title = "Range Properties";

	public StsSemblanceRangeProperties()
	{
	}

	public StsSemblanceRangeProperties(StsModel model, StsPreStackLineSet3dClass seismicClass, String fieldName)
	{
		super(title, fieldName);
	}

	public StsSemblanceRangeProperties(StsObject parentObject, StsModel model, String fieldName)
	{
		super(parentObject, null, title, fieldName);
//		initializeColorscale(model);
		initializeSemblanceZRange(model);
		initializeRangeLimits(model);
	}

	private void initializeRangeLimits(StsModel model)
	{
		StsProject project = model.getProject();
		if(project.velocityUnitsChanged(true))
		{
			velocityMin = project.convertVelocity(velocityMin, model);
			velocityMax = project.convertVelocity(velocityMax, model);
			velocityStep = project.convertVelocity(velocityStep, model);
			velocityMinLimit = project.convertVelocity(velocityMinLimit, model);
			velocityMaxLimit = project.convertVelocity(velocityMaxLimit, model);
			velocityIncMinLimit = project.convertVelocity(velocityIncMinLimit, model);
			velocityIncMaxLimit = project.convertVelocity(velocityIncMaxLimit, model);
		}
//        StsPreStackLineSet lineSet = (StsPreStackLineSet)parentObject;
//        zMinLimit = lineSet.zMin;
//        zMaxLimit = lineSet.zMax;
        percentRangeIncMinLimit = percentRange/1000;
		percentRangeIncMaxLimit = percentRange/10;
	}

	public void initializeSemblanceZRange(StsModel model)
	{
		StsProject project = model.getProject();
		zMin = project.getZorTMin();
		zMax = project.getZorTMax();
		zInc = project.getZorTInc();  //changed, we now want semblance at every sample default!! SWC
        zMinLimit = zMin;
        zMaxLimit = zMax;
    }


	public void initializeBeans()
	{
		StsModel model = StsObject.getCurrentModel();
		//initializeRangeLimits(model);
		if(propertyBeans == null)
		{
//            ipsBean = new StsFloatFieldBean(this, "inchesPerSecond", .1f, 50.0f, "Inches per Second:");
		    autoSaturateBean = new StsBooleanFieldBean(this, "autoSaturateColors", "Auto-Saturate Velocity Colors:");
		    velocityMinBean = new StsFloatFieldBean(this, "velocityMin", velocityMinLimit * 1000, velocityMaxLimit * 1000, "Velocity Min. (ft-m/sec):", true);
		    velocityMinBean.fixStep(velocityStep * 1000);
		    velocityMaxBean = new StsFloatFieldBean(this, "velocityMax", velocityMinLimit * 1000, velocityMaxLimit * 1000, "Velocity Max. (ft-m/sec):", true);
		    velocityMaxBean.fixStep(velocityStep * 1000);
		    velocityIncBean = new StsFloatFieldBean(this, "velocityStep", velocityIncMinLimit * 10 * 1000, velocityIncMaxLimit * 1000, "Velocity Inc. (ft-m/sec):", true);
		    velocityIncBean.fixStep(10 * velocityIncMinLimit * 1000);
		    percentRangeBean = new StsFloatFieldBean(this, "percentRange", 0.0f, 100.0f, "Residual Plot, Percent Range:", true);
		    percentRangeIncBean = new StsFloatFieldBean(this, "percentRangeInc", percentRange/1000, percentRange/10, "Percent Range Inc:", true);
		    vvsPercentRangeIncBean = new StsFloatFieldBean(this, "vvsPercentRangeInc", 0.001f, 99f, "VVS Panel % Inc:", true);
		    vvsPercentRangeIncBean.setToolTipText("Percent difference in velocity between VVS panels");
	  //    zMinBean = new StsFloatFieldBean(this, "zMin", zMinLimit, zMaxLimit, "Time/Depth Min:", true);   //not necessary, always use full time range while picking velocities SWC 8/26/09
      //	zMaxBean = new StsFloatFieldBean(this, "zMax", zMinLimit, zMaxLimit, "Time/Depth Max:", true);
      //	zIncBean = new StsFloatFieldBean(this, "zInc", true, "Time/Depth Inc:");
	  //    zIncBean.setRangeFixStep( model.project.getZorTMin(), model.project.getZorTMax(), model.project.getZorTInc());
		    velPickThresholdBean = new StsFloatFieldBean(this, "velPickThreshold", 0f, 3000f, "Velocity Pick Threshold (ms)", true);
		    velPickThresholdBean.setToolTipText("Minimum time between velocity picks (ms)");
		    velocityTypeBean = new StsComboBoxFieldBean(this, "velocityTypeString", "Velocity Display Type:", StsParameters.VEL_STRINGS);
		    propertyBeans = new StsFieldBean[]
		   {
	      //        ipsBean,
		            autoSaturateBean,
		            velocityTypeBean,
		            velocityMinBean,
		            velocityMaxBean,
		            velocityIncBean,
		            percentRangeBean,
		            percentRangeIncBean,
		            vvsPercentRangeIncBean,
		            velPickThresholdBean,
		            zMinBean,
		            zMaxBean,
		         //	zIncBean
		    };
		}
	}

	public float getZMin() { return zMin; }
	public float getZMax() { return zMax; }
	public float getZInc() { return zInc; }
	public float getPercentRange() { return percentRange; }
	public float getPercentRangeInc() { return percentRangeInc; }
    public float getVvsPercentRangeInc() { return vvsPercentRangeInc; }

    public void setZMin(float min) { zMin = min; setRangeChanged(true); }
	public void setZMax(float max) { zMax = max; setRangeChanged(true); }
	public void setZInc(float inc) { zInc = inc; setRangeChanged(true); }
	public void setPercentRange(float percent)
	{
		percentRange = percent;
		percentRangeIncMinLimit = percentRange/1000;
		percentRangeIncMaxLimit = percentRange/10;
		percentRangeIncBean.setRange(percentRangeIncMinLimit, percentRangeIncMaxLimit);
        setRangeChanged(true);
    }

    public void setPercentRangeInc(float percentInc) 
    { 
        percentRangeInc = percentInc; 
        setRangeChanged(true);
    }
    public void setVvsPercentRangeInc(float percentInc) 
    { 
        vvsPercentRangeInc = percentInc; 
        setRangeChanged(true);
    }

/*
    public void setInchesPerSecond(float inchesPerSecond, boolean setval)
    {
        if(this.inchesPerSecond == inchesPerSecond) return;
        if(setval) super.setInchesPerSecond(inchesPerSecond);
        if(ipsBean != null) ipsBean.setValue(inchesPerSecond);
    }
*/
/*
	public boolean getSetTracesRescaleRequired()
	{
		StsSemblanceRangeProperties savedProperties = (StsSemblanceRangeProperties)this.savedProperties;
		rangeChanged = 	zMin != savedProperties.zMin ||
						zMax != savedProperties.zMax ||
						zInc != savedProperties.zInc ||
						velocityMin != savedProperties.velocityMin ||
						velocityMax != savedProperties.velocityMax ||
						velocityInc != savedProperties.velocityInc;

		recompute = velocityInc != savedProperties.velocityInc ||
					zInc != savedProperties.zInc ||
					windowWidth != savedProperties.windowWidth ||
					velocityMin < savedProperties.velocityMin ||
					velocityMax > savedProperties.velocityMax ||
					semblanceType != savedProperties.semblanceType;
		return rangeChanged || recompute;
	}
*/
/*
	public void saveState()
	{
		savedProperties = (StsSemblanceRangeProperties)StsToolkit.copyObjectNonTransients(this);
	}
*/
/*
    public void displayDialog()
	{
		StsToolkit.createDialog(getPanel());
	}
*/
    
    public void setVelocityTypeString(String typeString)
    {
        StsPreStackVelocityModel velModel = StsPreStackLineSetClass.currentProjectPreStackLineSet.velocityModel;
        if (velModel == null) return;
        byte velocityType = StsParameters.getVelocityTypeFromString(typeString);
        if(velModel.getVelocityType() == velocityType) return;
        velModel.setVelocityType(velocityType);
        currentModel.viewObjectChangedAndRepaint(this, velModel);
    }
    
    public String getVelocityTypeString()
    {
        StsPreStackVelocityModel velModel = StsPreStackLineSetClass.currentProjectPreStackLineSet.velocityModel;
        if (velModel == null) return null;
        return velModel.getVelocityTypeString();
    }
    
    public boolean isAutoSaturateColors()
    {
        return autoSaturateColors;
    }

    public void setAutoSaturateColors(boolean autoSaturateColors)
    {
        this.autoSaturateColors = autoSaturateColors;
    }
    
    public float getVelPickThreshold() { return velPickThreshold; }
	public void setVelPickThreshold(float velPickThreshold) 
	{ 
		this.velPickThreshold = velPickThreshold;
	}
}

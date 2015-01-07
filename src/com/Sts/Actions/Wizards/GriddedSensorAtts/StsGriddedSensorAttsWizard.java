package com.Sts.Actions.Wizards.GriddedSensorAtts;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.Types.StsPoint;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Types.StsSeismicBoundingBox;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;

public class StsGriddedSensorAttsWizard extends StsWizard
{
    private StsSelectSensors selectSensors;
    private StsVolumeDefinition defineVolume;
    private StsSensorVVolume sensorVirtualVolume;
    private StsGriddedSensorAttsProcess processVolume;

    String volName = "SenVolName";
    StsSeismicBoundingBox seismicBoundingBox = null;
    public Object[] selectedSensors = new Object[0];
    StsSeismicVolume[] volumes = null;
    public StsSensorVirtualVolume sensorVolume = null;

    private StsWizardStep[] mySteps =
    {
        selectSensors = new StsSelectSensors(this),
        defineVolume = new StsVolumeDefinition(this),    
        sensorVirtualVolume = new StsSensorVVolume(this),
        processVolume = new StsGriddedSensorAttsProcess(this)
    };

    public StsGriddedSensorAttsWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 600);
        addSteps(mySteps);

        model.getProject().setZDomainString(StsParameters.TD_DEPTH_STRING);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        StsVirtualVolume[] vvolumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes();
        volumes = (StsSeismicVolume[]) StsMath.arrayAddArray(volumes, vvolumes);
        if(volumes ==  null || volumes.length == 0)
        {
            if(!StsYesNoDialog.questionValue(frame,"No seismic volumes exist, if you plan on loading any\n it is recommended you do so prior to generating gridded attributes.\n\n Do you want to continue?"))
                return false;
        }
        
        dialog.setTitle("Define Gridded Sensor Attributes");
        if(!super.start()) return false;

        checkAddToolbar();

        return true;
    }

   	public boolean checkPrerequisites()
	{
		StsDynamicSensorClass sensorClass = (StsDynamicSensorClass) model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
		StsObject[] sensors = sensorClass.getSensors();
		if (sensors.length == 0) {
			reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
			return false;
		}
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
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
        if(currentStep == selectSensors)
        {
            if(volumes.length == 0)
            {
                new StsMessage(this.frame,StsMessage.WARNING,"No seismic volumes exist, will need to define grid dimensions.\n");
                // Create a seismic volume from user input
                gotoStep(defineVolume);
            }
            else
            {
                seismicBoundingBox = volumes[0];
    			gotoStep(sensorVirtualVolume);
            }
        }
        else
    	{
    	    gotoNextStep();
    	}
    }

    public void setVolumeName(String name) { this.volName = name; }
    public String getVolumeName() { return volName; }

    public StsSeismicVolume[] getSeismicVolumes()
    {
        return (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
    }

    public StsVirtualVolume[] getVirtualVolumes()
    {
        StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class);
        return vvClass.getVirtualVolumes();
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

    public boolean buildVolume(StsProgressPanel ppanel)
    {
        String name = null;
        Object[] sensors = (Object[])getSelectedSensors();
        StsTimeCurve selectedCurve = getAttributeCurve();
        byte accumType = getAccumType();
        float xyOffset = getXyOffset();
        float zOffset = getZOffset();
        byte shapeType = getShapeType();
        float zStep = getZStep();

        //
        // Construct a name if one is not supplied
        //
        name = getVolumeName();
        if((name == null) || (name.length() <= 0))
        {
            name = "SenVolName";
        }
        ppanel.appendLine("Constructing sensor virtual volume:" + name);

        sensorVolume = StsSensorVirtualVolume.constructor(name, sensors, ppanel, getBoundingBox(), selectedCurve, accumType, zStep, xyOffset, zOffset, shapeType);
        if (sensorVolume == null)
        {
            new StsMessage(frame, StsMessage.ERROR, "Failed to construct virtual volume.\n");
            return false;
        }
        ppanel.appendLine("Successfully created sensor virtual volume:" + name);

        sensorVolume.addToProject(false);
        sensorVolume.addToModel();

        model.viewObjectChangedAndRepaint(this, sensorVolume);

        ppanel.appendLine("Press Finish>> Button to complete virtual volume creation.");
        return true;
    }

    public StsSeismicBoundingBox getBoundingBox()
    {
        return seismicBoundingBox;
    }

    public void setBoundingBox(StsSeismicBoundingBox bb)
    {
        seismicBoundingBox = bb;
    }

    public StsTimeCurve getAttributeCurve()
    {
        return sensorVirtualVolume.getAttribute();
    }
    public byte getAccumType()
    {
        return sensorVirtualVolume.getAccumType();
    }
    public float getXyOffset()
    {
        return sensorVirtualVolume.getXyOffset();
    }
    public float getZOffset()
    {
        return sensorVirtualVolume.getZOffset();
    }
    public float getZStep()
    {
        return sensorVirtualVolume.getZStep();
    }
    public byte getShapeType()
    {
        return sensorVirtualVolume.getShapeType();
    }
}

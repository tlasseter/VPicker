package com.Sts.Actions.Wizards.Alarms;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.Sounds.StsSound;

public class StsAlarmsWizard extends StsWizard
{
    private StsSelectSensors selectSensors;
    private StsDefineSurfaceAlarm surfaceDefine;
    private StsDefineWellAlarm wellDefine;
    private StsDefineValueAlarm valueDefine;

    transient String name = "Alarm";
    transient byte alarmType = StsAlarm.WELL;
    transient byte queueType = StsAlarm.AUDIO;
    transient String soundFile = StsSound.BUZZ;
    transient float offset = 100.0f;
    transient float threshold = 5.0f; // Percentage
    transient float value = 0.0f;
    transient boolean inside = true;
    transient Object[] selectedSensors = new Object[0];
    transient StsSurface selectedSurface = null;
    transient byte refDirection = StsSurfaceAlarm.ABOVE;
    transient StsWell selectedWell = null;

    private StsWizardStep[] mySteps =
    {
        selectSensors = new StsSelectSensors(this),
        surfaceDefine = new StsDefineSurfaceAlarm(this),
        wellDefine = new StsDefineWellAlarm(this),
        valueDefine = new StsDefineValueAlarm(this)
    };

    public StsAlarmsWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 400);
        addSteps(mySteps);
        setHelpSet();
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Alarms");
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        StsAlarm alarm = null;

        // Build Defined Alarm
        if(alarmType == StsAlarm.SURFACE)
        {
            if(selectedSurface != null)
                alarm = new StsSurfaceAlarm(soundFile, selectedSurface, offset, refDirection);
        }
        else if(alarmType == StsAlarm.WELL)
        {
            if(selectedWell != null)
                alarm = new StsWellAlarm(soundFile, selectedWell, offset, inside);
        }
        else if(alarmType == StsAlarm.VALUE)           // ToDo: Need to add user supplied name
            alarm = new StsValueAlarm("Value_Alarm", soundFile, value, threshold, inside);

        // Add alarm to all selected sensors
        if(alarm != null)
        {
            for(int i=0; i<selectedSensors.length; i++)
                ((StsSensor)selectedSensors[i]).addAlarm(alarm);
        }
        return super.end();
    }

    public void previous()
    {
        if(currentStep != selectSensors)
            gotoStep(selectSensors);
    }

    public void next()
    {
        if(alarmType == StsAlarm.SURFACE)
        {
            if(hasSurfaces())
                gotoStep(surfaceDefine);
            else
                new StsMessage(frame, StsMessage.ERROR, "No surfaces found in project, unable to create surface alarm.");
        }
        else if(alarmType == StsAlarm.WELL)
        {
            if(hasWells())
                gotoStep(wellDefine);
            else
                new StsMessage(frame, StsMessage.ERROR, "No wells found in project, unable to create well alarm.");
        }
        else if(alarmType == StsAlarm.VALUE)
            gotoStep(valueDefine);
    }

    public Object getSelectedSensors() { return selectedSensors; }
    public void setSelectedSensors(Object sensor)
    {
        if(sensor == null) return;
        selectedSensors = selectSensors.getSelectedSensors();
    }

    public Object getSelectedSurface() { return selectedSurface; }
    public void setSelectedSurface(Object surface)
    {
        if(surface == null) return;
        selectedSurface = surfaceDefine.getSelectedSurface();
    }

    public Object getSelectedWell() { return selectedWell; }
    public void setSelectedWell(Object well)
    {
        if(well == null) return;
        selectedWell = wellDefine.getSelectedWell();
    }
    public String getDirection() { return StsSurfaceAlarm.REF_DIRECTION[refDirection]; }
    public void setDirection(String type)
    {
        if(type.equalsIgnoreCase(StsSurfaceAlarm.REF_DIRECTION[StsSurfaceAlarm.BELOW]))
            refDirection = StsSurfaceAlarm.BELOW;
        else
            refDirection = StsSurfaceAlarm.ABOVE;
    }

    public String getAlarmType() { return StsAlarm.ALARM_TYPE_STRINGS[alarmType]; }
    public void setAlarmType(String type) { alarmType = StsAlarm.getAlarmType(type); }

    public String getQueueType()  { return StsAlarm.QUEUE_TYPE_STRINGS[queueType]; }
    public void setQueueType(String type) { queueType = StsAlarm.getQueueType(type); }

    public float getOffset()  { return offset; }
    public void setOffset(float offset) { this.offset = offset; }
    public float getValue()  { return value; }
    public void setValue(float val) { this.value = val; }
    public float getThreshold()  { return threshold; }
    public void setThreshold(float val) { this.threshold = val; }
    public boolean getInside()  { return inside; }
    public void setInside(boolean inside) { this.inside = inside; }
    public String getSoundFile()  { return soundFile; }
    public void setSoundFile(String sound)
    {
        this.soundFile = sound;
        StsSound.play(sound);
    }

	public boolean checkPrerequisites()
	{
		StsDynamicSensorClass sensorClass = (StsDynamicSensorClass) model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
		StsObject[] sensors = sensorClass.getSensors();
		if (sensors.length == 0) {
			reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
			return false;
		}
        return true;
    }

    public boolean hasSurfaces()
    {
        int nSurfaces = ((StsSurface[])model.getCastObjectList(StsSurface.class)).length;
		if(nSurfaces < 1)
        {
			return false;
		}
        return true;
	}
    public boolean hasWells()
    {
        int nWells = ((StsWell[])model.getCastObjectList(StsWell.class)).length;
		if(nWells < 1)
        {
			return false;
		}
        return true;
	}

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsAlarmsWizard vvWizard = new StsAlarmsWizard(actionManager);
        vvWizard.start();
    }
}
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System
package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPreStackWizard extends StsSeismicWizard implements Runnable
{
    protected StsPreStackMultiVolFormat multiVolumeFormat;
    protected StsPostStackTraceDefinition traceDef;
    protected StsPreStackSurveyDefinition surveyDef;
    protected StsPreStackBatch filesToProcess;
    protected StsPreStackProcess volumeProcess;

    private boolean ignoreMultiVolume = true;

    public StsPreStackWizard(StsActionManager actionManager)
    {
        this(actionManager, 675, 768);
    }

    public StsPreStackWizard(StsActionManager actionManager, int width, int height)
    {
        super(actionManager, width, height, StsSEGYFormat.PRESTACK_RAW);
    }

    protected void addSteps()
    {
        addSteps
        (
            new StsWizardStep[]
            {
                    volumeSelect      = new StsPreStackSelect(this),
                    fileFormat        = new StsPreStackFileFormat3d(this),
                    multiVolumeFormat = new StsPreStackMultiVolFormat(this),
                    traceDef          = new StsPreStackTraceDefinition3d(this),
                    surveyDef         = new StsPreStackSurveyDefinition(this),
                    filesToProcess    = new StsPreStackBatch(this),
                    volumeProcess     = new StsPreStackProcess(this)
            }
        );
    }

    public boolean start()
    {
        return start("Process PreStack Volumes");
    }

    public byte getSegyFormatDataType() { return StsSEGYFormat.PRESTACK_RAW; }
    
    protected StsSeismicBoundingBox createNewStsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String outputDirectory, StsSEGYFormat format)
    {
        return StsPreStackSegyLine3d.constructor(seismicWizard, file, outputDirectory, frame, format);
    }

    public StsPostStackAnalyzer getAnalyzer(StsProgressPanel panel, StsTablePanelNew volumeStatusPanel)
    {
        return StsPreStack3dAnalyzer.getAnalyzer(this, panel, volumeStatusPanel);
    }

    public void previous()
    {
        if(segyVolumesList.size() == 0)
            gotoStep( volumeSelect);
        
        if(currentStep == filesToProcess)
            gotoStep(volumeSelect);
        else if(currentStep == traceDef)
        {
            if( segyVolumesList.size() > 0)
            {
                StsPreStackSegyLine segyVolume = (StsPreStackSegyLine)segyVolumesList.get(0);
                if(!segyFormat.isMultiVolume(getIsLittleEndian(), segyVolume.getBinaryHeader()))
                    gotoStep(fileFormat);
                else
                    gotoPreviousStep();
            }
            else
                gotoStep( volumeSelect);            
        }
        else
        {
            gotoPreviousStep();
        }
    }

    public void next()
    {
        if(currentStep == volumeSelect)
        {
            if(segyVolumesList.size() ==  0)
            {
                if(segyVolumesToProcess != null && segyVolumesToProcess.length > 0)
                {
                    gotoStep(filesToProcess);
                    return;
                }
                new StsMessage(frame, StsMessage.ERROR, "No volumes selected or name specified: select or cancel.");
                return;
            }
            else if(volumeSelect.panel.getVolumeName() == null)
            {
                new StsMessage(frame, StsMessage.ERROR, "No name specified: enter name or cancel.");
                return;
            }
            gotoNextStep();
            // reanalyzeTraces(true, true, getProgressPanel(), true, false, true);
            return;
        }

        if(currentStep == fileFormat)
        {
            if(segyVolumesList.size() >  0)
            {
                StsPreStackSegyLine segyVolume = (StsPreStackSegyLine)segyVolumesList.get(0);
                if(segyVolume.segyData.isMultiVolumeFile)
                {
                    ignoreMultiVolume = false;
                    gotoNextStep();
                }
                else
                {
                    gotoStep(traceDef);
                    // reanalyzeTraces(false, false, getProgressPanel(), true, true, false);
                }
                return;
            }
        }
        if(currentStep == multiVolumeFormat)
        {
            gotoNextStep();
            return;
        }
        if(currentStep == traceDef)
        {
            enableFinish();
            gotoNextStep();

            //reanalyzeTraces(false, false, null, true, true, false);

            return;
        }
        if(currentStep == surveyDef)
        {
            saveSEGYFormats();
            filesToProcess.addVolumesForProcessing();
            gotoNextStep();
            return;
        }

        if(currentStep == filesToProcess)
        {
            volumeProcess.constructPanel();
            gotoNextStep();
            return;
        }

        gotoNextStep();
    }

    public boolean isStatusOk( StsSeismicBoundingBox segyVolume)
    {
        if( currentStep == fileFormat)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_FILE_OK);
        if( currentStep == traceDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_TRACES_OK);
        if( currentStep == surveyDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_GEOMETRY_OK);
        return false;
    }

    public boolean getNMOed()
    {
        StsSeismicBoundingBox[] volumes = getSelectedVolumes();
        if(volumes.length == 0) return false;
        return ((StsPreStackSegyLine)volumes[0]).isNMOed;
    }

    public void setNMOed(boolean value)
    {
        StsSeismicBoundingBox[] volumes = getSelectedVolumes();
        if(volumes == null) return;
        for(int i = 0; i < volumes.length; i++)
        {
            ((StsPreStackSegyLine)volumes[i]).setIsNMOed(value);
        }
    }


    public float getVelocityInc()
    {
        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return 0.0f;
        return (float)segyFormat.getBinaryHdrValue("VELINC", segyVolume.getBinaryHeader(), getIsLittleEndian()) - 1.0f;
    }

    public int getNumCdps()
    {
        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return 0;
        return (int)segyFormat.getBinaryHdrValue("NUMCDPS", segyVolume.getBinaryHeader(), getIsLittleEndian());
    }

    public int getNumVolumes()
    {
        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return 0;
        return (int)segyFormat.getBinaryHdrValue("NUMVOLS", segyVolume.getBinaryHeader(), getIsLittleEndian());
    }

    public String getMultiVolType()
    {
        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return "";
        return StsSEGYFormat.multiVolumeTypes[(int)segyFormat.getBinaryHdrValue("MVOLTYPE", segyVolume.getBinaryHeader(), getIsLittleEndian())];
    }

    public void setMultiVolType(String type)
    {
    }

    public String getTimeValuePair(int num)
    {
        String timeName = "TIME" + num;
        String valueName = "VALUE" + num;
        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return "";
        float time = (float)segyFormat.getBinaryHdrValue(timeName, segyVolume.getBinaryHeader(), getIsLittleEndian());
        float value = (float)segyFormat.getBinaryHdrValue(valueName, segyVolume.getBinaryHeader(), getIsLittleEndian());
        return "(" + time + "," + value + ")";
    }

    public String[] getTimeValuePairs()
    {
        String[] tvPairs = new String[10];
        for(int i = 0; i < 10; i++)
        {
            tvPairs[i] = getTimeValuePair(i + 1);
        }
        return tvPairs;
    }

    public void setTimeValuePairs(String[] pairs) {}

    public boolean getIgnoreMultiVolume() { return ignoreMultiVolume; }
    public void setIgnoreMultiVolume(boolean value) { ignoreMultiVolume = value; }

    public static void main(String[] args)
    {
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsModel model = new StsModel();
                model.setProject(new StsProject());
                StsActionManager actionManager = new StsActionManager(model);
                StsPreStackWizard volumeWizard = new StsPreStackWizard(actionManager);
                volumeWizard.start();
            }
        };
        StsToolkit.runWaitOnEventThread(winThread);
    }
}

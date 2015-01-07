package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSegyVspWizard extends StsPostStackWizard implements Runnable
{
	public StsSegyVspWizard(StsActionManager actionManager)
	{
        super(actionManager, 675, 768, StsSEGYFormat.VSP);
    }

    protected StsSeismicBoundingBox createNewStsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String outputDirectory, StsSEGYFormat format)
    {
        return StsSegyVsp.constructor(seismicWizard, file, outputDirectory, frame, format);
    }

    public StsPostStackAnalyzer getAnalyzer(StsProgressPanel panel, StsTablePanelNew volumeStatusPanel)
    {
         return StsPostStackVspAnalyzer.getAnalyzer(this, panel, volumeStatusPanel);
    }

    protected void addSteps()
    {
        addSteps
        (
            new StsWizardStep[]
            {
                volumeSelect   = new StsPostStackSelect(this),
                fileFormat     = new StsSegyVspFileFormat(this),
                traceDef       = new StsSegyVspTraceDefinition(this),
                filesToProcess = new StsSegyVspBatch(this),
                volumeProcess  = new StsSegyVspProcess(this)
            }
        );
    }

    public boolean start()
    {
        return start("Process VSP data");
    }

    public byte getSegyFormatDataType() { return StsSEGYFormat.VSP; }

    public void previous()
    {
        if(segyVolumesList.size() == 0)
            gotoStep( volumeSelect);
        
        if(currentStep == filesToProcess)
        {
            gotoStep(volumeSelect);
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
            }
            else if(volumeSelect.panel.getVolumeName() == null)
            {
                new StsMessage(frame, StsMessage.ERROR, "No name specified: enter name or cancel.");
            }
            else
            {
                gotoNextStep();
            }
//            reanalyzeTraces(true, true, getProgressPanel(), true, false, true);
        }
        else if(currentStep == traceDef)
        {
            saveSEGYFormats();
            filesToProcess.addVolumesForProcessing();
            gotoNextStep();
        }
        else if(currentStep == filesToProcess)
        {
            volumeProcess.constructPanel();
            gotoNextStep();
        }
        else
            gotoNextStep();
    }
/*
    public void finish()
    {
        //next();
        super.finish();
        new StsMessage(getModel().win3d, StsMessage.INFO, "You will need to run the Load Poststack PostStack3d Step to visualize the data.");
    }
*/
    public boolean isStatusOk( StsSeismicBoundingBox segyVolume)
    {
        if( currentStep == fileFormat)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_FILE_OK);
        if( currentStep == traceDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_TRACES_OK);
        return false;
    }

    // added for StsPostStackWizard API consistency
    public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
    {
        return null;
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


    public static void main(String[] args)
    {
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsModel model = StsModel.constructor();
                StsActionManager actionManager = new StsActionManager(model);
                StsSegyVspWizard volumeWizard = new StsSegyVspWizard(actionManager);
                volumeWizard.start();
            }
        };
        StsToolkit.runWaitOnEventThread(winThread);
    }
}

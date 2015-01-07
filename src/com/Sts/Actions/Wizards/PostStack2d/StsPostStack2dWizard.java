package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

public class StsPostStack2dWizard extends StsPostStackWizard 
{
    StsPostStackSurveyDefinition2d surveyDef;

    public StsPostStack2dWizard(StsActionManager actionManager)
    {
      super(actionManager, 600, 768, StsSEGYFormat.POSTSTACK2D);
    }

    public StsPostStackAnalyzer getAnalyzer(StsProgressPanel panel, StsTablePanelNew volumeStatusPanel)
     {
         return StsPostStack2dAnalyzer.getAnalyzer(this, panel, volumeStatusPanel);
     } 

    public boolean start()
    {
        return start("Process PostStack 2d Volumes");
    }

    public byte getSegyFormatDataType() { return StsSEGYFormat.POSTSTACK2D; }

    protected void addSteps()
    {
       addSteps(new StsWizardStep[]
        {
            volumeSelect   = new StsPostStackSelect(this),
            fileFormat     = new StsPostStackFileFormat2d(this),
            traceDef       = new StsPostStackTraceDefinition2d(this),
            surveyDef      = new StsPostStackSurveyDefinition2d(this),
//            rangeEdit      = new StsPostStackRangeEdit(this),
            filesToProcess = new StsPostStackBatch(this),
            volumeProcess  = new StsPostStackProcess2d(this)
        });
    }
    public void previous()
    {
        if(segyVolumesList.size() == 0)
            gotoStep( volumeSelect);

        if(currentStep == filesToProcess)
            gotoStep(volumeSelect);
        else if(currentStep == rangeEdit)
        {
            if(getOverrideGeometry())
                gotoPreviousStep();
            else
                gotoStep(traceDef);
        }
        else
            gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == volumeSelect)
        {
            if(segyVolumesList.size() == 0)
            {
                if(segyVolumesToProcess != null && segyVolumesToProcess.length > 0)
                {
                    gotoStep(filesToProcess);
                    return;
                }
                new StsMessage(frame, StsMessage.ERROR, "No volumes selected or name specified: select or cancel.");
            }
            else if(volumeSelect.panel.getVolumeName() == null)
                new StsMessage(frame, StsMessage.ERROR, "No name specified: enter name or cancel.");
            else
            {
//                setVolumesSegyFormat();
            	setSelectedVolumesToAll();
                gotoNextStep();
            }
        }

        else if(currentStep == fileFormat)
        {
            gotoStep(traceDef);
        }
        else if(currentStep == this.traceDef)
        {
            gotoNextStep();
        }
        // TODO construct rangeEdit2dPanel and step
    /*
        else if(currentStep == surveyDef)
            gotoNextStep();
        else if(currentStep == rangeEdit)
        {
            saveSEGYFormats();
            filesToProcess.addVolumesForProcessing();
            gotoNextStep();
        }
    */
        else if(currentStep == surveyDef)
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

    protected boolean hasOutputFiles( StsSeismicBoundingBox volume)
    {
        return ((StsSegyLine2d)volume).hasOutputFiles();
    }

    public boolean isStatusOk( StsSeismicBoundingBox segyVolume)
    {
        if( currentStep == fileFormat)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_FILE_OK);
        if( currentStep == traceDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_TRACES_OK);
        if( currentStep == surveyDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_GEOMETRY_OK);
        if( currentStep == rangeEdit)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_GEOMETRY_OK);
        return false;
    }

    public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
    {
        return surveyDef.panel.getSurveyDefinitionBoundingBox();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }
    
    protected StsSeismicBoundingBox createNewStsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String outputDirectory, StsSEGYFormat format)
    {
        return StsSegyLine2d.constructor(seismicWizard, file, outputDirectory, frame, format);
    }

    public boolean cropEnabled() { return false;}
}

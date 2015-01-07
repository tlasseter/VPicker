package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import javax.swing.*;

public class StsPostStack3dWizard extends StsPostStackWizard implements Runnable
{
    protected StsPostStackSurveyDefinition3d surveyDef;
 
    public StsPostStack3dWizard(StsActionManager actionManager)
    {
        super(actionManager, 750, 900, StsSEGYFormat.POSTSTACK);
    }

    public StsPostStackAnalyzer getAnalyzer(StsProgressPanel panel, StsTablePanelNew volumeStatusPanel)
    {
        return StsPostStack3dAnalyzer.getAnalyzer(this, panel, volumeStatusPanel);
    }

    protected void addSteps()
    {
        addSteps
        (
            new StsWizardStep[]
            {
                volumeSelect   = new StsPostStackSelect(this),
                fileFormat     = new StsPostStackFileFormat3d(this),
                traceDef       = new StsPostStackTraceDefinition3d(this),
                surveyDef      = new StsPostStackSurveyDefinition3d(this),
                rangeEdit      = new StsPostStackRangeEdit(this),
                filesToProcess = new StsPostStackBatch(this),
                volumeProcess  = new StsPostStackProcess3d(this)
            }
        );
    }

    public boolean start()
    {
        return start("Process PostStack 3d Volumes");
    }

    public byte getSegyFormatDataType() { return StsSEGYFormat.POSTSTACK; }

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
            gotoStep(traceDef);
        else if(currentStep == this.traceDef)
        {
//            enableFinish();
            if(getOverrideGeometry())
                gotoNextStep();
            else
                gotoStep(rangeEdit);
        }
        else if(currentStep == surveyDef)
            gotoNextStep();
        else if(currentStep == rangeEdit)
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

    public void cancel()
    {
        if(currentStep == volumeProcess)
            volumeProcess.panel.cancel();
        else if(segyVolumesToProcess != null && segyVolumesToProcess.length > 0)
        {
            int answer = JOptionPane.showConfirmDialog(dialog, "There are volumes awaiting processing. Are you sure you want to cancel?", "Cancel Wizard", JOptionPane.YES_NO_OPTION);
            if(answer == JOptionPane.NO_OPTION) return;
        }
        super.cancel();
    }

    public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
    {
        return surveyDef.panel.getSurveyDefinitionBoundingBox();
    }
}
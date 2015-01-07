package com.Sts.Actions.Wizards.PreStack2d;

import com.Sts.Actions.Wizards.PreStack3d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;

public class StsPreStack2dWizard extends StsPreStackWizard
{
   public StsPreStack2dWizard(StsActionManager actionManager)
   {
	  super(actionManager, 600, 750);
   }

   protected void addSteps()
   {
	   addSteps(new StsWizardStep[]
				{
				volumeSelect   = new StsPreStackSelect(this),
				fileFormat     = new StsPreStackFileFormat2d(this),
				traceDef       = new StsPreStackTraceDefinition2d(this),
				filesToProcess = new StsPreStackBatch(this),
				volumeProcess  = new StsPreStackProcess(this)
	   }
		   );
   }
   
   public boolean end()
   {
       model.setActionStatus(getClass().getName(), StsModel.STARTED);
       return super.end();
   }
   
   public void next()
   {
	  if (currentStep == traceDef)
	  {
		  saveSEGYFormats();
		  filesToProcess.addVolumesForProcessing();
          enableFinish();
          gotoNextStep();
          return;
	  }
      else
        super.next();
    }

    public boolean isStatusOk( StsSeismicBoundingBox segyVolume)
    {
        if( currentStep == fileFormat)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_FILE_OK);
        if( currentStep == traceDef)
            return (segyVolume.status >= StsSeismicBoundingBox.STATUS_TRACES_OK);
        return false;
    }

    protected StsSeismicBoundingBox createNewStsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String outputDirectory, StsSEGYFormat format)
    {
        return StsPreStackSegyLine2d.constructor(seismicWizard, file, outputDirectory, frame, format);
    }

   public boolean getOverrideGeometry()
   {
	  return false;
   }
}

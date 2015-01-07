package com.Sts.Actions.Wizards.MicroseismicPreStack;

import com.Sts.Actions.Wizards.PreStack2d.*;
import com.Sts.Actions.Wizards.PreStack3d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsMicroPreStackWizard extends StsPreStackWizard
{
   long timeLong = System.currentTimeMillis();
   long gatherDuration = 60000;
   int byteLocation = 0;
   int xByteLocation = 81;
   int yByteLocation = 85;
   boolean useHeaders = false;
   boolean hdrCoordinates = false;
   
   StsAbstractFileSet fileSet;
   TreeMap filenameTreeMap;
   
   StsDefineTimeLimits defineLimits;
   StsMicroPreStackLoad microSetProcess;
   StsOrderFiles orderFiles;
   
   public StsMicroPreStackWizard(StsActionManager actionManager)
   {
	  super(actionManager, 600, 750);	  
   }

   protected void addSteps()
   {
	   addSteps(new StsWizardStep[]	{ volumeSelect   = new StsPreStackSelect(this),
									  fileFormat     = new StsPreStackFileFormat2d(this),			   
			   						  defineLimits = new StsDefineTimeLimits(this),
			   						  orderFiles = new StsOrderFiles(this),
			   						  microSetProcess = new StsMicroPreStackLoad(this) } );
   }
   
   public boolean start()
   { 
	   if(!super.start()) return false;	   
	   return true;
   }

   public void setStartTime(String timeStg)
   {
	   timeLong = StsDateFieldBean.convertToLong(timeStg);
	   defineLimits.panel.setCheckFields();
   }
   public String getStartTime()
   {
       String timeStg = StsDateFieldBean.convertToString(timeLong);
       return timeStg;
   }
   public long getStartTimeLong()
   {
	   return timeLong;
   }
   public void setGatherDuration(long duration)
   {
	   gatherDuration = duration;
	   defineLimits.panel.setCheckFields();
   }
   public long getGatherDuration() { return gatherDuration; }
   public int getXByteLocation() { return xByteLocation; }
   public void setXByteLocation(int location)
   {
	   xByteLocation = location;
	   getSelectedSegyDatasets()[0].getSegyFormat().getTraceRec(StsSEGYFormat.SHT_X).setLoc(xByteLocation);  
	   defineLimits.panel.setCheckFields();	   
   }
   
   public int getYByteLocation() { return yByteLocation; }
   public void setYByteLocation(int location)
   {
	   yByteLocation = location;
	   getSelectedSegyDatasets()[0].getSegyFormat().getTraceRec(StsSEGYFormat.SHT_Y).setLoc(yByteLocation);
	   defineLimits.panel.setCheckFields();
   }
   public int getByteLocation() { return byteLocation; }
   
   public boolean getUseHeaders() { return useHeaders; }
   public void setUseHeaders(boolean use)
   {
	   useHeaders = use;
	   defineLimits.panel.setUseHeaders(use);
	   defineLimits.panel.setCheckFields();
   }
   public boolean getCoordinates() { return hdrCoordinates; }
   public void setCoordinates(boolean hdr)
   {
	   hdrCoordinates = hdr;
	   defineLimits.panel.setCoordinates(hdr);
	   defineLimits.panel.setCheckFields();
   }   
   
   public boolean end()
   {
       model.setActionStatus(getClass().getName(), StsModel.STARTED);
       return super.end();
   }
   
   public void next()
   {
       if(currentStep == volumeSelect)
       {
    	   if(volumeSelect.panel.getVolumeName() == null)
           {
               new StsMessage(frame, StsMessage.ERROR, "No name specified: enter name or cancel.");
               return;
           }
           gotoNextStep();
           return;
       }

       if(currentStep == fileFormat)
       {
           if(segyVolumesList.size() >  0)
               gotoStep(defineLimits);
           else
        	   new StsMessage(frame, StsMessage.WARNING, "Select at least one file first.");
           return;
       }
       if(currentStep == defineLimits)
       {
    	   if(getUseHeaders())
    	   {
    		   microSetProcess.constructPanel();
    		   enableFinish();
    		   gotoNextStep();
    		   return;
    	   }
       }
       if(currentStep == orderFiles)
       {
    	   microSetProcess.constructPanel();
    	   enableFinish();
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
   
   public long computeStartTime()
   {
	   
	   return 0l;
   }
   
   public boolean createMicroseismicSet(StsProgressPanel panel)
   {
	   float zMin = StsParameters.largeFloat;
	   float zMax = -StsParameters.largeFloat;
	   float dataMin = StsParameters.largeFloat;
	   float dataMax = -StsParameters.largeFloat;
	   
	   StsPreStackMicroseismicSet microSet = new StsPreStackMicroseismicSet(false);
	   microSet.setName("microSet-0");
	   microSet.setMicroSegyFormat(segyFormat);
	   StsSeismicBoundingBox[] boxes = this.getSegyVolumes();	   
	   panel.initialize(boxes.length);
       for (int n = 0; n < boxes.length; n++)
       {
    	   long start = StsDateFieldBean.convertToLong(boxes[n].getBornDate());
    	   long end = start + getGatherDuration();
    	   microSet.addMicroseismicGather(boxes[n].getSegyDirectory() + boxes[n].getSegyFilename(), start, end);
       	
    	   if(zMin > boxes[n].getZMin())
    		   zMin = boxes[n].getZMin();
    	   if(zMax < boxes[n].getZMax())
    		   zMax = boxes[n].getZMax();
    	   if(dataMin > boxes[n].getDataMin())
    		   dataMin = boxes[n].getDataMin();
    	   if(dataMax < boxes[n].getDataMax())
    		   dataMax = boxes[n].getDataMax();
    	   
    	   panel.setValue(n+1);
    	   panel.setDescription("Loaded file : " + (n+1) + " of " + boxes.length);
    	   panel.appendLine("Completed loading file: " + boxes[n].getSegyFilename());
       }

       if (!microSet.addToProject())
       {
    	   new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add microseismic set to project.");
    	   new StsException(StsException.WARNING, "Failed to add prestack microseismic to project.");
       }
       microSet.initialize(zMin,zMax,dataMin,dataMax);       
	   return true;
   }
}

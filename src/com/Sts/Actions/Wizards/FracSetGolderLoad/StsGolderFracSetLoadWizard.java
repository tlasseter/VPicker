
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetGolderLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsActionManager;
import com.Sts.MVC.StsModel;
import com.Sts.Utilities.StsMath;

public class StsGolderFracSetLoadWizard extends StsWizard // implements Runnable
{
    public StsGolderFracSetSelect fracSetSelect = null;
    public StsGolderFracSetLoad fracSetLoad = null;
    private StsFile[] fracSetFiles = null;

    private StsWizardStep[] mySteps =
    {
        fracSetSelect = new StsGolderFracSetSelect(this),
        fracSetLoad = new StsGolderFracSetLoad(this)
    };

    public StsGolderFracSetLoadWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 600); 
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("FracSet Load Wizard");
        initialize();
        disableFinish();
        if(!super.start()) return false;
        return true;
    }

    public void initialize()
    {
//        hUnits = model.getProject().getXyUnits();
//        vUnits = model.getProject().getDepthUnits();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);    	
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
//        if(currentStep == FracSetSelect)
//        {
//            if(defineTimeRequired() == true)
//                gotoStep(FracSetTime);
//            else if(defineStaticLocationRequired())
//                gotoStep(FracSetStatic);
//            else if(defineRelativeLocationRequired())
//            	gotoStep(FracSetRelative);
//            else
//               gotoLoad();
//        }
//        else if(currentStep == FracSetTime)
//        {
//            if(defineStaticLocationRequired())
//               gotoStep(FracSetStatic);
//            else if(defineRelativeLocationRequired())
//            	gotoStep(FracSetRelative);            
//            else
//               gotoLoad();
//        }
//        else if(currentStep == FracSetStatic)
//        {
//        	if(defineRelativeLocationRequired())
//            	gotoStep(FracSetRelative); 
//        	else
//        		gotoLoad();
//        }
//        else
            gotoLoad();

    }

    public void gotoLoad()
    {
 //       prepareFileSets();
        fracSetLoad.constructPanel();
        gotoStep(fracSetLoad);
    }
    
//    public byte getFracSetType() { return FracSetSelect.panel.getFracSetType(); }
//    private boolean defineTimeRequired()
//    {
//        StsFracSetFile[] selectedFiles = getFracSetFiles();
//        for(int i=0; i<selectedFiles.length; i++)
//        {
//            if((selectedFiles[i].timeType != StsFracSetKeywordIO.TIME_AND_DATE) &&
//            	(selectedFiles[i].timeType != StsFracSetKeywordIO.TIME_OR_DATE))
//            {
//               return true;
//            }
//        }
//        return false;
//    }
//    
//    private boolean defineStaticLocationRequired()
//    {
//        StsFracSetFile[] selectedFiles = getFracSetFiles();
//        for(int i=0; i<selectedFiles.length; i++)
//        {
//            if(selectedFiles[i].positionType == StsFracSet.STATIC)
//            {
//               return true;
//            }
//        }
//        return false;
//    }
//    
//    private boolean defineRelativeLocationRequired()
//    {
//        StsFracSetFile[] selectedFiles = getFracSetFiles();
//        for(int i=0; i<selectedFiles.length; i++)
//        {
//            if(selectedFiles[i].relativeCoordinates)
//            {
//               return true;
//            }
//        }
//        return false;
//    }    
    
//    public void prepareFileSets()
//    {
//       // Loop through each FracSet file and create its own FracSet set.
//       StsGolderFracSetFile[] selectedFiles = getFracSetFiles();
//       for(int i=0; i<selectedFiles.length; i++)
//       {
//           //StsFracSetKeywordIO.parseAsciiFilename(FracSetFiles[i].file.getFilename());
//           String[] fileNames = new String[] { selectedFiles[i].file.getFilename() };
//           //StsFracSetImport.addFracSetFilenameSets(fileNames, StsFracSetImport.ASCIIFILES);
//       }
//    }

    public void finish()
    {
        super.finish();
    }

//    public double[] getPosition()
//    {
//        return null;
//    }
//    public double[] getRelativePosition()
//    {
//        return null;
//    }  
    
    public String getFracSetName()
    {
        return "FracSetName";
    }
    
    public boolean getArchiveIt()
    {
        return false;
    }

//    public String getAsciiStartTime()
//    {
//        Date date = new Date(0L);
//        SimpleDateFormat format = new SimpleDateFormat(model.project.getTimeDateFormatString());
//        String time = format.format(date);
//        try
//        {
////            time = FracSetTime.panel.getStartTime();
//            if(time == null)
//            {
//                new StsMessage(frame, StsMessage.ERROR,
//                    "Invalid date & time value (" + format.format(date) + ") in time series file selection step.\n" +
//                    "\n   Solution: Return to step and re-enter valid start time.");
//                return null;
//            }
//            date = format.parse(time);
//        }
//        catch (Exception e)
//        {
//            StsMessageFiles.logMessage("Failed to create date, setting to 01-01-71 00:00:00.0");
//            return "01-01-71 00:00:00.0";
//        }
//        return time;
//    }

    public void addFile(StsFile file)
    {
        fracSetFiles = (StsFile[]) StsMath.arrayAddElement(fracSetFiles, file);
    }

    public void removeFile(StsFile file)
    {
    	fracSetFiles = (StsFile[])StsMath.arrayDeleteElement(fracSetFiles, file);
    }

    public StsFile[] getFracSetFiles() { return fracSetFiles; }

    static String cleanTimeString(String ts)
    {
        if(ts.length() == 19)
            return ts;
        if(ts.indexOf(".") == -1)
            ts = ts + ".0";
        if(ts.indexOf(".") == 7)
            ts = "0" + ts;
        if((ts.indexOf(".") < 7) || (ts.indexOf(".") > 8))
            return null;
        ts = "2000-01-01 " + ts;
        return ts;
    }

}

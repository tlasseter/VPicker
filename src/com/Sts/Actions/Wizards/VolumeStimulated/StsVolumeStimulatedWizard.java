//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsVolumeStimulatedWizard extends StsWizard
 {
     public StsFractureSet[] selectedFractureSets = null;
     public Object[] selectedSubVolumeObjects = null;
     StsProgressPanel progressPanel;

     public StsDefineVolumetrics defineVolumetrics = null;
     public StsSelectFractureSets selectFractureSets = null;
     public StsSelectSubVolumes selectSubVolumes = null;
     public StsSelectType stimulatedType = null;
     public StsSelectVolume selectVolume = null;
     final static byte VOLUME = 0;
     final static byte FRACTURES = 1;

     final static int CUFT = 0;
     final static int CUMTR = 1;
     final static int CUYD = 2;
     final static int ACREFT = 3;
     final static String[] unitStrings = {"Cubic Feet", "Cubic Meters", "Cubic Yards", "Acre Feet"};
     int units = 0;

     int numFractureSets = 0;

     float areaScale = 1.0f;
     float effectivePerm = 0.05f;
     float effectivePorosity = 2.0f;
     float slabThickness = 1.0f;

     public StsVolumeStimulatedWizard(StsActionManager actionManager)
     {
         super(actionManager, 300, 500);
         addSteps
             (
                 new StsWizardStep[]
                     {
                         stimulatedType = new StsSelectType(this),
                         selectSubVolumes = new StsSelectSubVolumes(this),
                         selectFractureSets = new StsSelectFractureSets(this),
                         selectVolume = new StsSelectVolume(this)
                         //defineVolumetrics = new StsDefineVolumetrics(this)
                     }
             );
     }

     public boolean start()
     {
         disableFinish();
         System.runFinalization();
         System.gc();
         dialog.setTitle("Volume Stimulated");
         units = CUFT;
         if (model.getProject().getXyUnits() == StsParameters.DIST_METER)
             units = CUMTR;
         return super.start();
     }

     public void previous()
     {
         if(currentStep == selectVolume)
            gotoStep(selectSubVolumes);
         gotoPreviousStep();
     }

     public void next()
     {
         /*
         if(currentStep == selectFractureSets)
         {
             if(StsMessage.questionValue(frame, "Do you want to simulate the stimulated volume?"))
                 gotoNextStep();
             else
                 finish();
         }
         else
         */
         if (currentStep == selectSubVolumes)
         {
             // Turn off all subvolumes
             StsSubVolumeClass subVolumeClass = (StsSubVolumeClass) model.getStsClass(StsSubVolume.class);
             if (subVolumeClass != null && subVolumeClass.getIsVisible())
             {
                 StsSubVolume[] availableSubVolumes = subVolumeClass.getSubVolumes();
                 for (int i = 0; i < availableSubVolumes.length; i++)
                     availableSubVolumes[i].setIsApplied(false);
             }
             // Activate the selected ones
             selectedSubVolumeObjects = selectSubVolumes.panel.getSelectedObjects();
             if (selectedSubVolumeObjects != null)
                 for (int i = 0; i < selectedSubVolumeObjects.length; i++)
                     ((StsSubVolume) selectedSubVolumeObjects[i]).setIsApplied(true);
             if(stimulatedType.getType() == VOLUME)
                 gotoStep(selectVolume);
             else
                 gotoStep(selectFractureSets);
         }
         else
         gotoNextStep();
     }

     public void setFractureSets(StsFractureSet[] sets)
     {
         selectedFractureSets = sets;
         setStimulated(true);
     }

     public StsFractureSet[] getFractureSets() { return selectedFractureSets; }

     public void setSubVolumeObjects(Object sv)
     {
         if (selectSubVolumes == null)
             selectedSubVolumeObjects = null;
         else
             selectedSubVolumeObjects = selectSubVolumes.panel.getSelectedObjects();
     }

    /* Does not appear to work, one above does.
    public void setSubVolumeObjects(Object[] subVolumeObjects)
    {
        selectedSubVolumeObjects = subVolumeObjects;
        if(selectedSubVolumeObjects == null) return;
        for(int n = 0; n < selectedSubVolumeObjects.length; n++)
        {
            StsSubVolume subVolume = (StsSubVolume)selectedSubVolumeObjects[n];
            subVolume.setIsApplied(true);
        }
    } 
    */   
    public Object getSubVolumeObjects() { return selectedSubVolumeObjects; }

     public StsSubVolume[] getSubVolumes()
     {
         return (StsSubVolume[]) StsMath.arraycopy(selectedSubVolumeObjects, StsSubVolume.class);
     }

     public void finish()
     {
         setStimulated(false);
         super.finish();
     }

     public void cancel()
     {
         setStimulated(false);
         super.cancel();
     }

     public int getUnitsType()
     {
         return units;
     }

     public String getUnits() { return unitStrings[units]; }

     public void setUnits(String unitStg)
     {
         for (int i = 0; i < unitStrings.length; i++)
             if (unitStg.equalsIgnoreCase(unitStrings[i]))
                 units = i;
     }

     private void setStimulated(boolean draw)
     {
         if (selectedFractureSets != null)
         {
             for (int i = 0; i < selectedFractureSets.length; i++)
                 selectedFractureSets[i].setDrawStimulated(draw);
         }
     }

     public float getSlabThickness() { return slabThickness; }

     public void setSlabThickness(float thickness)
     {
         slabThickness = thickness;
     }

     public float getPorosity() { return effectivePorosity; }

     public void setPorosity(float por)
     {
         effectivePorosity = por;
     }

     public float getPermeability() { return effectivePerm; }

     public void setPermeability(float perm)
     {
         effectivePerm = perm;
     }

     public float getAreaScale() { return areaScale; }

     public void setAreaScale(float scale)
     {
         areaScale = scale;
     }

     public boolean checkPrerequisites()
     {
         numFractureSets = model.getNObjects(StsFractureSet.class);
         int numSensorVolumes = model.getNObjects(StsSensorVirtualVolume.class);
         if ((numFractureSets <= 0) && (numSensorVolumes <= 0))
         {
             reasonForFailure = "Volume stimulated can be run on Fracture Sets and Sensor Virtual Volumes. Please interpret fractures, construct volumes or load seismic volumes first.";
             return false;
         }
         else
         {
             return true;
         }
     }

     public boolean analyzeFractureSets()
     {
         float zMin = StsParameters.largeFloat;
         float zMax = -StsParameters.largeFloat;

         byte zDomainProject = model.getProject().getZDomain();
         if(zDomainProject == model.getProject().TD_TIME)
         {
             new StsMessage(frame, StsMessage.WARNING, "Cannot compute volume in time, please cancel wizard, switch to Depth and try again.");
             return false;
         }
         Main.logUsage();

         StsFractureSetClass fractureSetClass = (StsFractureSetClass) model.getStsClass(StsFractureSet.class);
         Iterator fractureIter = fractureSetClass.getFractureIterator(true);
         while (fractureIter.hasNext())
         {
             StsFracture fracture = (StsFracture) fractureIter.next();
             float z = fracture.getZMin();
             zMin = Math.min(zMin, z);
             z = fracture.getZMax();
             zMax = Math.max(zMax, z);
         }

         // Loop through each of the cursor slices to compute the total volume.
         float originalZ = model.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
         StsRotatedGridBoundingBox fractureBoundingBox = fractureSetClass.getFractureBoundingBox();
         StsFractureCursorSection fractureCursorSection = (StsFractureCursorSection) model.getCursor3d().getDisplayableSection(StsCursor3d.ZDIR, StsFractureSet.class);
         if (fractureCursorSection == null) return false;

         zMin = fractureCursorSection.cursorBoundingBox.getNearestSliceZCoor(zMin);
         zMax = fractureCursorSection.cursorBoundingBox.getNearestSliceZCoor(zMax);
         float zInc = fractureCursorSection.cursorBoundingBox.getZInc();

         int nFilledCells = 0;
         int nSlices = (int)((zMax - zMin)/zInc) + 1;
         progressPanel.setMaximum(nSlices);
         int cnt = 0;
         for (float z = zMin; z <= zMax; z += zInc)
         {
             nFilledCells += fractureCursorSection.getNumberOfFilledCells(z);
             progressPanel.setDescription("Processing slice " + z + "....");
             progressPanel.setValue(cnt++);
         }
         progressPanel.setValue(nSlices);
         progressPanel.setDescription("Analysis complete");
         model.getCursor3d().getDisplayableSection(StsCursor3d.ZDIR, StsFractureSet.class).setDirCoordinate(originalZ);
         float cellVolume = fractureBoundingBox.getGridCellVolume();
         float totalVolume = cellVolume * nFilledCells;

         // Convert units to selected
         int fromUnits = CUFT;
         if (model.getProject().getXyUnits() == StsParameters.DIST_METER)
             fromUnits = CUMTR;
         totalVolume = convertVolume(totalVolume, fromUnits, units);

         selectFractureSets.setAnalysisMessage(totalVolume + " " + unitStrings[units] + " stimulated.");
         Main.logUsage();

         return true;
     }

     public boolean analyzeVolume()
     {
         float zMin = StsParameters.largeFloat;
         float zMax = -StsParameters.largeFloat;
         byte zDomainProject = model.getProject().getZDomain();
         if(zDomainProject == model.getProject().TD_TIME)
         {
             new StsMessage(frame, StsMessage.WARNING, "Cannot compute volume in time, please cancel wizard, switch to Depth and try again.");
             return false;
         }
         Main.logUsage();
         StsSeismicVolume volume = selectVolume.getVolume();
         if(volume == null)
         {
             new StsMessage(frame, StsMessage.WARNING, "Please select a volume to compute on and try again.");
             return false;
         }
         zMin = volume.getZMin();
         zMax = volume.getZMax();
         float zInc = volume.getZInc();
         float originalZ = model.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
         int nFilledCells = 0;
         int nSlices = (int)((zMax - zMin)/zInc) + 1;
         progressPanel.setMaximum(nSlices);
         nFilledCells = volume.getNumberOfFilledCells(progressPanel);
         progressPanel.setValue(nSlices);
         progressPanel.setDescription("Analysis complete");
         float cellVolume = volume.getGridCellVolume();
         float totalVolume = cellVolume * nFilledCells;
         int fromUnits = CUFT;
         if (model.getProject().getXyUnits() == StsParameters.DIST_METER)
             fromUnits = CUMTR;
         totalVolume = convertVolume(totalVolume, fromUnits, units);
         selectVolume.setAnalysisMessage(totalVolume + " " + unitStrings[units] + " stimulated.");
         Main.logUsage();
         return true;
     }
     public float convertVolume(float volume, int fromUnits, int toUnits)
     {
         if (fromUnits == toUnits) return volume;
         switch (fromUnits)
         {
             case CUFT:
                 if (toUnits == CUYD)
                     return volume * 0.037f;
                 else if (toUnits == ACREFT)
                     return volume * 0.000023f;
                 else if (toUnits == CUMTR)
                     return volume * 0.02832f;
             case CUMTR:
                 if (toUnits == CUYD)
                     return volume * 1.308f;
                 else if (toUnits == ACREFT)
                     return volume * 0.0008123f;
                 else if (toUnits == CUFT)
                     return volume * 35.31f;
             default:
                 return volume;
         }
     }

     public int getNumFractureSets() { return numFractureSets; }

     public void analyzeFractureSets(StsProgressPanel panel)
     {
         progressPanel = panel;
         progressPanel.resetProgressBar();
         progressPanel.setValue(0);
         Runnable runnable = new Runnable()
         {
             public void run()
             {
                 analyzeFractureSets();
             }
         };
         Thread processThread = new Thread(runnable);
         processThread.start();
     }
     public void analyzeVolume(StsProgressPanel panel)
     {
         progressPanel = panel;
         progressPanel.resetProgressBar();
         progressPanel.setValue(0);
         Runnable runnable = new Runnable()
         {
             public void run()
             {
                 analyzeVolume();
             }
         };
         Thread processThread = new Thread(runnable);
         processThread.start();
     }

     static public void main(String[] args)
     {
         StsModel model = StsModel.constructor();
         StsActionManager actionManager = new StsActionManager(model);
         StsVolumeStimulatedWizard volumetricsWizard = new StsVolumeStimulatedWizard(actionManager);
         volumetricsWizard.start();
     }
 }

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

public class StsWellWizard extends StsWizard
{
    int nWells;
    StsPoint[] topPoints = null;
    byte[] types = null;
    boolean projectBound = false;

    static final byte WEBJARFILE = StsWellImport.WEBJARFILE;
    static final byte JARFILE = StsWellImport.JARFILE;
    static final byte BINARYFILES = StsWellImport.BINARYFILES;
    static final byte ASCIIFILES = StsWellImport.ASCIIFILES;
    static final byte LASFILES = StsWellImport.LASFILES;
    static final byte UTFILES = StsWellImport.UTFILES;
    static final byte GEOGRAPHIXFILES = StsWellImport.GEOGRAPHIX;

    public static final byte S2S_WELLS = 0;
    public static final byte GEOGRAPHIX_WELLS = 1;

    static final int NUM_FILE_TYPES = 5;
    public static final byte HEADER = 0;
    public static final byte SURVEY = 1;
    public static final byte LOGS = 2;
    public static final byte TD = 3;
    public static final byte TOPS = 4;

    StsFileType selectFileType = null;
    StsWellHeadDefinition wellheadDef = null;
    StsWellSelect selectWells = null;
    StsFileType defineFileType = null;
    StsHeaderFile defineHeaderFile = null;
    StsOtherFiles defineOtherFiles = null;
    StsColumnDefine defineColumns = null;
    StsDefineUtWells defineUtWells = null;
    StsWellLoad loadWells = null;
    StsWizardStep[] mySteps;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte binaryHorzUnits = StsParameters.DIST_NONE;
    byte binaryVertUnits = StsParameters.DIST_NONE;

    double[][] wellheadXYZs = null;
    String[] wellNames = null;
    boolean[] applyKbs = null;
    int[] selectedWells = null;

    private StsAbstractFile[] wellFiles = new StsFile[0];
    static public final boolean debug = false;

    public StsWellWizard(StsActionManager actionManager)
    {
        super(actionManager, 600, 700);
        StsWellImport.initialize(model);
        StsKeywordIO.initialize(model);
        addSteps(true);
        StsWellImport.setReloadAscii(false);
    }

    public StsWellWizard(StsActionManager actionManager, boolean displayReloadAscii)
    {
        super(actionManager);
        StsWellImport.initialize(model);
        StsKeywordIO.initialize(model);
        addSteps(displayReloadAscii);
    }

    private void addSteps(boolean displayReloadAscii)
    {
        selectFileType = new StsFileType(this);
        wellheadDef = new StsWellHeadDefinition(this);
        selectWells = new StsWellSelect(this, displayReloadAscii);
        defineUtWells = new StsDefineUtWells(this);
        loadWells = new StsWellLoad(this);
        mySteps = new StsWizardStep[] {selectFileType, selectWells, wellheadDef, defineUtWells, loadWells};
        addSteps(mySteps);
    }

    public boolean start()
    {

        System.runFinalization();
        System.gc();
        dialog.setTitle("Load Wells");
        initialize();
        this.disableFinish();
        return super.start();
    }

    public void initialize()
    {
        loadWells.setWellFactory(new StsWellFactory());
        hUnits = model.getProject().getXyUnits();
        vUnits = model.getProject().getDepthUnits();
    }

    public boolean end()
    {
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectWells)
        {
            StsWellImport.setCurrentDirectory(selectWells.getSelectedDirectory());
            StsWellImport.setApplyDatumShift(selectWells.panel.getApplyDatumShift());
            String[] selectedWellnames = selectWells.getSelectedWells();
            readHeaders();
            if(selectedWellnames.length == 0)
            {
                new StsMessage(dialog, StsMessage.INFO, "Must select at least one well");
                return;
            }
            constructAvailableFiles();

            if(getFileType() == GEOGRAPHIX_WELLS)
                gotoStep(wellheadDef);
            else if(StsWellImport.selectedOfType(UTFILES, selectedWellnames))
            {
                defineUtWells.setSelectedWellnames(selectedWellnames);
                gotoStep(defineUtWells);
            }
            else
            {
                loadWells.setSelectedWellnames(selectedWellnames);
                gotoStep(loadWells);
            }
        }
        else if(currentStep == wellheadDef)
        {
            String[] selectedWellnames = selectWells.getSelectedWells();
            loadWells.setSelectedWellnames(selectedWellnames);
            if (StsWellImport.selectedOfType(UTFILES, selectedWellnames))
            {
                defineUtWells.setSelectedWellnames(selectedWellnames);
                gotoStep(defineUtWells);
            }
            else
            {
                gotoStep(loadWells);
            }
        }
        else if (currentStep == defineUtWells)
        {
            StsPoint[] topPoints = defineUtWells.getTopHoleLocations();
            byte[] types = defineUtWells.getTypes();
            StsWellImport.setTopHoleLocations(topPoints, types);
            String[] selectedWellnames = selectWells.getSelectedWells();
            loadWells.setSelectedWellnames(selectedWellnames);
            gotoStep(loadWells);
        }
        else
        {
            gotoNextStep();
        }
    }

    public byte getFileType() {  return selectFileType.getType(); }
    public void finish()
    {
        success = true;
        super.finish();
    }

    public void loadMoreWells()
    {
        gotoFirstStep();
        selectWells.panel.resetSelection();
    }
    public boolean getArchiveIt() { return selectWells.panel.getArchiveIt(); }
    public float getDatumShift() { return selectWells.panel.getDatumShift(); }

    public void constructAvailableFiles()
    {
        StsWellImport.initializeWellFilenameSets();
        // User overriding standard naming conventions.
        if(StsWellImport.getOverrideFilter())
        {
            // Manually construct the well sets.
           StsWellImport.constructWellFilenameSets(selectWells.panel.getUserFilter(), LASFILES);
        }
        // Webstart Download Jar File
        else if (Main.isWebStart && Main.isJarFile)
        {
            StsWellImport.constructWellFilenameSets(WEBJARFILE);
        }

        // Load from jar files
        else if (Main.isJarFile)
        {
            // Jar Files
            StsWellImport.constructWellFilenameSets(JARFILE);
        }
        // Load from ASCII/Binary files
        else if(selectFileType.getType() == GEOGRAPHIX_WELLS)
        {
            StsWellImport.constructWellFilenameSets(GEOGRAPHIXFILES);             
        }
        else
        {
            StsWellImport.constructWellFilenameSets(LASFILES);
            StsWellImport.constructWellFilenameSets(ASCIIFILES);
            StsWellImport.constructWellFilenameSets(UTFILES);
        }

        // Compress the Wellset list.
        StsWellImport.compressWellSet();
    }

    public boolean addFile(StsAbstractFile file, String outputDirectory)
    {
        if(selectFileType.getType() == GEOGRAPHIX_WELLS)
        {
            if(wellFiles.length > 0)
            {
                new StsMessage(frame,StsMessage.WARNING,"Currently only able to load one GeoGraphix file at a time.");
                return false;
            }
        }
        wellFiles = (StsFile[]) StsMath.arrayAddElement(wellFiles, file);
        return true;
    }

    public void removeFile(StsAbstractFile file)
    {
        wellFiles = (StsAbstractFile[]) StsMath.arrayDeleteElement(wellFiles, file);
    }

    public void removeFiles()
    {
        if(wellFiles == null) return;
        for(int i = 0; i<wellFiles.length; i++)
            wellFiles = (StsFile[]) StsMath.arrayDeleteElement(wellFiles, wellFiles[i]);

        wellFiles = null;
    }
    public StsAbstractFile[] getSelectedFiles() { return wellFiles; }

    public String[] getSelectedWells()
    {
        if((!projectBound) || (selectedWells == null))
            return wellNames;
        else
        {
            String[] names = new String[selectedWells.length];
            for(int i=0; i<selectedWells.length; i++)
                names[i] = wellNames[selectedWells[i]];
            return names;
        }
    }

    public JDialog processFileDialog(StsProgressPanel panel)
    {
        JDialog dialog;
        dialog = new JDialog(frame, "Loading Realtime..." + getName(), false);
        dialog.setPreferredSize(new Dimension(300,50));
        dialog.add(panel);
        dialog.pack();
        dialog.setLocationRelativeTo(frame);
        dialog.setVisible(true);
        return dialog;
    }

    public boolean readHeaders()
    {
        String[] wellnames = null;
        wellNames = null;

        for(int i=0; i<wellFiles.length; i++)
        {
            if(getFileType() == GEOGRAPHIX_WELLS)
            {
                StsGeographixKeywordIO.readFileHeader(selectWells.getSelectedDirectory(), wellFiles[i].getFilename());
                wellnames = StsGeographixKeywordIO.getFileHeader().getWellnames();
                wellheadXYZs = StsGeographixKeywordIO.getFileHeader().getWellheads();
                applyKbs = StsGeographixKeywordIO.getFileHeader().getApplyKbCorrection();
            }
            else
            {
                wellnames = new String[] {wellFiles[i].name};
            }
            if(wellnames != null)
                wellNames = (String[]) StsMath.arrayAddArray(wellNames, wellnames);
        }
        return true;
    }
    
    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsWellWizard wellWizard = new StsWellWizard(actionManager);
        wellWizard.start();
    }
    public int getNumberRowsToSkip(int type)
    {
         if(type != StsWellWizard.HEADER)
            return defineOtherFiles.panel.getNumRowsToSkip(type);
        else
            return defineHeaderFile.panel.getNumberRowsToSkip();
    }
    public String[] getLogFiles()
    {
          return defineOtherFiles.panel.getLogFiles();
    }
    public String getDefaultFilePath()
    {
        return defineHeaderFile.panel.getHeaderFilePath();
    }

    public int[] getColOrder(byte type)
    {
        return defineColumns.panel.getColOrder(type);
    }
    public String getFilename(int type)
    {
        switch(type)
        {
            case HEADER:
                return defineHeaderFile.panel.getHeaderFile();
            case SURVEY:
                return defineOtherFiles.panel.getSurveyFile();
            case TOPS:
                return defineOtherFiles.panel.getTopsFile();
            case TD:
                return defineOtherFiles.panel.getTdFile();
            default:
                new StsMessage(frame, StsMessage.WARNING, "Unsupported file type");
        }
        return null;
    }

    public double getX()
    {
        if(wellheadDef.getSelectedIndices() == null)
            return 0.0f;
        if(selectedWells != null)
            return wellheadXYZs[selectedWells[wellheadDef.getSelectedIndices()[0]]][0];
        else
            return wellheadXYZs[wellheadDef.getSelectedIndices()[0]][0];
    }
    public double getY()
    {
        if(wellheadDef.getSelectedIndices() == null)
            return 0.0f;
        if(selectedWells != null)
            return wellheadXYZs[selectedWells[wellheadDef.getSelectedIndices()[0]]][1];
        else
            return wellheadXYZs[wellheadDef.getSelectedIndices()[0]][1];
    }
    public double getZ()
    {
        if(wellheadDef.getSelectedIndices() == null)
            return 0.0f;
        if(selectedWells != null)
            return wellheadXYZs[selectedWells[wellheadDef.getSelectedIndices()[0]]][2];
        else
            return wellheadXYZs[wellheadDef.getSelectedIndices()[0]][2];
    }
    public boolean getApplyKb()
    {
        if(wellheadDef.getSelectedIndices() == null)
            return false;
        if(selectedWells != null)
            return applyKbs[selectedWells[wellheadDef.getSelectedIndices()[0]]];
        else
            return applyKbs[wellheadDef.getSelectedIndices()[0]];
    }
    public void setX(double val)
    {
        int[] indices = wellheadDef.getSelectedIndices();
        if(indices == null)
            return;
        for(int i=0; i<indices.length; i++)
        {
            if(selectedWells != null)
                wellheadXYZs[selectedWells[indices[i]]][0] = val;
            else
                wellheadXYZs[indices[i]][0] = val;
        }
    }
    public void setY(double val)
    {
        int[] indices = wellheadDef.getSelectedIndices();
        if(indices == null)
            return;
        for(int i=0; i<indices.length; i++)
        {
            if(selectedWells != null)
                wellheadXYZs[selectedWells[indices[i]]][1] = val;
            else
                wellheadXYZs[indices[i]][1] = val;
        }
    }
    public void setZ(double val)
    {
        int[] indices = wellheadDef.getSelectedIndices();
        if(indices == null)
            return;
        for(int i=0; i<indices.length; i++)
        {
            if(selectedWells != null)
                wellheadXYZs[selectedWells[indices[i]]][2] = val;
            else
                wellheadXYZs[indices[i]][2] = val;
        }
    }
    public void setApplyKb(boolean val)
    {
        int[] indices = wellheadDef.getSelectedIndices();
        if(indices == null)
            return;
        for(int i=0; i<indices.length; i++)
        {
            if(selectedWells != null)
                applyKbs[selectedWells[indices[i]]] = val;                
            else
                applyKbs[indices[i]] = val;
        }
    }
    public void setIsProjectBound(boolean value)
    {
        projectBound = value;
        if(projectBound)
        {
            StsProject project = model.getProject();
            selectedWells = new int[wellNames.length];
            double projectXOrigin = project.getXOrigin();
            double projectYOrigin = project.getYOrigin();
            double xMin = project.getXOrigin() - wellheadDef.panel.getHalo();
            double yMin = project.getYOrigin() - wellheadDef.panel.getHalo();
            double xMax = xMin + (project.getXMax() - project.getXMin()) + wellheadDef.panel.getHalo();
            double yMax = yMin + (project.getYMax() - project.getYMin()) + wellheadDef.panel.getHalo();
            StsBoundingBox unrotatedProjectBoundingBox = project.unrotatedBoundingBox;
            int cnt = 0;
            for(int i=0; i<wellNames.length; i++)
            {
                float x = (float)(wellheadXYZs[i][0] - projectXOrigin);
                float y = (float)(wellheadXYZs[i][1] - projectYOrigin);
                if(unrotatedProjectBoundingBox.isInsideXY(x, y))
                    selectedWells[cnt++] = i;
                //if((wellheadXYZs[i][0] > xMin) && (wellheadXYZs[i][0] < xMax) && (wellheadXYZs[i][1] > yMin) && (wellheadXYZs[i][1] < yMax))
            }
            selectedWells = (int[])StsMath.trimArray(selectedWells,cnt);
        }
        else
            selectedWells = null;
    }
    public boolean getIsProjectBound() { return projectBound;}
}

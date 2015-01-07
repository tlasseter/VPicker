//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Export;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.Simulation.*;
import com.Sts.UI.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.text.*;
import java.util.*;

public class StsEclipseOutput
{
    StsModel model;
    public StsStatusUI status;
    StsGridBoundingBox modelBoundingBox;
    StsBuiltModel builtModel;
    StsEclipseModel eclipseModel;
    StsZone[] zones;
    private StsBlock[] blocks;
    int nFaultGrids;
    ArrayList<StsFaultGrid> faultGrids = new ArrayList<StsFaultGrid>();
    TreeSet<StsBlock.BlockCellColumn.GridCell> debugCellList;
    EclipseWriter dataWriter;
    EclipseWriter gridWriter;
    EclipseWriter editWriter;
    EclipseWriter scheduleWriter;
    public EclipseWriter outputWriter;
    NumberFormat numberFormat;
    StsGridDefinition grid;
    // int rowMin, rowMax, colMin, colMax; // replace these with a gridRectangle (transient!)
    int nCellRows, nCellCols;
    int nTotalCells;
    double xMin, yMin, xInc, yInc;
    double zMin, zMax;
    int nCelllLayers;
    int nCorrectEclipseRows;

    double kxMultiplier = 1.0;
    double kyMultiplier = 1.0;
    double kzMultiplier = .01;

    ArrayList<StsNNC> nncList;

    public float progress = 0.0f;

    public boolean hasProperties = false;

    protected int numberMergedCells = 0;
    private int numberTruncatedCells = 0;
    public float mergeCellMinFraction = 0.5f;
    protected float actualMergedCellMinFraction = 1.0f;
    protected int numberUnmergedSmallCells = 0;

    public boolean adjustCells = true;

    static public StsEclipseOutput eclipseOutputInstance;
    
    static public float tranFactor = 117.2771f;
    static public float tranMultiplier = 1.0f/tranFactor;
    static int fieldsPerLine = 10;
    static final boolean debug = true;
    static final float nullValue = StsParameters.nullValue;
    static float fillZ = 0.0f;

    /** used in time-constant calculations */
    static float waterCompressibility = 4.67E-5f;
    static float minPoreVolume = 0.0f;

    static public final byte CELL_PLUS_DX = StsBlock.CELL_PLUS_DX;
    static public final byte CELL_PLUS_DY = StsBlock.CELL_PLUS_DY;
    static public final byte CELL_PLUS_DZ = StsBlock.CELL_PLUS_DZ;
    static public final byte CELL_MINUS_DX = StsBlock.CELL_MINUS_DX;
    static public final byte CELL_MINUS_DY = StsBlock.CELL_MINUS_DY;
    static public final byte CELL_MINUS_DZ = StsBlock.CELL_MINUS_DZ;

    static final float largeFloat = StsParameters.largeFloat;

    static final String model2constant = "2constant";
    static final String model3constant = "3constant";
    static final String model3ortho = "3ortho";
    static final String model1 = "1";
    static final String model2 = "2";
    static final String model3 = "3";
    static final String model4 = "4";
    static final String model5 = "5";
    static final String model6 = "6";

    static final String[] modelNameStrings = new String[] { model2constant, model3constant, model3ortho, model1, model2, model3, model4, model5, model6,  };
    static final int[] modelNameNumbers = new int[] { 21, 31, 32, 1, 2, 3, 4, 5, 6 };

    private StsEclipseOutput(StsModel model)
    {
        this.model = model;
        builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        zones = builtModel.getZones();
        blocks = builtModel.getBlocks();
        modelBoundingBox = builtModel.modelBoundingBox;
        xMin = 0.0;
        yMin = 0.0;

        xInc = modelBoundingBox.xInc;
        yInc = modelBoundingBox.yInc;
        zMin = modelBoundingBox.zMin;
        zMax = modelBoundingBox.zMax;

        StsBlock.minPoreVolume = minPoreVolume;
    }

    static public StsEclipseOutput getInstance(StsModel model)
    {
        if(eclipseOutputInstance == null)
            eclipseOutputInstance = new StsEclipseOutput(model);
        return eclipseOutputInstance;
    }

    public void constructLayers()
    {
        nCelllLayers = 0;
        for (StsZone zone : zones)
        {
            StsObjectRefList subZones = zone.getSubZones();
            int nSubZones = subZones.getSize();
            for (int sz = 0; sz < nSubZones; sz++)
            {
                StsSubZone subZone = (StsSubZone) subZones.getElement(sz);
                nCelllLayers += subZone.getNLayers();
            }
        }
    }

    /*
        public void initializeFaultGrids()
        {
            StsObject[] sections = model.getObjectList(StsSection.class);
            int nSections = sections.length;
            faultGrids = new StsFaultGrid[nSections];
            int nFaultGrids = 0;
            for(int n = 0; n < nSections; n++)
            {
                StsSection section = (StsSection)sections[n];
                if(section.isFault())
                    faultGrids[nFaultGrids++] = new StsFaultGrid(section);
            }
            faultGrids = (StsFaultGrid[])StsMath.trimArray(faultGrids, nFaultGrids);
        }
    */
    /** Blocks are laid out in sequential row-order with a null row between blocks to keep them from being interconnected automatically.
     *  For each block, set the Eclipse starting cellRowMin and find the longest row which determines the max number of columns in the Eclipse grid.
     */
    private void constructBlockLayout()
    {
        int nBlocks = blocks.length;
        nCellRows = 0;
        nCellCols = 0;
        for (int n = 0; n < nBlocks; n++)
        {
            StsBlock block = blocks[n];
            block.eclipseRowMin = nCellRows + 1;
            nCellRows += block.nCellRows;
            nCellCols = Math.max(nCellCols, block.nCellCols);
            // insert a null row between blocks
            if(n < nBlocks-1) nCellRows++;
            block.setLayerRange(nCelllLayers);
        }
        nTotalCells = nCellRows * nCellCols * nCelllLayers;
    }

    private void fillBlockCellGrids()
    {
        for (StsBlock block : blocks)
            block.fillInsideCells();
    }

    /*
        private void addBlockGridCells()
        {
            int nBlocks = blocks.length;
            for(int n = 0; n < nBlocks; n++)
                addBlockGridCells(blocks[n]);
        }

        private void addBlockGridCells(StsBlock block)
        {
            StsBlockGrid[] blockGrids = block.getBlockGrids();
            for(StsBlockGrid blockGrid : blockGrids)
            {
                StsCellTypeGrid cellGrid = blockGrid.getEdgeLoop().getCellTypeGrid();
                if(cellGrid == null) continue;
                block.addBlockGridCells(cellGrid);
            }
        }
    */
    /** For each faulted section, construct a faultGrid which has polygons on each side */
    private boolean addFaultGridPolygons()
    {
        StsObject[] sections = model.getObjectList(StsSection.class);
        int nSections = sections.length;
        for (int n = 0; n < nSections; n++)
        {
            StsSection section = (StsSection) sections[n];
            if(!section.isFault()) continue;
            faultGrids.add(new StsFaultGrid(section));
        }
        return true;
    }


    public boolean createEclipseFiles(String stemname)
    {
        try
        {
            StsDecimalFormat floatFormat = new StsDecimalFormat(8);
            initializeConversionFactors();
            constructEclipseModel();
            if(stemname.equals("default"))
                stemname = builtModel.name + "." + eclipseModel.name;
            dataWriter = new EclipseWriter(eclipseModel.inputPathname + stemname + ".DATA", floatFormat, 5);
            gridWriter = new EclipseWriter(eclipseModel.inputPathname + stemname + ".GRID", floatFormat, 5);
            editWriter = new EclipseWriter(eclipseModel.inputPathname + stemname + ".EDIT", floatFormat, 5);
            scheduleWriter = new EclipseWriter(eclipseModel.inputPathname + stemname + ".SCH", floatFormat, 5);
            outputWriter = new EclipseWriter(eclipseModel.inputPathname + stemname + ".OUT", floatFormat, 5);
            // clearBlockProperties();
            constructProperties();
            constructLayers();
            constructBlockLayout();

            setProgress(5.0f);
            int nBlocks = blocks.length;
            boolean persistent = false;
            for (int n = 0; n < nBlocks; n++)
            {
                StsBlock block = blocks[n];
                block.constructCellColumns(nCelllLayers);
                block.constructBasicProperties(zones, eclipseModel);
                status.setProgress(progress);
            }
            setBlockPropertyRanges();
            StsBlock.addCellColumnPolygons(model, blocks);
            setProgress(10.0f);
            addFaultGridPolygons();
            setProgress(20.0f);
            StsBlock.constructBlockCellColumns(model, blocks);
            StsZone.constructLayerGrids(zones);
            setProgress(25.0f);
            StsBlock.constructLayerColumnGrids(model, blocks);
            StsBlock.constructEdgeGridCells(model, blocks);
            setProgress(35.0f);
            StsBlock.computeCellProperties(zones, blocks, eclipseModel);
            initializeOutput();
            outputRunSpec(dataWriter);
            setProgress(44.0f);
//            computePoreVolume(4);
            computeDepths(4);
//            computeBlockTranX(4);
//            computeBlockTranY(4);
//            computeBlockTranZ(4);
            constructNNCs(4);
            createTruncatedCellListComputeFractions();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "createEclipseFile", e);
            return false;
        }
    }

    protected void runOutputEclipseFile()
    {
        StsToolkit.runRunnable
        (
            new Runnable()
            {
                public void run()
                {
                    outputEclipseFile();
                }
            }
        );
    }
    public boolean outputEclipseFile()
    {
        try
        {
            outputTruncationInfo();
            merge();
            computePoreVolume(4);
            // hackRemoveSmallPoreVolumes();
            computeActiveCellIndexMaps(2);
            // countActiveCells();
            computeBlockTranX(4);
            computeBlockTranY(4);
            computeBlockTranZ(4);
            if (!outputCoord(gridWriter, 4)) return false;
            if (!outputZCorn(gridWriter, 4)) return false;
            if(!outputBooleanProperty(gridWriter, eclipseModel.actnum,  2)) return false;
            if(!outputProperty(gridWriter, eclipseModel.porosity,  2)) return false;
            if(!outputProperty(gridWriter, eclipseModel.permX,  2)) return false;
            if(!outputProperty(gridWriter, eclipseModel.permY,  2)) return false;
            if(!outputProperty(gridWriter, eclipseModel.permZ,  2)) return false;
            if(!outputNNCs(gridWriter, 2)) return false;
            
            if(!outputProperty(editWriter, eclipseModel.depth, 2)) return false;
            if(!outputProperty(editWriter, eclipseModel.poreVolume, 2)) return false;
            if(!outputProperty(editWriter, eclipseModel.tranX, 2, tranMultiplier)) return false;
            if(!outputProperty(editWriter, eclipseModel.tranY, 2, tranMultiplier)) return false;
            if(!outputProperty(editWriter, eclipseModel.tranZ, 2, tranMultiplier)) return false;
            if(!outputSchedule(scheduleWriter)) return false;
            status.setProgress(100);
            dataWriter.close();
            gridWriter.close();
            editWriter.close();
            scheduleWriter.close();
            outputWriter.close();
            outputDataFiles();
        //    clearTransients();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "createEclipseFile", e);
            return false;
        }
    }

    static public void cancel()
    {
        clearTransients();
    }

    static public void clearTransients()
    {
        for(StsBlock block : eclipseOutputInstance.blocks)
            block.clearTransientArrays();
        eclipseOutputInstance = null;
    }

    private void outputDataFiles()
    {
        outputDataFile(eclipseModel.porosity);
        outputDataFile(eclipseModel.permX);
        outputDataFile(eclipseModel.permY);
        outputDataFile(eclipseModel.permZ);
        outputDataFile(eclipseModel.poreVolume);
        outputDataFile(eclipseModel.depth);
        outputDataFile(eclipseModel.actnum);
        outputDataFile(eclipseModel.tranX);
        outputDataFile(eclipseModel.tranY);
        outputDataFile(eclipseModel.tranZ);
        outputDataFile(eclipseModel.indexMap);
    }

    private void constructEclipseModel()
    {
        eclipseModel = new StsEclipseModel(builtModel);
    }

    private void outputDataFile(StsPropertyType propertyType)
    {
        for(StsBlock block : blocks)
            eclipseModel.adjustRange(propertyType, block);
        propertyType.rangeChanged();
        for(StsBlock block : blocks)
            eclipseModel.outputDataFile(propertyType, block);
    }

    private void setProgress(float progress)
    {
        this.progress = progress;
        status.setProgress(progress);
    }

    private void initializeConversionFactors()
    {
    }

    private void initializeOutput()
    {
        nCorrectEclipseRows = 0;
        for(StsBlock block : blocks)
            nCorrectEclipseRows += block.nCellRows;
        // add separator rows between blocks
        nCorrectEclipseRows += blocks.length - 1;
    }
/*
    private void clearBlockProperties()
    {
        for (StsBlock block : blocks)
            block.clearProperties();

    }
*/
    private void constructProperties()
    {
        String projectName = model.getProject().getName();
        int nModelNameStrings = modelNameStrings.length;
        int modelNumber = -1;

        for(int n = 0; n < nModelNameStrings; n++)
        {
            if(StsStringUtils.stringContainsString(projectName, modelNameStrings[n]))
            {
                modelNumber = modelNameNumbers[n];
                break;
            }
        }
        if(modelNumber == -1)
        {
            StsException.systemError(this, "constructProperties", "Project name " + projectName + " does not contain model number (1,2,3,4,5,or 6)\n" +
                "No properties will be generated");
            return;
        }

        switch(modelNumber)
        {
            case 21: // 2a
                constructConstantProperties(0.10f, 100.0f, 100.0f);
                // constructLayerProperties();
                // constructFaultBlockProperties(0.1f, 1.0f, 1.0f);
                break;
            case 31: // 3a
                constructConstantProperties(0.10f, 100.0f, 100.0f);
                break;
            case 32:
                constructLayerProperties();
                constructFaultBlockProperties(0.1f, 1.0f, 1.0f);
                break;
            case 1: 
            case 2: // 2
                constructLayerProperties();
                break;
            case 3:
            case 5:
                constructLayerProperties();
                constructFaultBlockProperties(0.1f, 1.0f, 1.0f);
                break;
            case 4:
                constructAnisotropicLayerProperties();
                constructFaultBlockProperties(0.1f, 5.0f, 10.0f);
                break;
            case 6:
                constructLayer6Properties();
                constructFaultBlockProperties(0.1f, 1.0f, 1.0f);
                break;
        }
    }

    private void constructConstantProperties(float porosity, float permH, float permV)
    {
        for(StsZone zone : zones)
        {
            ArrayList<StsPropertyVolume> zoneProperties = zone.getZoneProperties();

            StsPropertyVolume porosityVolume = new StsPropertyVolume(StsPropertyType.porosity, zone);
            porosityVolume.setConstant(porosity);
            zoneProperties.add(porosityVolume);

            StsPropertyVolume permxVolume = new StsPropertyVolume(StsPropertyType.permX, zone);
            permxVolume.setConstant(permH);
            zoneProperties.add(permxVolume);

            StsPropertyVolume permyVolume = new StsPropertyVolume(StsPropertyType.permY, zone);
            permyVolume.setConstant(permH);
            zoneProperties.add(permyVolume);

            StsPropertyVolume permzVolume = new StsPropertyVolume(StsPropertyType.permZ, zone);
            permzVolume.setConstant(permV);
            zoneProperties.add(permzVolume);
        }
    }

    private void constructLayerProperties()
    {
        float[] porosities = new float[]{.25f, .1f, .1f, .25f, .25f, .1f, .25f, .1f, .25f, .1f};
        float[] perm = new float[]{1000f, 1f, 1f, 1000f, 1000f, 1f, 1000f, 1f, 1000f, 1f};
        boolean persistent = false;
        for(StsZone zone : zones)
        {
            ArrayList<StsPropertyVolume> zoneProperties = zone.getZoneProperties();

            StsPropertyVolume porosityVolume = new StsPropertyVolume(StsPropertyType.porosity, zone);
            porosityVolume.setVariableZ(porosities);
            zoneProperties.add(porosityVolume);

            StsPropertyVolume permxVolume = new StsPropertyVolume(StsPropertyType.permX, zone);
            permxVolume.setVariableZ(perm);
            zoneProperties.add(permxVolume);

            StsPropertyVolume permyVolume = new StsPropertyVolume(StsPropertyType.permY, zone);
            permyVolume.setVariableZ(perm);
            zoneProperties.add(permyVolume);

            StsPropertyVolume permzVolume = new StsPropertyVolume(StsPropertyType.permZ, zone);
            permzVolume.setVariableZ(perm);
            zoneProperties.add(permzVolume);
        }
    }

    private void setBlockPropertyRanges()
    {

        ArrayList<StsPropertyVolume> propertyVolumes = eclipseModel.getBlockPropertyVolumes();
        for(StsPropertyVolume propertyVolume : propertyVolumes)
        {
            if(propertyVolume.name == StsPropertyType.permX.name)
                propertyVolume.setRange(0.0f, 1000.0f);
            else if(propertyVolume.name == StsPropertyType.permY.name)
                propertyVolume.setRange(0.0f, 1000.0f);
            else if(propertyVolume.name == StsPropertyType.permZ.name)
                propertyVolume.setRange(0.0f, 1000.0f);
            else if(propertyVolume.name == StsPropertyType.porosity.name)
                propertyVolume.setRange(0.0f, 1.0f);
        }
    }

    private void constructFaultBlockProperties(float porosity, float permH, float permV)
    {
        StsBlock block = blocks[1];
        String outputPathname = eclipseModel.outputPathname;
        StsPropertyVolume porosityVolume = new StsPropertyVolume(eclipseModel.porosity, block, outputPathname);
        porosityVolume.setConstant(porosity);
        eclipseModel.addBlockPropertyVolume(porosityVolume);

        StsPropertyVolume permxVolume = new StsPropertyVolume(eclipseModel.permX, block, outputPathname);
        permxVolume.setConstant(permH);
        eclipseModel.addBlockPropertyVolume(permxVolume);

        StsPropertyVolume permyVolume = new StsPropertyVolume(eclipseModel.permY, block, outputPathname);
        permyVolume.setConstant(permH);
        eclipseModel.addBlockPropertyVolume(permyVolume);

        StsPropertyVolume permzVolume = new StsPropertyVolume(eclipseModel.permZ, block, outputPathname);
        permzVolume.setConstant(permV);
        eclipseModel.addBlockPropertyVolume(permzVolume);
    }

    private void constructAnisotropicLayerProperties()
    {
        float[] porosities = new float[]{.25f, .1f, .1f, .25f, .25f, .1f, .25f, .1f, .25f, .1f};
        float[] permxy = new float[]{1000f, 10f, 10f, 1000f, 1000f, 10f, 1000f, 10f, 1000f, 10f};
        float[] permz = new float[]{10f, .1f, .1f, 10f, 10f, .1f, 10f, .1f, 10f, .1f};

        for(StsZone zone : zones)
        {
            ArrayList<StsPropertyVolume> zoneProperties = zone.getZoneProperties();

            StsPropertyVolume porosityVolume = new StsPropertyVolume(StsPropertyType.porosity, zone);
            porosityVolume.setVariableZ(porosities);
            zoneProperties.add(porosityVolume);

            StsPropertyVolume permxVolume = new StsPropertyVolume(StsPropertyType.permX, zone);
            permxVolume.setVariableZ(permxy);
            zoneProperties.add(permxVolume);

            StsPropertyVolume permyVolume = new StsPropertyVolume(StsPropertyType.permY, zone);
            permyVolume.setVariableZ(permxy);
            zoneProperties.add(permyVolume);

            StsPropertyVolume permzVolume = new StsPropertyVolume(StsPropertyType.permZ, zone);
            permzVolume.setVariableZ(permz);
            zoneProperties.add(permzVolume);
        }
    }

    private void constructLayer6Properties()
    {
        float[] porosities = new float[]{.25f, .1f, .25f, .1f, .25f, .1f, .25f, .1f, .25f, .1f};
        float[] perms = new float[]{1000f, 1f, 1000f, 1f, 1000f, 1f, 1000f, 1f, 1000f, 1f};
        boolean persistent = false;
        int nZones = zones.length;
        int nLayer = 0;
        for(int n = 0; n < nZones; n++)
        {
            StsZone zone = zones[n];
            int nSubZones = zone.getNSubZones();
            ArrayList<StsPropertyVolume> zoneProperties = zone.getZoneProperties();
            StsPropertyVolume porosityVolume = new StsPropertyVolume(StsPropertyType.porosity, zone);
            porosityVolume.setVariableZ(porosities, nLayer, nSubZones);
            zoneProperties.add(porosityVolume);
            StsPropertyVolume permxVolume = new StsPropertyVolume(StsPropertyType.permX, zone);
            permxVolume.setVariableZ(perms, nLayer, nSubZones);
            zoneProperties.add(permxVolume);
            StsPropertyVolume permyVolume = new StsPropertyVolume(StsPropertyType.permY, zone);
            permyVolume.setVariableZ(perms, nLayer, nSubZones);
            zoneProperties.add(permyVolume);
            StsPropertyVolume permzVolume = new StsPropertyVolume(StsPropertyType.permZ, zone);
            permzVolume.setVariableZ(perms, nLayer, nSubZones);
            zoneProperties.add(permzVolume);
            nLayer += nSubZones;
        }
    }

    public boolean outputRunSpec(EclipseWriter eclipseWriter)
    {
        eclipseWriter.println("DIMENS");
        eclipseWriter.println("-- NX   NY   NZ");
        eclipseWriter.println(nCellCols + " " + nCellRows + " " + nCelllLayers);
        eclipseWriter.endRecord();
        return true;
    }

    public boolean outputCoord(EclipseWriter eclipseWriter, float dProgress)
    {
        int nBlock, eclipseRow, col, blockRow;
        try
        {
            eclipseWriter.initialize(6);
            StsMessageFiles.logMessage("Writing COORD section...");
            eclipseWriter.println("COORD"); // COORD Keyword

            double[] xyzTop = new double[3];
            double[] xyzBot = new double[3];
            xyzTop[2] = zMin;
            xyzBot[2] = zMax;
            xyzTop[1] = yMin;
            xyzBot[1] = xyzTop[1];
            int nBlocks = blocks.length;
            dProgress /= nBlocks;
            eclipseRow = 1;
            for (nBlock = 0; nBlock < nBlocks; nBlock++, progress += dProgress)
            {
                StsBlock block = blocks[nBlock];
                int blockIndex = block.getIndex();
                int blockRowMin = block.getBlockRowMin();
                int blockRowMax = block.getBlockRowMax();
                eclipseWriter.comment("block: " + blockIndex + " block grid rows " + block.nCellRows + " to 0" + " nCols " + nCellCols);
                for (blockRow = blockRowMax; blockRow >= blockRowMin; blockRow--, eclipseRow++)
                {
                    eclipseWriter.comment("COORD block: " + blockIndex +  " block row: " + blockRow + " Eclipse row " + eclipseRow );
                    xyzTop[0] = 0.0f;
                    xyzBot[0] = 0.0f;
                    for (col = 0; col <= nCellCols; col++)
                    {
                        eclipseWriter.print(xyzTop);
                        eclipseWriter.print(xyzBot);
                        xyzTop[0] += xInc;
                        xyzBot[0] += xInc;
                    }

                    xyzTop[1] += yInc;
                    xyzBot[1] += yInc;
                }
                status.setProgress(progress);
            }
            eclipseWriter.endRecord();
            if(debug)
            {
                int nRows = 0;
                for(StsBlock block : blocks)
                    nRows += block.nCellRows+1;
                int correctCount = 6*nRows*(nCellCols+1);
                debugCheckCount(eclipseWriter, "outputCoord", " ", correctCount);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEclipseOutput.outputCoord() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /** Writes out the cell geometry.  There are nLayer layers of cells.  For each layer we write out the top and bottom
     *  planes of data
     * @param eclipseWriter PrintWriter for data output
     * @param dProgress progress % change for this operation
     * @return  true if successful
     */
    public boolean outputZCorn(EclipseWriter eclipseWriter, float dProgress)
    {
        eclipseWriter.initialize(fieldsPerLine);
        StsMessageFiles.logMessage("Writing ZCORN section...");
        eclipseWriter.println("ZCORN");

        dProgress /= nCelllLayers;
        for (int nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
        {
            if(!writeZCornPlane(eclipseWriter, nLayer, "top")) return false;
            if(!writeZCornPlane(eclipseWriter, nLayer+1, "bottom")) return false;
            status.setProgress(progress);
        }
        eclipseWriter.endRecord();
        if(debug)
        {
            int correctCount = nCelllLayers *nCorrectEclipseRows*nCellCols*8;
            debugCheckCount(eclipseWriter, "outputZCorn", "", correctCount);
        }
        return true;
    }

    /** Writes a single plane of data which might be the top or bottom of a layer of cells.
     *  Rows in the plane are written out back to front in decreasing S2S row order, increasing Eclipse row order, decreasing in Y value.
     *  For each row in the plane, we write out the back row of points for the cells and then the front row of points.  For each cell,
     *  we write out the left and right points.  So we have a total of four points written out for each cell.
     *
     * @param eclipseWriter PrintWriter to which data is written
     * @param nLayer current cell layer
     * @param layerPosition position of this plane in the layer, typically "top" or "bottom"
     * @return true if successful
     */
    private boolean writeZCornPlane(EclipseWriter eclipseWriter, int nLayer, String layerPosition)
    {
        int nBlock = -1, row = -1;
        int blockRow;
        int nBlocks = blocks.length;
        int eclipseCellRow = 1;
        try
        {
            for (nBlock = 0; nBlock < nBlocks; nBlock++)
            {
                StsBlock block = blocks[nBlock];
                int blockIndex = block.getIndex();
                float[][] layerGrid = block.getLayerGrid(nLayer);
                int nBlockCellRows = block.nCellRows;
                int nBlockCellCols = block.nCellCols;
                int rowMax = block.getBlockRowMax();
                eclipseWriter.comment("block: " + blockIndex + " layer: " + nLayer + " " + " nBlockCellRows: " + nBlockCellRows + " nBlockCellCols" + nBlockCellCols);
                fillZ = getFirstFillZ(layerGrid, nBlockCellRows, nBlockCellCols);
                eclipseWriter.comment("ZCORN layer " + nLayer + " " + layerPosition + " default fill Z: " + fillZ);
                for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--, eclipseCellRow++)
                {
                    // for each cell in block, write back left and right corner z values
                    rowComment(eclipseWriter, "ZCORN", layerPosition + " back", blockIndex, nLayer, row, row, eclipseCellRow);
                    if(!writeZCornRow(eclipseWriter, layerGrid, nBlock, nLayer, blockRow+1, nBlockCellCols)) return false;
                    // for each cell in block, write front left and right corner z values
                    rowComment(eclipseWriter, "ZCORN", layerPosition + " front", blockIndex, nLayer, row , row, eclipseCellRow);
                    if(!writeZCornRow(eclipseWriter, layerGrid, nBlock, nLayer, blockRow, nBlockCellCols)) return false;
                }
                // add a null cell row separator (except for last block): top and bottom rows with 2*nCellCols nullValues in each
                if(nBlock < nBlocks-1)
                {
                    rowComment(eclipseWriter, "ZCORN", "Block Row Separator", blockIndex, nLayer, row, eclipseCellRow);
                    eclipseWriter.printRepeatValues(fillZ, 2*nCellCols);
                    eclipseWriter.printRepeatValues(fillZ, 2*nCellCols);
                    eclipseCellRow++;
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputZCorn", " Could not write ZCORN plane at: block: " + nBlock + " layer: " + nLayer + " row: " + row, e);
            return false;
        }
    }

    private float getFirstFillZ(float[][] layerGrid, int nRows, int nCols)
    {
        for (int i = nRows - 1; i >= 0; i--)
            for (int j = 0; j < nCols; j++)
                if(layerGrid[i][j] != nullValue) return layerGrid[i][j];
        StsException.systemError(this, "getFirstFillZ", "failed to find fillZ value...using zero.");
        return 0.0f;
    }

    private boolean writeZCornRow(EclipseWriter eclipseWriter, float[][] layerGrid, int nBlock, int nLayer, int row, int nCols)
    {
        int j = -1;
        try
        {
            for (j = 0; j < nCols; j++)
            {
                checkWriteZCorn(eclipseWriter, layerGrid, row, j);
                checkWriteZCorn(eclipseWriter, layerGrid, row, j+1);
            }
            // finish cell front row out to max number of cell columns with null values
            for (; j < nCellCols; j++)
            {
                eclipseWriter.print(fillZ);
                eclipseWriter.print(fillZ);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "writeZCornRow",
                " Could not write ZCORN at: block: " + nBlock + " layer: " + nLayer + " row: " + row + " col: " + j, e);
            return false;
        }
    }

    private void checkWriteZCorn(EclipseWriter eclipseWriter, float[][] layerGrid, int row, int col)
    {
        float z = layerGrid[row][col];
        if(z == nullValue)
            z = fillZ;
        else
        {
            fillZ = z;
            if(fillZ == 0.0f)
                StsException.systemError(this, "checkWriteZCorn", "fillZ set to zero.");
        }
        eclipseWriter.print(z);
    }

    public boolean constructNNCs(float dProgress)
    {
        if(faultGrids.size() == 0) return true;
        StsMessageFiles.logMessage("Computing NNCs...");

        int nFaultGrids = faultGrids.size();
        dProgress /= nFaultGrids;
        for(StsFaultGrid faultGrid : faultGrids)
        {
            if(!faultGrid.constructNNCs(blocks, this, status, dProgress)) return false;
            status.setProgress(progress += dProgress);
        }  
        return true;
    }

    public boolean outputNNCs(EclipseWriter eclipseWriter, float dProgress)
    {
        if(faultGrids.size() == 0) return true;
        StsMessageFiles.logMessage("Writing NNC section...");
        eclipseWriter.println("NNC");
        eclipseWriter.comment("S2S index is from top front left. Indexes are IJKB where I is row, J is col, K is layer, and B is block. Indexes start from zero.");
        eclipseWriter.comment("Eclipse index is from top back left. Indexes are IJK where I is col, J is row, and K is layer. Indexes start from one.");
        eclipseWriter.comment("S2S to Eclipse indexing transformation: each S2S block has a rowMin and colMin, and a number of rowsPerBlock.");
        eclipseWriter.comment("    S2S blocks are output from the back row to the front and a blank row is inserted between blocks.");

        outputWriter.println("Number of fault NNCs: " + nncList.size());
        for(StsNNC nnc : nncList)
        {
            try
            {
                nnc.computeTrans();
            }
            catch(Exception e)
            {
                StsException.systemError(this, "outputNNCs", "Failed to compute Fault NNC: " + nnc.detailString());
            }
            eclipseWriter.comment("Fault NNC: " + nnc.toString());
            eclipseWriter.comment("           " + nnc.detailString());
            String rightBlockString = nnc.getEclipseRightBlockIndexString();
            String leftBlockString = nnc.getEclipseLeftBlockIndexString();
            eclipseWriter.println(rightBlockString + "    " + leftBlockString + "    " + nnc.trans*tranMultiplier + "/");
        }

        int nNeighborNNCs = 0;
        for(StsBlock block : blocks)
        {
            if(block.neighborNncList == null) continue;
            Collections.sort(block.neighborNncList);
            eclipseWriter.comment("Interblock NNCs for block: " + block.getIndex());
            nNeighborNNCs += block.neighborNncList.size();
            for(StsNNC nnc : block.neighborNncList)
            {
                try
                {
                    nnc.computeTrans();
                }
                catch(Exception e)
                {
                    StsException.systemError(this, "outputNNCs", "Failed to compute neighbor NNC: " + nnc.detailString());
                }
                if(nnc.trans < 10.0f) continue;
                eclipseWriter.comment("Neighbor NNC: " + nnc.toString());
                String rightBlockString = nnc.getEclipseRightBlockIndexString();
                String leftBlockString = nnc.getEclipseLeftBlockIndexString();
                eclipseWriter.println(rightBlockString + "    " + leftBlockString + "    " + nnc.trans*tranMultiplier + "/");
            }
        }
        outputWriter.println("Number of Neighbor NNCs: " + nNeighborNNCs);
        eclipseWriter.endRecord();
        progress += dProgress;
        status.setProgress(progress);
        return true;
    }

    public int getNumberMergedCells()
    {
        return numberMergedCells;
    }

    public void setNumberMergedCells(int numberMergedCells)
    {
        this.numberMergedCells = numberMergedCells;
    }

    public float getMergeCellMinFraction()
    {
        return mergeCellMinFraction;
    }

    public void setMergeCellMinFraction(float mergeCellMinFraction)
    {
        this.mergeCellMinFraction = mergeCellMinFraction;

    }

    public float getActualMergedCellMinFraction()
    {
        return actualMergedCellMinFraction;
    }

    public void setActualMergedCellMinFraction(float actualMergedCellMinFraction)
    {
        this.actualMergedCellMinFraction = actualMergedCellMinFraction;
    }

    public int getNumberUnmergedSmallCells()
    {
        return numberUnmergedSmallCells;
    }

    public void setNumberUnmergedSmallCells(int numberUnmergedSmallCells)
    {
        this.numberUnmergedSmallCells = numberUnmergedSmallCells;
    }

    public int getNumberTruncatedCells()
    {
        return numberTruncatedCells;
    }

    public void setNumberTruncatedCells(int numberTruncatedCells)
    {
        this.numberTruncatedCells = numberTruncatedCells;
    }

    class TransmissibilityComparator implements Comparator<StsNNC>
    {
        public TransmissibilityComparator(){ }

        public int compare(StsNNC nnc1, StsNNC nnc2)
        {
            if (nnc1.trans > nnc2.trans) return -1;
            else if (nnc1.trans < nnc2.trans) return 1;
            else return 0;
        }
    }

   // this creates ZCORN output and actNumArray which is output subsequently
    public boolean outputDZ(EclipseWriter eclipseWriter, float dProgress)
    {
        int n = -1, blockRow = -1, blockCol = -1, nLayer = -1;
        int row;
        int eclipseCellRow;
        float dz;
        try
        {
            eclipseWriter.initialize(fieldsPerLine);
            StsMessageFiles.logMessage("Writing DZ section...");
            eclipseWriter.println("DZ");

            dProgress /= nCelllLayers;
            int nBlocks = blocks.length;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                eclipseCellRow = 1;
                for (n = 0; n < nBlocks; n++)
                {
                    StsBlock block = blocks[n];
                    int blockIndex = block.getIndex();
                    float[][] topLayerGrid = block.getLayerGrid(nLayer);
                    float[][] botLayerGrid = block.getLayerGrid(nLayer+1);
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    eclipseWriter.comment("block: " + blockIndex + " layer: " + nLayer + " " + " nBlockCellRows: " + nBlockCellRows + " nBlockCellCols" + nBlockCellCols);
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--, eclipseCellRow++)
                    {
                        rowComment(eclipseWriter, "DZ", blockIndex, nLayer, row , row, eclipseCellRow);

                        // for each cell in block, write top left dz value (more accurate to compute average, but we aren't using these values anyways)
                        for (blockCol = 0; blockCol < nBlockCellCols; blockCol++)
                        {
                            dz = botLayerGrid[blockRow][blockCol] - topLayerGrid[blockRow][blockCol];
                            eclipseWriter.checkPrint(dz);
                        }
                        // finish cell top row out to max number of cell columns with null values
                        for (; blockCol < nCellCols; blockCol++)
                            eclipseWriter.checkPrint(0.0f);
                    }
                    // add a null cell row separator (except for last block): top and bottom rows with 2*nCellCols nullValues in each
                    if(n < nBlocks-1)
                    {
                        rowComment(eclipseWriter, "DZ", "Block Row Separator", blockIndex, nLayer, row, eclipseCellRow);
                        eclipseCellRow++;
                    }
                    eclipseWriter.checkLineFeed();
                }
                status.setProgress(progress);
            }
            eclipseWriter.endRecord();
            if(debug) debugCheck(eclipseWriter, "outputDZ");

            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputZCorn",
                " Could not write DZ at: block: " + n + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private void debugCheck(StsEclipseOutput.EclipseWriter eclipseWriter, String methodName)
    {
        debugCheck(eclipseWriter, methodName, "", 1);
    }

    private void debugCheck(StsEclipseOutput.EclipseWriter eclipseWriter, String methodName, String comment)
    {
        debugCheck(eclipseWriter, methodName, comment, 1);
    }

    private void debugCheck(StsEclipseOutput.EclipseWriter eclipseWriter, String methodName, String comment, int countPerCell)
    {
        int correctCount = nCelllLayers *nCorrectEclipseRows*nCellCols*countPerCell;
        debugCheckCount(eclipseWriter, methodName, comment, correctCount);
    }

    private void debugCheckCount(EclipseWriter eclipseWriter, String methodName, String comment, int totalCount)
    {
        if(totalCount != eclipseWriter.totalValueCount)
            StsException.systemError(this, methodName, comment + " Correct value count should be " + totalCount + " but is " + eclipseWriter.totalValueCount);
    }

    private void rowComment(EclipseWriter eclipseWriter, String dataGroup, int blockIndex, int nLayer, int blockRow, int row, int eclipseCellRow)
    {
        eclipseWriter.checkLineFeed();
        eclipseWriter.comment(dataGroup + " block: " + blockIndex + " block row (J-1): " + blockRow + " row: " + row + " layer (K-1): " + nLayer + " Eclipse row " + eclipseCellRow +  " Eclipse K " + (nLayer + 1) );
    }

    private void rowComment(EclipseWriter eclipseWriter, String dataGroup, String comment, int blockIndex, int nLayer, int blockRow, int row, int eclipseCellRow)
    {
        eclipseWriter.checkLineFeed();
        eclipseWriter.comment(dataGroup + " " + comment + " block: " + blockIndex + " layer: " + nLayer + " block row: " + blockRow + " row: " + row + " Eclipse J " + eclipseCellRow);
    }

    private void rowComment(EclipseWriter eclipseWriter, String dataGroup, String comment, int blockIndex, int nLayerGrid, int row, int eclipseCellRow)
    {
        eclipseWriter.checkLineFeed();
        eclipseWriter.comment(dataGroup + " " + comment + " block: " + blockIndex + " layer: " + nLayerGrid + " row: " + row + " Eclipse row " + eclipseCellRow);
    }
    public boolean outputProperty(EclipseWriter eclipseWriter, StsPropertyType eclipsePropertyType, float dProgress)
    {
        return outputProperty(eclipseWriter, eclipsePropertyType, dProgress, 1.0f);
    }

    public boolean outputProperty(EclipseWriter eclipseWriter, StsPropertyType propertyType, float dProgress, float multiplier)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1;
        StsPropertyVolume propertyVolume;
        int row;
        int col;
        int eclipseCellRow;
        try
        {
            StsMessageFiles.logMessage("Writing " + propertyType.eclipseName + " section...");
            eclipseWriter.println(propertyType.eclipseName);
            eclipseWriter.initialize(fieldsPerLine);
            int nBlocks = blocks.length;
            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                eclipseCellRow = 1;
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    propertyVolume = eclipseModel.getBlockPropertyVolume(block, propertyType);
                    if(propertyVolume == null)
                    {
                        StsException.systemError(this, "outputProperty", "failed to find propertyVolume for propertyVolume " + propertyType.eclipseName + " block " + block.getName());
                        return false;
                    }
                    int blockIndex = block.getIndex();
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    int colMin = block.getBlockColMin();
                    eclipseWriter.comment("block: " + blockIndex + " layer: " + nLayer + " " + " nBlockCellRows: " + nBlockCellRows + " nBlockCellCols" + nBlockCellCols);
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--)
                    {
                        rowComment(eclipseWriter, propertyType.eclipseName, blockIndex, nLayer, blockRow, row, eclipseCellRow);

                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float value = propertyVolume.getBlockValue(blockRow, blockCol, nLayer);
                            if(value != nullValue)
                                eclipseWriter.checkPrint(value*multiplier);
                            else
                                eclipseWriter.checkPrint(0.0f);
                        }
                        // finish out block cell columns beyond block with null values
                        for(; blockCol < nCellCols; blockCol++)
                            eclipseWriter.print(0);
                        eclipseCellRow++;
                    }
                    if(b < nBlocks-1)
                    {
                        rowComment(eclipseWriter, propertyType.eclipseName, "Block Row Separator", blockIndex, nLayer, row, eclipseCellRow);
                        eclipseWriter.printRepeatValues(0, nCellCols);
                        eclipseCellRow++;
                    }
                    eclipseWriter.checkLineFeed();
                }
                status.setProgress(progress);
            }
            eclipseWriter.endRecord();
            if(debug) debugCheck(eclipseWriter, "outputProperty", "Property " + propertyType.eclipseName);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputProperty",
                " Could not write " + propertyType + " for block: " + b + " at layer: " + nLayer + " blockRow: " + blockRow + " blockCol: " + blockCol, e);
            return false;
        }
    }

    public boolean outputBooleanProperty(EclipseWriter eclipseWriter, StsPropertyType eclipsePropertyType, float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        StsPropertyVolume propertyVolume;
        int row;
        int col;
        int eclipseCellRow;
        try
        {
            StsMessageFiles.logMessage("Writing " + eclipsePropertyType.eclipseName + " section...");
            eclipseWriter.println(eclipsePropertyType.eclipseName);
            eclipseWriter.initialize(fieldsPerLine);
            int nBlocks = blocks.length;
            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                eclipseCellRow = 1;
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    propertyVolume = eclipseModel.getBlockPropertyVolume(block, eclipsePropertyType);
                    if(propertyVolume == null)
                    {
                        StsException.systemError(this, "outputProperty", "failed to find propertyVolume for property " + eclipsePropertyType.eclipseName + " block " + block.getName());
                        return false;
                    }
                    int blockIndex = block.getIndex();
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    int colMin = block.getBlockColMin();
                    eclipseWriter.comment("block: " + blockIndex + " layer: " + nLayer + " " + " nBlockCellRows: " + nBlockCellRows + " nBlockCellCols" + nBlockCellCols);
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--)
                    {
                        rowComment(eclipseWriter, eclipsePropertyType.eclipseName, blockIndex, nLayer, blockRow, row, eclipseCellRow);

                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float value = propertyVolume.getBlockValue(blockRow, blockCol, nLayer);
                            if(value == 0.0f || value == nullValue)
                                eclipseWriter.print(0);
                            else
                                eclipseWriter.print(1);
                        }
                        // finish out block cell columns beyond block with null values
                        for(; blockCol < nCellCols; blockCol++)
                            eclipseWriter.print(0);
                        eclipseCellRow++;
                    }
                    if(b < nBlocks-1)
                    {
                        rowComment(eclipseWriter, eclipsePropertyType.eclipseName, "Block Row Separator", blockIndex, nLayer, row, eclipseCellRow);
                        eclipseWriter.printRepeatValues(0, nCellCols);
                        eclipseCellRow++;
                    }
                    eclipseWriter.checkLineFeed();
                }
                status.setProgress(progress);
            }
            eclipseWriter.endRecord();
            if(debug) debugCheck(eclipseWriter, "outputProperty", "Property " + eclipsePropertyType.eclipseName);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputBooleanProperty",
                " Could not write " + eclipsePropertyType + " for block: " + b + " at layer: " + nLayer + " blockRow: " + blockRow + " blockCol: " + blockCol, e);
            return false;
        }
    }

    public boolean outputConstantProperty(EclipseWriter eclipseWriter, String eclipsePropertyName, String comment, double value, int nValues)
    {
        StsMessageFiles.logMessage("Writing " + eclipsePropertyName + " section...");
        eclipseWriter.println(eclipsePropertyName);
        eclipseWriter.comment(eclipsePropertyName + " constant values " + comment);
        eclipseWriter.println(nValues + "*" + value);
        eclipseWriter.endRecord();
        return true;
    }

    private boolean computePoreVolume(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block PORV propertyVolumes...");
            int nBlocks = blocks.length;
            int[] nBlockActiveCells = new int[nBlocks];
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.poreVolume, StsPropertyVolume.PROP_VARIABLE_XYZ);
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.actnum, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }

            dProgress /= nCelllLayers;
            int nActiveCell = 0;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockPoreVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.poreVolume);
                    StsPropertyVolume blockActnumVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.actnum);

                    StsPropertyVolume porosityVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.porosity);
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMin = block.getBlockRowMin();
                    int colMin = block.getBlockColMin();

                    for (blockRow = 0, row = rowMin; blockRow < nBlockCellRows; blockRow++, row++)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            StsBlock.BlockCellColumn.GridCell gridCell = block.getBlockGridCell(blockRow, blockCol, nLayer);
                            if(gridCell == null) continue;

                            float poreVolume = gridCell.poreVolume;
                            if(poreVolume <= 0.0f)
                                block.deleteBlockGridCell(blockRow, blockCol, nLayer);
                            else
                            {
                                blockActnumVolume.setBlockValue(blockRow, blockCol, nLayer, 1.0f);
                                blockPoreVolume.setBlockValue(blockRow, blockCol, nLayer, poreVolume);
                                nBlockActiveCells[b]++;
                                nActiveCell++;
                            }
                        }
                    }
                }
                status.setProgress(progress);
            }
            StsMessageFiles.logMessage("Number of cells with non-zero pore volume: " + nActiveCell);
            outputWriter.println("Summary of active cells");
            outputWriter.println("Number of cells with non-zero pore volume: " + nActiveCell);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputTranX",
                " Could not write TRANX at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }
   private boolean computeDepths(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block depth propertyVolumes...");
            int nBlocks = blocks.length;
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.depth, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }

            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockDepthVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.depth);

                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMin = block.getBlockRowMin();
                    int colMin = block.getBlockColMin();
                    float[][] topLayerGrid = block.getLayerGrid(nLayer);
                    float[][] botLayerGrid = block.getLayerGrid(nLayer+1);
                    float topFillZ = getFirstFillZ(topLayerGrid, nBlockCellRows, nBlockCellCols);
                    float botFillZ = getFirstFillZ(botLayerGrid, nBlockCellRows, nBlockCellCols);
                    fillZ = (topFillZ + botFillZ)/2;
                    for (blockRow = 0, row = rowMin; blockRow < nBlockCellRows; blockRow++, row++)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float depth;
                            StsBlock.BlockCellColumn.GridCell gridCell = block.getBlockGridCell(blockRow, blockCol, nLayer);
                            if(gridCell == null)
                                depth = fillZ;
                            else
                            {
                                depth = gridCell.getCellDepth();
                                fillZ = depth;
                            }
                            blockDepthVolume.setBlockValue(blockRow, blockCol, nLayer, depth);
                        }
                    }
                }
                status.setProgress(progress);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeDepths",
                " Could not compute depth at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private boolean outputDepths(EclipseWriter eclipseWriter, float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        int eclipseCellRow;
        try
        {
            eclipseWriter.initialize(fieldsPerLine);
            StsMessageFiles.logMessage("Writing DEPTH section...");
            eclipseWriter.println("DEPTH");

            int nBlocks = blocks.length;

            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.depth, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }


            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                eclipseCellRow = 1;
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume depthVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.depth);
                    int blockIndex = block.getIndex();
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    int colMin = block.getBlockColMin();
                    float[][] topLayerGrid = block.getLayerGrid(nLayer);
                    float[][] botLayerGrid = block.getLayerGrid(nLayer+1);
                    float topFillZ = getFirstFillZ(topLayerGrid, nBlockCellRows, nBlockCellCols);
                    float botFillZ = getFirstFillZ(botLayerGrid, nBlockCellRows, nBlockCellCols);
                    fillZ = (topFillZ + botFillZ)/2;
                    eclipseWriter.comment("block: " + blockIndex + " layer: " + nLayer + " " + " nBlockCellRows: " + nBlockCellRows + " nBlockCellCols" + nBlockCellCols);
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--)
                    {
                        rowComment(eclipseWriter, "DEPTH", blockIndex, nLayer, blockRow, row, eclipseCellRow);

                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float cellDepth = block.computeCellDepth(blockRow, blockCol, nLayer, topLayerGrid, botLayerGrid, fillZ);
                            fillZ = cellDepth;
                            eclipseWriter.checkPrint(cellDepth);
                            depthVolume.setBlockValue(blockRow, blockCol, nLayer, cellDepth);
                        }
                        // finish out block cell columns beyond block with null values
                        for(; blockCol < nCellCols; blockCol++)
                            eclipseWriter.print(0.0);
                        eclipseCellRow++;
                    }
                    if(b < nBlocks-1)
                    {
                        rowComment(eclipseWriter, "DEPTH", "Block Row Separator", blockIndex, nLayer, row, eclipseCellRow);
                        eclipseWriter.printRepeatValues(0, nCellCols);
                        eclipseCellRow++;
                    }
                    eclipseWriter.checkLineFeed();
                }
                status.setProgress(progress);
            }
            eclipseWriter.endRecord();
            if(debug) debugCheck(eclipseWriter, "outputDepths");
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputDepths",
                " Could not write DEPTH at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private boolean computeBlockTranX(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            dataWriter.initialize(fieldsPerLine);
            StsMessageFiles.logMessage("Creating Block TRANX propertyVolumes...");

            int nBlocks = blocks.length;
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.tranX, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }

            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockVolumeTranx = eclipseModel.getBlockPropertyVolume(block, eclipseModel.tranX);
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMin = block.getBlockRowMin();
                    int colMin = block.getBlockColMin();
                    for (blockRow = 0, row = rowMin; blockRow < nBlockCellRows; blockRow++, row++)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols-1; blockCol++, col++)
                        {
                            if(StsBlock.debugIJK(row, col, nLayer))
                            {
                                StsException.systemDebug(this, "computeBlockTranX");
                            }
                            float trans = block.computeTranXPlus(blockRow, blockCol, nLayer);
                            blockVolumeTranx.setBlockValue(blockRow, blockCol, nLayer, trans);
                        }
                    }
                }
                status.setProgress(progress);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "constructBlockTranX",
                " Could not write TRANX at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private boolean computeBlockTranY(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block TRANY propertyVolumes...");

            int nBlocks = blocks.length;
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.tranY, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }

            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockVolumeTrany = eclipseModel.getBlockPropertyVolume(block, eclipseModel.tranY);
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMin = block.getBlockRowMin();
                    int colMin = block.getBlockColMin();
                    for (blockRow = 0, row = rowMin; blockRow < nBlockCellRows-1; blockRow++, row++)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float trans = block.computeTranYPlus(blockRow, blockCol, nLayer);
                            blockVolumeTrany.setBlockValue(blockRow, blockCol, nLayer, trans);
                        }
                    }
                }
                status.setProgress(progress);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputTranX",
                " Could not write TRANX at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private boolean computeBlockTranZ(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block TRANZ propertyVolumes...");

            int nBlocks = blocks.length;
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                eclipseModel.createBlockPropertyVolume(block, eclipseModel.tranZ, StsPropertyVolume.PROP_VARIABLE_XYZ);
            }
            dProgress /= nCelllLayers;
            for (nLayer = 0; nLayer < nCelllLayers -1; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockVolumeTranZ = eclipseModel.getBlockPropertyVolume(block, eclipseModel.tranZ);
                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMin = block.getBlockRowMin();
                    int colMin = block.getBlockColMin();
                    for (blockRow = 0, row = rowMin; blockRow < nBlockCellRows-1; blockRow++, row++)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            float trans = block.computeTranZPlus(blockRow, blockCol, nLayer);
                            blockVolumeTranZ.setBlockValue(blockRow, blockCol, nLayer, trans);
                        }
                    }
                }
                status.setProgress(progress);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputTranX",
                " Could not write TRANX at: block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }
/*
    private void mergeCells(float dProgress)
    {
        // if(minEdgeCellTimeConstant > 0.25f*minCellTimeConstant) return;
        MergeCellPanel mergeCellPanel = new MergeCellPanel();
                    new StsOkApplyCancelDialog(model.win3d, mergeCellPanel, "Merge Cell Operation", true);
        status.setProgress(progress += dProgress);
    }
*/
    private void createTruncatedCellListComputeFractions()
    {
        for(StsBlock block : blocks)
            block.createTruncatedCellListComputeTruncatedFraction(mergeCellMinFraction);

        computeTruncatedCellFractions();
    }

    private void computeTruncatedCellFractions()
    {
        actualMergedCellMinFraction = blocks[0].actualMergedCellMinFraction;
        numberTruncatedCells = blocks[0].getNTruncatedCells();
        numberMergedCells = blocks[0].numberMergedCells;
        numberUnmergedSmallCells = blocks[0].numberUnmergedSmallCells;
        int nBlocks = blocks.length;
        for(int n = 1; n < nBlocks; n++)
        {
            actualMergedCellMinFraction = Math.min(actualMergedCellMinFraction, blocks[n].actualMergedCellMinFraction);
            numberTruncatedCells += blocks[n].getNTruncatedCells();
            numberMergedCells += blocks[n].numberMergedCells;
            numberUnmergedSmallCells += blocks[n].numberUnmergedSmallCells;
        }
    }

    private void outputTruncationInfo()
    {
        outputWriter.println("Cell truncation summary: ");
        outputWriter.println("    mergeCellMinFraction " + mergeCellMinFraction);
        outputWriter.println("    numberMergedCells " + numberMergedCells);
        outputWriter.println("    numberTruncatedCells " + numberTruncatedCells);
        outputWriter.println("    actualMergedCellMinFraction " + actualMergedCellMinFraction);
        outputWriter.println("    numberUnmergedSmallCells " + numberUnmergedSmallCells);
    }
    
    private void merge()
    {
        constructNNCs();
        combineNNCs();
        debugCellList = new TreeSet<StsBlock.BlockCellColumn.GridCell>(StsBlock.getTruncatedFractionComparator());
        for(StsBlock block : blocks)
            block.merge(debugCellList, eclipseModel, this);
    }

    private void constructNNCs()
    {
        nncList = new ArrayList<StsNNC>();
        for(StsFaultGrid faultGrid : faultGrids)
        {
            for(StsNNC nnc : faultGrid.nncList)
            {
                try
                {
                    nnc.computeTrans();
                }
                catch(Exception e)
                {
                    StsException.systemError(this, "constructNNCs", "Failed to compute NNC: " + nnc.detailString());
                }
                StsNNC addNnc;
                if(nnc.parentNNC == null)
                    addNnc = nnc;
                else
                    addNnc = nnc.parentNNC;

                // if(addNnc.trans >= 0.001f)
                addNnc.adjustArea();
                nncList.add(addNnc);
            }
        }
    }

    private void combineNNCs()
    {
        nncList = new ArrayList<StsNNC>();
        for(StsFaultGrid faultGrid : faultGrids)
        {
            for(StsNNC nnc : faultGrid.nncList)
            {
                StsNNC addNnc;
                if(nnc.parentNNC == null)
                    addNnc = nnc;
                else
                    addNnc = nnc.parentNNC;

                // if(addNnc.trans >= 0.001f)
                addNnc.adjustArea();
                nncList.add(addNnc);
            }
        }
        // sort NNCs by lowest (between left and right) IJKB order in reverse sequence (i.e., B==0 first, B==1, next, largest I last
        try
        {
            Collections.sort(nncList);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "combineNNCs", "Sort failed", e);
        }
        Iterator<StsNNC> iterator = nncList.iterator();
        StsNNC combineNNC = null;
        while (iterator.hasNext())
        {
            StsNNC nnc = iterator.next();
            // if NNCs have the same connection, sum the transmissiblity
            // and remove this redundant connection
            if (combineNNC != null)
            {
                if (nnc.sameNNC(combineNNC))
                {
                    combineNNC.trans += nnc.trans;
                    combineNNC.area += nnc.area;
                    iterator.remove();
                }
                else
                    combineNNC = nnc;
            }
            else
                combineNNC = nnc;
        }
    }

    private boolean computeActiveCellIndexMaps(float dProgress)
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block PORV propertyVolumes...");
            int nBlocks = blocks.length;
            int[] nBlockActiveCells = new int[nBlocks];
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                StsPropertyVolume blockIndexVolume = eclipseModel.createBlockPropertyVolume(block, eclipseModel.indexMap, StsPropertyVolume.PROP_VARIABLE_XYZ);
                blockIndexVolume.initializeValues(-1.0f);
            }

            dProgress /= nCelllLayers;
            int nActiveCell = 0;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++, progress += dProgress)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockActnumVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.actnum);
                    StsPropertyVolume blockIndexVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.indexMap);

                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    int colMin = block.getBlockColMin();
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            boolean isActive = blockActnumVolume.getValue(row, col, nLayer) == 1.0f;
                            if(isActive)
                            {
                                blockIndexVolume.setBlockValue(blockRow, blockCol, nLayer, nActiveCell);
                                StsBlock.BlockCellColumn.GridCell gridCell = block.getBlockGridCell(blockRow, blockCol, nLayer);
                                if(!adjustCells && gridCell.parentCell != null)
                                {
                                    StsBlock.BlockCellColumn.ParentCell parentCell = gridCell.parentCell;
                                    ArrayList<StsBlock.BlockCellColumn.GridCell> children = parentCell.children;
                                    for(StsBlock.BlockCellColumn.GridCell child : children)
                                    {
                                        int[] blockIJK = child.getBlockIJK();
                                        blockIndexVolume.setBlockValue(blockIJK, nActiveCell);
                                    }
                                }
                                nActiveCell++;
                            }
                        }
                    }
                }
                status.setProgress(progress);
            }
            StsMessageFiles.logMessage("Number of active cells after merging: " + nActiveCell);
            outputWriter.println("Number of active cells after merging: " + nActiveCell);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeActiveCellIndexMaps",
                " block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }
     private boolean countActiveCells()
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Creating Block PORV propertyVolumes...");
            int nBlocks = blocks.length;

            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                StsPropertyVolume blockIndexVolume = eclipseModel.createBlockPropertyVolume(block, eclipseModel.indexMap, StsPropertyVolume.PROP_VARIABLE_XYZ);
                blockIndexVolume.initializeValues(-1.0f);
            }

            int nActiveCells = 0;
            int nPoreVolumeCells = 0;
            for (nLayer = 0; nLayer < nCelllLayers; nLayer++)
            {
                for (b = 0; b < nBlocks; b++)
                {
                    StsBlock block = blocks[b];
                    StsPropertyVolume blockActnumVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.actnum);
                    StsPropertyVolume blockPoreVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.poreVolume);

                    int nBlockCellRows = block.nCellRows;
                    int nBlockCellCols = block.nCellCols;
                    int rowMax = block.getBlockRowMax();
                    int colMin = block.getBlockColMin();
                    for (blockRow = nBlockCellRows - 1, row = rowMax-1; blockRow >= 0; blockRow--, row--)
                    {
                        for (blockCol = 0, col = colMin; blockCol < nBlockCellCols; blockCol++, col++)
                        {
                            boolean isActive = blockActnumVolume.getValue(row, col, nLayer) == 1.0f;
                            if(isActive)nActiveCells++;
                            float poreVolume = blockPoreVolume.getValue(row, col, nLayer);
                            if(poreVolume > 0) nPoreVolumeCells++;
                        }
                    }
                }
                status.setProgress(progress);
            }
            StsException.systemDebug(this, "countActiveCells","Number of active cells after merging: " + nActiveCells);
            StsException.systemDebug(this, "countActiveCells","Number of non-zero pore volumes: " + nPoreVolumeCells);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeActiveCellIndexMaps",
                " block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
            return false;
        }
    }

    private void removeSmallPoreVolumes()
    {
        int b = -1, blockRow = -1, blockCol = -1, nLayer = -1, z;
        int row;
        int col;
        try
        {
            StsMessageFiles.logMessage("Hack removing small pore volumes...");
            int nBlocks = blocks.length;
            int[] nBlockActiveCells = new int[nBlocks];

            int nActiveCell = 0;
            int nRemoved = 0;
            for (b = 0; b < nBlocks; b++)
            {
                StsBlock block = blocks[b];
                StsPropertyVolume blockActnumVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.actnum);
                StsPropertyVolume blockPoreVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.poreVolume);
                float[] actNums = blockActnumVolume.getValues();
                float[] poreVolumes = blockPoreVolume.getValues();
                int nValues = poreVolumes.length;
                for(int n = 0; n < nValues; n++)
                {
                    if(actNums[n] == 0.0f) continue;

                    if(poreVolumes[n] < minPoreVolume)
                    {
                        poreVolumes[n] = 0.0f;
                        actNums[n] = 0.0f;
                        nRemoved++;
                    }
                    else
                    {
                        nActiveCell++;
                        nBlockActiveCells[b]++;
                    }
                }
            }
            StsMessageFiles.logMessage("Number of small pore volumes removed: " + nRemoved);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "hackRemoveSmallPoreVolumes",
                " block: " + b + " layer: " + nLayer + " row: " + blockRow + " col: " + blockCol, e);
        }
    }
/*
    private void createIndexMapVolumes()
    {
        int[] activeIndex = new int[nTotalCells];

        for(StsBlock block : blocks)
        {
            StsPropertyVolume indexVolume = eclipseModel.createBlockPropertyVolume(block, StsPropertyVolume.indexMap, StsPropertyVolume.PROP_VARIABLE_XYZ, null);
            block.createIndexMap(indexVolume, nCellRows);
        }

        for(StsBlock block : blocks)
        {
            StsPropertyVolume indexVolume = eclipseModel.createBlockPropertyVolume(block, StsPropertyVolume.indexMap, StsPropertyVolume.PROP_VARIABLE_XYZ, null);
            block.createIndexMap(indexVolume, nCellRows);
        }
    }
*/
    static public StsJPanel getEclipseMergeCellPanel(StsModel model)
    {
        if(eclipseOutputInstance == null) StsEclipseOutput.getInstance(model);
        return eclipseOutputInstance.getMergeCellPanel();
    }

    protected StsJPanel getMergeCellPanel()
    {
        return new MergeCellPanel();
    }

    class MergeCellPanel extends StsJPanel
    {
        StsGroupBox mergeCellsBox = new StsGroupBox("Truncated Cell Merge");

        MergeCellPanel()
        {
            mergeCellsBox.add(new StsFloatFieldBean(StsEclipseOutput.this, "mergeCellMinFraction", "Target merge cell min fraction"));
            mergeCellsBox.add(new StsIntFieldBean(StsEclipseOutput.this, "numberTruncatedCells", false, "Number of truncated cells"));
            mergeCellsBox.add(new StsIntFieldBean(StsEclipseOutput.this, "numberMergedCells", false, "Number of merged cells"));
            mergeCellsBox.add(new StsIntFieldBean(StsEclipseOutput.this, "numberUnmergedSmallCells", false, "Number of small cells"));
            mergeCellsBox.add(new StsFloatFieldBean(StsEclipseOutput.this, "actualMergedCellMinFraction", false, "Actual Merge cell min fraction"));
            mergeCellsBox.add(new StsButton("Apply", "Apply truncation size", this, "apply"));
            add(mergeCellsBox);
        }

        public void apply()
        {
            checkClearParentCells();
            checkMergeCells();
            if(!adjustCells) constructParentNNCs();
            recomputeTruncatedCellFractions();
            updateBeans();
        }
    }

    private void checkClearParentCells()
    {
        for(StsBlock block : blocks)
            block.checkClearParentCells();
    }

    public void checkMergeCells()
    {
        for(StsBlock block : blocks)
            block.checkMergeCells(this);
    }

    public void constructParentNNCs()
    {
    if(faultGrids.size() == 0) return;
    StsMessageFiles.logMessage("Computing Parent NNCs...");
        for(StsFaultGrid faultGrid : faultGrids)
        {
            ArrayList<StsNNC> nncList = faultGrid.nncList;
            for(StsNNC nnc : nncList)
                nnc.constructParentNNC();
        }
    }

    private void recomputeTruncatedCellFractions()
    {
        for(StsBlock block : blocks)
        {
            if(!adjustCells)
                block.recomputeParentCellTruncatedCellFractions(mergeCellMinFraction);
            else
                block.recomputeAdjustedTruncatedCellFractions(mergeCellMinFraction);
        }
        computeTruncatedCellFractions();
    }

    public boolean outputProperty(EclipseWriter eclipseWriter, String label, StsPropertyVolumeOld property,
                                  double multiplier, int fractionDigits)
    {
        int i = -1, j = -1, k = -1;

        try
        {
            if (property == null)
            {
                StsMessageFiles.logMessage("No property volume available for: " + label);
                return true;
            }

            eclipseWriter.initialize(fieldsPerLine);

            float[] values = property.getValues();

            StsMessageFiles.logMessage("Writing " + label + " section...");
            eclipseWriter.println(label);
            int cellRowMin = 0;
            int cellColMin = 0;
            for (k = 0; k < nCelllLayers; k++)
            {
                for (i = nCellRows; i >= cellRowMin; i--)
                {
                    eclipseWriter.comment("layer: " + k + " cell row: " + i + " cols: " + cellColMin + " - " + nCellCols);
                    for (j = cellColMin; j < nCellCols; j++)
                    {
                        int index = k * nCellRows * nCellCols + i * nCellCols + j;
                        eclipseWriter.print(values[index] * multiplier);
                    }
                }
            }
            eclipseWriter.endRecord();
            return true;
        }
        catch (Exception e)
        {

            StsException.outputWarningException(this, "outputProperty",
                " Could not write Property " + label + " at layer: " + k + " row: " + i + " col: " + j, e);
            return false;
        }
    }

    private boolean outputSchedule(EclipseWriter eclipseWriter)
    {
        StsWell[] wells = (StsWell[])model.getCastObjectList(StsWell.class);
        eclipseWriter.initialize(20);
        StsMessageFiles.logMessage("Writing Well completion section...");

        eclipseWriter.println("WELSPECS");
        for(StsWell well : wells)
            if(!outputWellSpec(well, eclipseWriter)) return false;
        eclipseWriter.println("/");

        eclipseWriter.println("COMPDAT");
        for(StsWell well : wells)
            if(!outputCompDat(well, eclipseWriter)) return false;
        eclipseWriter.println("/");

        return true;

    }

    private boolean outputWellSpec(StsWell well, EclipseWriter eclipseWriter)
    {
        float[] xyz = well.getTopRotatedPoint().getXYZ();
        float rowF = modelBoundingBox.getRowCoor(xyz[1]);
        float colF = modelBoundingBox.getColCoor(xyz[0]);
        int[] rowCol = getEclipseRowColF(rowF, colF);
        eclipseWriter.comment("well " + well.getName() + " S2S: rowF " + rowF + " colF " + colF);
        eclipseWriter.println(well.getName() + " G1 " + rowCol[1] + " " + rowCol[0] + " 1* " + " Separator " + " /");
        return true;
    }

    int[] getEclipseRowCol(int row, int col)
    {
        for(StsBlock block : blocks)
        {
            int[] rowCol = block.getEclipseRowCol(row, col);
            if(rowCol != null) return rowCol;
        }
        return null;
    }

    int[] getEclipseRowColF(float rowF, float colF)
    {
        for(StsBlock block : blocks)
        {
            int[] rowCol = block.getEclipseRowColF(rowF, colF);
            if(rowCol != null) return rowCol;
        }
        return null;
    }
    /** for preproject, wells are completed through all layers */
    private boolean outputCompDat(StsWell well, EclipseWriter eclipseWriter)
    {
        eclipseWriter.println(well.getName() + " 2*" + " 1 " + nCelllLayers + " 3* " + " 0.2 " + " /");
        return true;
    }

    public NumberFormat getNumberFormat(){ return numberFormat; }

    public void setValuesPerLine(EclipseWriter eclipseWriter, int valuesPerLine)
    {
        eclipseWriter.valuesPerLine = valuesPerLine;
    }

    public StsBlock[] getBlocks()
    {
        return blocks;
    }

    static final int maxCharactersPerLine = 130;

    public class EclipseWriter extends PrintWriter
    {
        StsDecimalFormat numberFormat;
        int valuesPerLine = 10;
        int lineValueCount = 0;
        int totalValueCount = 0;
        int lineCharacterCount = 0;

        EclipseWriter(String pathname, StsDecimalFormat numberFormat, int valuesPerLine) throws IOException
        {
            super(new FileWriter(pathname, false));
            this.numberFormat = numberFormat;
            this.valuesPerLine = valuesPerLine;
        }

        EclipseWriter(String pathname, StsDecimalFormat numberFormat) throws IOException
        {
            super(new FileWriter(pathname, false));
            this.numberFormat = numberFormat;
        }

        public void initialize(int valuesPerLine)
        {
            this.valuesPerLine = valuesPerLine;
            lineValueCount = 0;
            totalValueCount = 0;
            lineCharacterCount = 0;
        }

        public void comment(String string)
        {
            checkLineFeed();
            println("-- " + string);
        }

        public void endRecord()
        {
            checkLineFeed();
            println("/");
        }

        public void print(double value)
        {
            String valueString = numberFormat.formatValue(value) + "  ";
            checkLineLength(valueString);
            super.print(valueString);
            incrementLineCount();
            totalValueCount++;
        }

        public void checkPrint(float value)
        {
            if(value == StsParameters.nullValue) value = 0;
            print(value);
        }

        public void checkPrintRepeatValues(float value, int nRepeats)
        {
            if(value == StsParameters.nullValue) value = 0;
            printRepeatValues(value, nRepeats);
        }

        public void print(float value)
        {
            String valueString = numberFormat.formatValue(value) + "  ";
            checkLineLength(valueString);
            // print(valueString);
            try
            {
                out.write(valueString);
            }
            catch(Exception e)
            {
                StsException.outputWarningException(this, "print", "value: " + value, e);
            }
            incrementLineCount();
            totalValueCount++;
        }

        public void printRepeatValues(float value, int nRepeats)
        {
            checkLineFeed();
            String repeatValueString = nRepeats + "*" + numberFormat.formatValue(value);
            checkLineLength(repeatValueString);
            println(repeatValueString);
            totalValueCount += nRepeats;
        }

        public void printRepeatValues(int value, int nRepeats)
        {
            checkLineFeed();
            String repeatValueString = nRepeats + "*" + value;
            checkLineLength(repeatValueString);
            println(repeatValueString);
            totalValueCount += nRepeats;
        }

        public void print(int value)
        {
            String valueString = value + "  ";
            checkLineLength(valueString);
            print(valueString);
            incrementLineCount();
            totalValueCount++;
        }

        void println(double[] values)
        {
            print(values);
        }

        void println(float[] values)
        {
            print(values);
        }

        void println(int[] values)
        {
            print(values);
        }

        void print(double[] values)
        {
            if (values == null) return;
            int nValues = values.length;
            for (int n = 0; n < nValues; n++)
            {
                String valueString = numberFormat.formatValue(values[n]) + "  ";
                checkLineLength(valueString);
                print(valueString);
                incrementLineCount();
            }
            totalValueCount += nValues;
        }

        void print(float[] values)
        {
            if (values == null) return;
            int nValues = values.length;
            for (int n = 0; n < nValues; n++)
            {
                String valueString = numberFormat.formatValue(values[n]) + "  ";
                checkLineLength(valueString);
                print(valueString);
                incrementLineCount();
            }
            totalValueCount += nValues;
        }

        void print(int[] values)
        {
            if (values == null) return;
            int nValues = values.length;
            for (int n = 0; n < nValues; n++)
            {
                String valueString = values[n] + "  ";
                checkLineLength(valueString);
                print(valueString);
                incrementLineCount();
            }
            totalValueCount += nValues;
        }

        void checkLineFeed()
        {
            if(lineValueCount != 0) lineFeed();
        }

        void incrementLineCount()
        {
            if (++lineValueCount == valuesPerLine) lineFeed();
        }

        void lineFeed()
        {
            lineValueCount = 0;
            lineCharacterCount = 0;
            println("");
        }

        private void checkLineLength(String newString)
        {
            lineCharacterCount += newString.length();
            if(lineCharacterCount > maxCharactersPerLine)
                lineFeed();
        }
    }
}

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

/** These are currently transient instances which are created when the prestack volume is constructed or loaded from existing files. */
public class StsPreStackVolume extends StsSeismicVolume
{
    // StsParameterFile which writes out this instance as an ASCII file cannot handle multidimensional arrays.
    // So we persist these, but have a convenience copy in planesOK array.
    private boolean[] rowPlanesOK = null;
    private boolean[] colPlanesOK = null;
    private boolean[] slicePlanesOK = null;

    transient private boolean[][] planesOK = null;

    transient protected StsSuperGather gather;
    transient public StsPreStackLineSet3d preStackLineSet3d;
    transient public boolean volumeOK = false;
    transient public boolean volumeChanged = false;

    static public byte[] byteTransparentTrace = null;
    static public float[] floatTransparentTrace = null;

    static final byte VELOCITY = StsPreStackLineSetClass.DISPLAY_VELOCITY;
    static final byte STACKED = StsPreStackLineSetClass.DISPLAY_STACKED;
    static final byte SEMBLANCE = StsPreStackLineSetClass.DISPLAY_SEMBLANCE;
    static public final String[] typenames = StsPreStackLineSetClass.DISPLAY_TYPE_STRINGS;

    static public final boolean debugTimer = false;

    public StsPreStackVolume()
    {
    }

    private StsPreStackVolume(StsPreStackLineSet3d preStackSeismic, byte type)
    {
        super(false);
        this.preStackLineSet3d = preStackSeismic;
        this.type = type;
        checkInitializePlaneOKFlags();
        initialize();
    }

    private StsPreStackVolume(StsPreStackLineSet3d preStackLineSet3d, String filename, StsModel model, byte type) throws FileNotFoundException, StsException
    {
        super(false);
        this.preStackLineSet3d = preStackLineSet3d;
        this.type = type;
        loadFromFile(preStackLineSet3d.stsDirectory, filename, model);
    }

    /*
       static private StsPreStackVolume initializeVolume(StsModel model, StsPreStackLineSet3d preStackSeismic, float dataMin, float dataMax, String stemname, byte type)
       {
            String attributeName = typenames[type];
            stemname = stemname + "." + attributeName;
            StsPreStackVolume attributeVolume = (StsPreStackVolume) model.getObjectWithName(StsPreStackVolume.class, stemname);
            if (attributeVolume != null) attributeVolume.delete();
            attributeVolume = new StsPreStackVolume(preStackSeismic, type);
            byte zDomain = preStackSeismic.getZDomain();
            attributeVolume.initializeVolume(model, preStackSeismic, dataMin, dataMax, true, stemname, zDomain, false);
            return attributeVolume;
       }
    */
    static public StsPreStackVolume checkLoadFromFilename(StsModel model, StsPreStackLineSet3d preStackLineSet3d, String filename, byte type)
    {
//        StsPreStackVolume volume;
        try
        {
            preStackLineSet3d.checkStsDirectoryForFilename(filename);
            return new StsPreStackVolume(preStackLineSet3d, filename, model, type);
        }
        catch(FileNotFoundException e)
        {
            new StsMessage(model.win3d, StsMessage.ERROR,  e.getMessage());
            return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackVolume.constructor() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public StsPreStackVolume checkLoadFromStemname(StsModel model, StsPreStackLineSet3d preStackLineSet3d, String stemname, byte type)
    {
        String filename = createHeaderFilename(group3d, stemname);
        return checkLoadFromFilename(model, preStackLineSet3d, filename, type);
    }

    /** for now, we will always build the velocity prestack volume; ultimately we will write it out at the end on exit.
     *  stack and semblance volumes are loaded from files if they are available
     */
    static public StsPreStackVolume loadConstruct(StsModel model, StsPreStackLineSet3d preStackLineSet3d, byte type, String volumeStemname, boolean loadFromFile)
    {
        StsPreStackVolume attributeVolume;
        String stemname = createStemname(volumeStemname, type);
        if(type == VELOCITY || !loadFromFile)
        {
            attributeVolume = construct(model, preStackLineSet3d, type, stemname);
        }
        else
        {
            attributeVolume = checkLoadFromStemname(model, preStackLineSet3d, stemname, type);
            if(attributeVolume == null) attributeVolume = construct(model, preStackLineSet3d, type, stemname);
        }
        return attributeVolume;
    }

     public Class getLoadSubClass() { return  StsPreStackVolume.class; }

     public byte getScaleType()
     {
         if(type == VELOCITY)
            return SCALE_VOLUME;
         else if(type == STACKED)
            return SCALE_TRACE_NORMALIZED;
         else // SEMBLANCE
            return SCALE_NORMALIZED;
     }

    private float[] getDataRange(byte type)
    {
        if(type == VELOCITY)
            return new float[]{preStackLineSet3d.velocityModel.dataMin, preStackLineSet3d.velocityModel.dataMax};
        else if(type == STACKED)
            // stacked data scaled by total data range
            return new float[]{preStackLineSet3d.dataMin, preStackLineSet3d.dataMax};
//            return new float[]{preStackLineSet3d.dataMin, preStackLineSet3d.dataMax};
        else // SEMBLANCE
            return new float[]{0.0f, 1.0f};
    }

    protected void loadFromFile(String pathname, StsModel model) throws FileNotFoundException, StsException
    {
        StsParameterFile.initialReadObjectFields(pathname, this, StsPreStackVolume.class, StsBoundingBox.class);
        setName(stemname);
        initialize();
        isVisible = true;
        getSeismicVolumeClass().setIsVisibleOnCursor(true);
    }

    static public StsPreStackVolume construct(StsModel model, StsPreStackLineSet3d preStackSeismic, byte type, String stemname)
    {
//         StsPreStackVolume preStackVolume = (StsPreStackVolume) model.getObjectWithName(StsPreStackVolume.class, stemname);
//         if (preStackVolume != null) preStackVolume.delete();
        StsPreStackVolume preStackVolume = new StsPreStackVolume(preStackSeismic, type);
        byte zDomain = preStackSeismic.getZDomain();
        float[] dataRange = preStackVolume.getDataRange(type);
        preStackVolume.initializeVolume(model, preStackSeismic, dataRange[0], dataRange[1], true, true, stemname, zDomain, false, "rw");
        return preStackVolume;
    }

    public boolean checkComputePreStackVolume()
    {
        if(volumeOK) return true;
        try
        {
            StsPreStackVolumeConstructor.computeVolumeWithDialog(currentModel, preStackLineSet3d, type, this);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "checkComputePreStackVolume", "Failed to compute.", e);
            setAllPlaneOKFlagsFalse();
            return false;
        }
    }

    public boolean writeHeaderFile()
    {
        try
        {
            String filename = getCreateHeaderFilename(stemname);
//            outputTime = System.currentTimeMillis();
            StsParameterFile.writeObjectFields(stsDirectory + filename, this, StsPreStackVolume.class, StsBoundingBox.class);
//            preStackLineSet3d.setOutputVolumeTime(type, outputTime);
            if(debug) System.out.println("Writing header file: " + stsDirectory + filename);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "writeHeaderFile", e);
            return false;
        }
    }

    public boolean closeChanged()
    {
        if(volumeChanged) writeHeaderFile();
        return volumeChanged;
    }

    public boolean setGLColorList(GL gl, boolean nullsFilled, int dir, int shader)
    {
        if(dir == ZDIR && currentAttribute != null && currentAttribute != nullAttribute)
            return currentAttribute.setGLColorList(gl, nullsFilled, shader);
        else
            return preStackLineSet3d.setGLColorList(gl, nullsFilled, type, shader);
    }

    String getTypeName() { return typenames[type]; }

    protected ByteBuffer getComputeColPlaneByteBuffer(int nPlane, boolean neighborsOnly)
    {
        StsAGCPreStackProperties agcProperties = preStackLineSet3d.agcProperties;
        boolean applyAGCPoststack = agcProperties.getApplyAGCPoststack();
        int windowWidth = 0;
        float maxAmplitude = 0.0f;
        if(applyAGCPoststack)
        {
            windowWidth = agcProperties.getWindowWidth(StsFilterProperties.POSTSTACK);
            maxAmplitude = Math.max(Math.abs(dataMin), Math.abs(dataMax));
        }

        MappedByteBuffer byteBlockBuffer = null;
        ByteBuffer lineBuffer = null; // buffer contains full line of traces
        int row = -1, col = -1;

        try
        {
            if(debugTimer)
            {
                StsSeismicTimer.clear();
                StsSeismicTimer.overallTimer.start();
                StsSeismicTimer.getOutputBlockTimer.start();
            }
//            filesMapBlocks[XDIR].unlockAllBlocks();
            lineBuffer = filesMapBlocks[XDIR].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE);
            if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
            if(isPlaneOK(XDIR, nPlane))
                return lineBuffer;

            if(debug)
                System.out.println("  ByteBuffer capacity: " + lineBuffer.capacity());
            byte[] data = new byte[nSlices];
            gather = new StsSuperGather(currentModel, preStackLineSet3d);

            preStackLineSet3d.checkTransparentTrace(nSlices);

            float maxZToStack = preStackLineSet3d.lineSetClass.getMaxZToStack();
            int nMaxDisplaySlices = Math.min(nSlices, preStackLineSet3d.getNearestBoundedSliceCoor(maxZToStack)) + 1;
            if(nMaxDisplaySlices < nSlices)
                System.arraycopy(preStackLineSet3d.byteTransparentTrace, nMaxDisplaySlices, data, nMaxDisplaySlices, nSlices - nMaxDisplaySlices);

            boolean[][] gridChanged = preStackLineSet3d.velocityModel.interpolation.gridChanged;
            boolean[][] isNeighbor = preStackLineSet3d.velocityModel.interpolation.isNeighbor;
            StsStatusArea statusArea = currentModel.win3d.statusArea;

            col = nPlane;
            int nUpdatedRows = 0;
            statusArea.setStatus("Calculating crossline " + getTypeName() + "...");
            statusArea.addProgress();
            statusArea.setMaximum(nRows);
            for(row = 0; row < nRows; row++)
            {
                statusArea.setProgress(row);
                if(isNeighbor[row][col] || !neighborsOnly)
                {
                    if(!gather.initializeSuperGather(row, col))
                    {
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                        continue;
                    }
                    float[] gatherData;
                    if(!applyAGCPoststack)
                        gatherData = computeTrace(gather, nMaxDisplaySlices);
                    else
                        gatherData = computeTraceAGC(gather, nMaxDisplaySlices, windowWidth, maxAmplitude);
                    if(gatherData == null)
                    {
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    else
                    {
                        nUpdatedRows++;
                        StsMath.scaleNormalize(gatherData);
                        scaleFloatsToBytes(gatherData, data);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(data);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    preStackLineSet3d.velocityModel.interpolation.setColStatus(row, col);
                }
                else // gridChanged && !isNeighbor && neighborsOnly
                {
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                    lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                }
            }
            if(!neighborsOnly) preStackLineSet3d.velocityModel.interpolation.setColStatus(col);
            statusArea.removeProgress();
            statusArea.clearStatus();

            if(debug) System.out.println("Updated " + nUpdatedRows + " out of " + nRows + " rows for col " + col);

            setPlaneOK(XDIR, nPlane, true);
            if(debug) System.out.println("Setting status to OK for dirNo " + XDIR + " plane " + nPlane);
//			byteBlockBuffer.force();
            lineBuffer.rewind();
            if(debugTimer)
            {
                StsSeismicTimer.printTimers("dirNo " + XDIR + " plane " + nPlane);
                StsSeismicTimer.overallTimer.stopPrint("Total time: ");
            }
            clearInputVolumes();
            // unlocking may be safer once the plane is no longer used in the display
//            filesMapBlocks[XDIR].unlockPlane(nPlane);
            return lineBuffer;
        }
        catch(Exception e)
        {
            if(lineBuffer != null)
                StsException.systemError(this, "getComputeColPlaneByteBuffer", "byteBuffer capacity: " + lineBuffer.capacity() + " position: " + lineBuffer.position());
            else
                StsException.systemError(this, "getComputeColPlaneByteBuffer", "null lineBuffer.");
            return null;
        }
    }

    public void unlockPlane(int dir, int nPlane)
    {
        filesMapBlocks[dir].unlockPlane(nPlane);   
    }

    protected ByteBuffer getComputeRowPlaneByteBuffer(int nPlane, boolean neighborsOnly)
    {
       StsAGCPreStackProperties agcProperties = preStackLineSet3d.agcProperties;
        boolean applyAGCPoststack = agcProperties.getApplyAGCPoststack();
        int windowWidth = 0;
        float maxAmplitude = 0.0f;
        if(applyAGCPoststack)
        {
            windowWidth = agcProperties.getWindowWidth(StsFilterProperties.POSTSTACK);
            maxAmplitude = Math.max(Math.abs(dataMin), Math.abs(dataMax));
        }

        MappedByteBuffer byteBlockBuffer = null;
        ByteBuffer lineBuffer = null; // buffer contains full line of traces
        int row = -1, col = -1;

        try
        {
            if(debugTimer)
            {
                StsSeismicTimer.clear();
                StsSeismicTimer.overallTimer.start();
                StsSeismicTimer.getOutputBlockTimer.start();
            }
//            filesMapBlocks[YDIR].unlockAllBlocks();
            lineBuffer = filesMapBlocks[YDIR].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE);
            if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
            if(isPlaneOK(YDIR, nPlane))
                return lineBuffer;

            if(debug)
                System.out.println("  ByteBuffer capacity: " + lineBuffer.capacity());
            byte[] data = new byte[nSlices];
            gather = new StsSuperGather(currentModel, preStackLineSet3d);

            preStackLineSet3d.checkTransparentTrace(nSlices);

            float maxZToStack = preStackLineSet3d.lineSetClass.getMaxZToStack();
            int nMaxDisplaySlices = Math.min(nSlices, preStackLineSet3d.getNearestBoundedSliceCoor(maxZToStack)) + 1;
            if(nMaxDisplaySlices < nSlices)
                System.arraycopy(preStackLineSet3d.byteTransparentTrace, nMaxDisplaySlices, data, nMaxDisplaySlices, nSlices - nMaxDisplaySlices);

            boolean[][] gridChanged = preStackLineSet3d.velocityModel.interpolation.gridChanged;
            boolean[][] isNeighbor = preStackLineSet3d.velocityModel.interpolation.isNeighbor;
            StsStatusArea statusArea = currentModel.win3d.statusArea;

            row = nPlane;
            int nUpdatedCols = 0;

            FloatBuffer floatBlockBuffer = null;
            if(preStackLineSet3d.isDataFloat)
            {
                floatBlockBuffer = fileMapRowFloatBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE).asFloatBuffer();
            }

            statusArea.setStatus("Calculating inline " + getTypeName() + "...");
            statusArea.addProgress();
            statusArea.setMaximum(nCols);
            for(col = 0; col < nCols; col++)
            {
                statusArea.setProgress(col);
                if(isNeighbor[row][col] || !neighborsOnly)
                {
                    if(!gather.initializeSuperGather(row, col))
                    {
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                        if(preStackLineSet3d.isDataFloat)
                            floatBlockBuffer.put(preStackLineSet3d.floatTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                        continue;
                    }
                    float[] gatherData;
                    if(!applyAGCPoststack)
                        gatherData = computeTrace(gather, nMaxDisplaySlices);
                    else
                        gatherData = computeTraceAGC(gather, nMaxDisplaySlices, windowWidth, maxAmplitude);

                    if(gatherData == null)
                    {
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                        if(preStackLineSet3d.isDataFloat)
                            floatBlockBuffer.put(preStackLineSet3d.floatTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    else
                    {
                        nUpdatedCols++;
                        StsMath.scaleNormalize(gatherData);
                        scaleFloatsToBytes(gatherData, data);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(data);
                        if(preStackLineSet3d.isDataFloat)
                            floatBlockBuffer.put(gatherData);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    preStackLineSet3d.velocityModel.interpolation.setRowStatus(row, col);
                }
                else // PLANE_STATUS_DIRTY && gridChanged && !isNeighbor && neighborsOnly
                {
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                    lineBuffer.put(preStackLineSet3d.byteTransparentTrace);
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                }
            }
            if(!neighborsOnly) preStackLineSet3d.velocityModel.interpolation.setRowStatus(row);
            statusArea.removeProgress();
            statusArea.clearStatus();

            if(debug) System.out.println("Updated " + nUpdatedCols + " out of " + nCols + " cols for row " + row);

            setPlaneOK(YDIR, nPlane, true);
            if(debug) System.out.println("Setting status to OK for dirNo " + YDIR + " plane " + nPlane);
//            byteBlockBuffer.force();
            lineBuffer.rewind();
            if(debugTimer)
            {
                StsSeismicTimer.printTimers("dirNo " + YDIR + " plane " + nPlane);
                StsSeismicTimer.overallTimer.stopPrint("Total time: ");
            }
            clearInputVolumes();
            // unlocking may be safer once the plane is no longer used in the display
//            filesMapBlocks[YDIR].unlockPlane(nPlane);
            return lineBuffer;
        }
        catch(Exception e)
        {
            if(lineBuffer != null)
                StsException.outputWarningException(this, "getComputeRowPlaneByteBuffer", "byteBuffer capacity: " + lineBuffer.capacity() + " position: " + lineBuffer.position(), e);
            else
                StsException.outputWarningException(this, "getComputeRowPlaneByteBuffer", "null lineBuffer.", e);
            return null;
        }
    }

    protected ByteBuffer getSlicePlaneByteBuffer(int nPlane, boolean neighborsOnly)
    {
        if(!isPlaneOK(ZDIR, nPlane)) return null;
        ByteBuffer lineBuffer = null; // buffer contains full line of traces
        try
        {
            if(debugTimer)
            {
                StsSeismicTimer.clear();
                StsSeismicTimer.overallTimer.start();
                StsSeismicTimer.getOutputBlockTimer.start();
            }
            lineBuffer = filesMapBlocks[ZDIR].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE);
            if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
            return lineBuffer;
        }
        catch(Exception e)
        {
            if(lineBuffer != null)
                StsException.outputWarningException(this, "getSlicePlaneByteBuffer", "byteBuffer capacity: " + lineBuffer.capacity() + " position: " + lineBuffer.position(), e);
            else
                StsException.outputWarningException(this, "getSlicePlaneByteBuffer", "null lineBuffer.", e);
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
        super.initialize(model, "rw");
        initializePlanesOK();
        return true;
    }

    public boolean initialize()
    {
        nRows = preStackLineSet3d.nRows;
        nCols = preStackLineSet3d.nCols;
        nSlices = preStackLineSet3d.nSlices;
        stsDirectory = preStackLineSet3d.stsDirectory;
        stemname = createStemname(preStackLineSet3d.stemname, type);
        isDataFloat = true;
        createVolumeFilenames(stemname);
        createFloatRowVolumeFilename(stemname);
        return initialize(currentModel, "rw");
    }

    static public String createStemname(String volumeStemname, byte type)
    {
        String typename = typenames[type];
        return currentModel.getProject().getName() + "." + volumeStemname + "." + typename;
    }

    public float[] computeTrace(StsSuperGather gather, int nValues)
    {
        if(type == VELOCITY)
            return preStackLineSet3d.velocityModel.computeInterpolatedFloatVelocityProfile(gather.superGatherRow, gather.superGatherCol, nValues);
        else if(type == SEMBLANCE)
            return gather.centerGather.computeSemblanceTrace(nValues);
        else
            return gather.centerGather.computeStackedTrace(nValues);
    }

    public float[] computeTraceAGC(StsSuperGather gather, int nValues, int windowWidth, float maxAmplitude)
    {
        if(type == VELOCITY)
            return preStackLineSet3d.velocityModel.computeInterpolatedFloatVelocityProfile(gather.superGatherRow, gather.superGatherCol, nValues);
        else if(type == SEMBLANCE)
            return gather.computeSemblanceTrace(nValues);
        else
        {
            float[] values = gather.centerGather.computeStackedTrace(nValues);
            return StsSeismicFilter.applyAGC(values, 0, 0, nValues, windowWidth, maxAmplitude);
        }
    }
/*
    public boolean scaleTrace(StsSuperGather gather, int nValues, byte[] byteData)
                  {
                      if(type == VELOCITY)
                          return preStackLineSet3d.velocityModel.computeInterpolatedByteVelocityProfile(gather.superGatherRow, gather.superGatherCol, byteData);
                      else if(type == SEMBLANCE)
                          return gather.computeSemblanceByteTrace(nValues, byteData);
                      else
                          return gather.centerGather.computeStackedFilteredByteTrace(nValues, byteData);
                  }
              */
/*
    public boolean scaleTrace(float[] floatValues, byte[] byteData)
    {
        if(type == VELOCITY)
            return preStackLineSet3d.velocityModel.scaleFloatsToBytes(floatValues, byteData);
        else if(type == SEMBLANCE)
            return scaleNormalizedFloatsToBytes(floatValues, byteData);
        else
            return normalizeAndScaleFloatsToBytes(floatValues, byteData);
    }
*/
    public ByteBuffer readByteBufferPlane(int dir, float dirCoordinate)
    {
        try
        {
            if(dir == ZDIR && currentAttribute != null && currentAttribute != nullAttribute)
                return currentAttribute.getByteBuffer();
            else
            {
                //ToDo: Need to scale to 254.
                int nPlane = this.getCursorPlaneIndex(dir, dirCoordinate);
                if(nPlane == -1) return null;
                switch(preStackLineSet3d.lineSetClass.stackOption)
                {
                    case StsPreStackLineSet.STACK_NONE:
                        return null;
                    case StsPreStackLineSet.STACK_NEIGHBORS:
                        ByteBuffer lineBuffer = getByteData(dir, nPlane, true);
                        if(lineBuffer == null)
                            return null;
                        else
                            return getNeighborBuffer(dir, nPlane, lineBuffer);
                    case StsPreStackLineSet.STACK_LINES:
                        return getByteData(dir, nPlane, false);
                    default:
                        return null;
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public void checkTransparentTrace(int nSlices)
    {
        if(byteTransparentTrace != null && byteTransparentTrace.length > nSlices) return;
        byteTransparentTrace = new byte[nSlices];
        floatTransparentTrace = new float[nSlices];
        for(int n = 0; n < nSlices; n++)
        {
            byteTransparentTrace[n] = -1;
            floatTransparentTrace[n] = 0;
        }
        // TODO: classInitialize floatTransparentTrace to appropriate value
    }


    private float[] getNeighborBuffer(int dir, int nPlane, float[] lineBuffer)
    {
        int row = -1, col = -1;
        try
        {
            checkTransparentTrace(nSlices);
            float[] floats = new float[nSlices];
            boolean[][] isNeighbor = preStackLineSet3d.velocityModel.interpolation.isNeighbor;
            int srcPos = 0;
            if(dir == XDIR)
            {
                col = nPlane;
                for(row = 0; row < nRows; row++)
                {
                    if(isNeighbor[row][col])
                        System.arraycopy(lineBuffer, srcPos, lineBuffer, srcPos, nSlices);
                    else
                        System.arraycopy(floatTransparentTrace, 0, lineBuffer, srcPos, nSlices);
                    srcPos += nSlices;
                }
            }
            else if(dir == YDIR)
            {
                row = nPlane;
                for(col = 0; col < nCols; col++)
                {
                    if(isNeighbor[row][col])
                        System.arraycopy(lineBuffer, srcPos, lineBuffer, srcPos, nSlices);
                    else
                        System.arraycopy(floatTransparentTrace, 0, lineBuffer, srcPos, nSlices);
                    srcPos += nSlices;
                }
            }
            return lineBuffer;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackLineSet3d.getNeighborBuffer() failed: null neighborBuffer.", e, StsException.WARNING);
            return null;
        }
    }

    private ByteBuffer getNeighborBuffer(int dir, int nPlane, ByteBuffer lineBuffer)
    {
        int row = -1, col = -1;
        ByteBuffer neighborBuffer = ByteBuffer.allocateDirect(lineBuffer.capacity());

        try
        {
            lineBuffer.rewind();
            neighborBuffer.rewind();
            checkTransparentTrace(nSlices);
            byte[] bytes = new byte[nSlices];
            boolean[][] isNeighbor = preStackLineSet3d.velocityModel.interpolation.isNeighbor;
            int position = 0;
            if(dir == XDIR)
            {
                //            progressBarDialog.setLabelText("Stacking Crossline #" + this.getColNumFromCol(index));
                //            progressBarDialog.setProgressMax(nRows);
                //            progressBarDialog.pack();
                //			  progressBarDialog.setVisible(true);
                col = nPlane;
                for(row = 0; row < nRows; row++)
                {
                    position += nSlices;
                    if(isNeighbor[row][col])
                    {
                        lineBuffer.get(bytes);
                        neighborBuffer.put(bytes);
                    }
                    else
                    {
                        neighborBuffer.put(byteTransparentTrace);
                        lineBuffer.position(position);
                    }
                }
            }
            else if(dir == YDIR)
            {
                //            progressBarDialog.setLabelText("Stacking Inline #" + this.getRowNumFromRow(index));
                //            progressBarDialog.setProgressMax(nCols);
                //            progressBarDialog.pack();
                //			  progressBarDialog.setVisible(true);

                row = nPlane;
                for(col = 0; col < nCols; col++)
                {
                    position += nSlices;
                    if(isNeighbor[row][col])
                    {
                        lineBuffer.get(bytes);
                        neighborBuffer.put(bytes);
                    }
                    else
                    {
                        neighborBuffer.put(byteTransparentTrace);
                        lineBuffer.position(position);
                    }
                }
            }
            neighborBuffer.rewind();
            return neighborBuffer;
        }
        catch(Exception e)
        {
            if(neighborBuffer != null)
                StsException.systemError("StsPreStackLineSet3d.getNeighborBuffer() failed. \n" +
                        "byteBuffer capacity: " + neighborBuffer.capacity() + " position: " + neighborBuffer.position());
            else
                StsException.outputException("StsPreStackLineSet3d.getNeighborBuffer() failed: null neighborBuffer.", e, StsException.WARNING);
            return null;
        }
    }

    /**
     * getByteData
     *
     * @param dir    int
     * @param nPlane int
     * @return ByteBuffer
     */
    public ByteBuffer getByteData(int dir, int nPlane, boolean neighborsOnly)
    {
        if(preStackLineSet3d.velocityModel == null)
            return null;
        if(!preStackLineSet3d.velocityModel.hasProfiles())
            return null;
        if(!preStackLineSet3d.velocityModel.checkRunInterpolation())
            return null;

        switch(dir)
        {
            case XDIR:
                return getComputeColPlaneByteBuffer(nPlane, neighborsOnly);
            case YDIR:
                return getComputeRowPlaneByteBuffer(nPlane, neighborsOnly);
            case ZDIR:
                return getSlicePlaneByteBuffer(nPlane, neighborsOnly);
            default:
                return null;
        }
    }

    public boolean isPlaneOK(int dir, int plane)
    {
        if(plane == -1) return true;
        return planesOK[dir][plane];
    }

    private void initializePlanesOK()
    {
        planesOK = new boolean[3][];
        planesOK[0] = colPlanesOK;
        planesOK[1] = rowPlanesOK;
        planesOK[2] = slicePlanesOK;
    }


    public void setPlaneOK(int dir, int plane, boolean isOK)
    {
        boolean wasOK = planesOK[dir][plane];
        if(wasOK != isOK)
        {
            planesOK[dir][plane] = isOK;
            volumeChanged = true;
//			StsArrayElementChangeCmd changeCmd = new StsArrayElementChangeCmd(this, new Boolean(wasOK), new Boolean(isOK), "planeOK", new int[] {dirNo, plane});
//			currentModel.addTransactionCmd("processVolume planeOK[" + dirNo + ", " + plane + "] set to" + isOK, changeCmd);
        }
    }

    /** check if planeOK flag arrays exist for this volume;
     *  if not initialize 3 sets of planeOK boolean flags for this volume.
     *  Return true if initialized.
     */
    public boolean checkInitializePlaneOKFlags()
    {
        boolean initialized = false;
        if(colPlanesOK == null)
        {
            colPlanesOK = new boolean[preStackLineSet3d.nCols];
            rowPlanesOK = new boolean[preStackLineSet3d.nRows];
            slicePlanesOK = new boolean[preStackLineSet3d.nSlices];
            initialized = true;
        }
        if(planesOK == null)
        {
            initializePlanesOK();
        }
        return initialized;
    }

    public void validateAllPlanes()
    {
//		StsDBMethodCmd cmd = new StsDBMethodCmd(this, "setAllPlaneOKFlagsTrue");
//		currentModel.addTransactionCmd("processVolume setAllPlaneOKFlagsTrue", cmd);
        setAllPlaneOKFlagsTrue();
        volumeChanged = true;
    }

    public void setAllPlaneOKFlagsFalse()
    {
        if(checkInitializePlaneOKFlags()) return;
        Arrays.fill(colPlanesOK, false);
        Arrays.fill(rowPlanesOK, false);
        Arrays.fill(slicePlanesOK, false);
        volumeChanged = true;
    }

    public void setAllPlaneOKFlagsTrue()
    {
        checkInitializePlaneOKFlags();
        Arrays.fill(colPlanesOK, true);
        Arrays.fill(rowPlanesOK, true);
        Arrays.fill(slicePlanesOK, true);
        volumeChanged = true;
    }

    public void setDirPlanesOK(int dir, boolean ok)
    {
        Arrays.fill(planesOK[dir], ok);
    }

    /** clear blocks in memory other than those in the output file set */
    protected void clearInputVolumes()
    {
        StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
        if(filesMapBlocks != null) memoryManager.clearOtherBlocks(filesMapBlocks);
	}
}
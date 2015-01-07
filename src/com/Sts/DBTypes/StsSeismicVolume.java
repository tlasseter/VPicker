package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.text.*;
import java.util.*;

public class StsSeismicVolume extends StsSeismic implements StsTreeObjectI, StsVolumeDisplayable
{
    transient StsFile out = null;
    transient OutputStream os = null;
    transient BufferedOutputStream bos = null;
    transient ByteArrayOutputStream baos = null;
    transient DataOutputStream ds = null;

    transient float exportScale = 1.0f;
    transient String exportTextHeader = null;

    transient StsBlocksMemoryManager blocksMemory = null;
    transient public StsCubeFileBlocks[] filesMapBlocks;
    transient public StsCubeFileBlocks fileMapRowFloatBlocks;

    /** These are textures of this seismic volume on surfaces.  Texture data is filled only when needed. */
    transient StsSurfaceTexture[] surfaceTextures;

    /** if true, the calculateDataPlaneMethod will be called to compute the data */
    transient boolean calculateDataPlane = false;
    transient StsMethod calculateDataPlaneMethod = null;

    //	transient InterpolatedTrace interpolatedTrace = null;

    transient int displayListNum = 0;

    /** Display lists should be used (controlled by View:Display Options) */
    transient boolean useDisplayLists;
    /** Display lists currently being used for surface geometry */
    transient boolean usingDisplayLists = true;

    transient IsoValueGridSet isoValueGrids = null;
    transient IsoValueGridSet isoValueGridsRow0 = null;
    transient IsoValueGridSet isoValueGridsRow1 = null;

    transient VoxelVertexArrays[] voxelVertexArraySets;
    static final int arraySize = 500000; // maximum number of vertices allowed in each VoxelVertexArray

    transient int isoSurfaceID = 0;
    transient int nIsoValueGrids = 0;
    /** don't display voxels smaller than this size (0: draw all) */
    transient int smallestVoxel = 0;
    /*
        static public String DISPLAY_ATTRIBUTE_SEISMIC = "Seismic";
        static public String[] displayAttributes = new String[]
            {
            DISPLAY_ATTRIBUTE_SEISMIC};
    */

    /** Has the basemap parameters changed? */
    transient boolean basemapDisplayChanged = false;
    /** What attribute is isVisible on the basemap */
    //	transient String displayAttribute = DISPLAY_ATTRIBUTE_SEISMIC;

    /**
     * A vector of the currently active voxelKeys which are defined on the colorscale.
     * A VoxelKey has a min and a max value in the colorscale and defines the range over which a voxel
     * is isVisible.  We actually only draw the min and the max surfaces of the voxel blobs.
     * To know what blob we are in, we fill a byte vector with -1 if the color is outside a voxelKey
     * range and the voxelKey index ( >= 0) if it is inside a voxelKey range.
     */
    transient StsColorscale.VoxelKey[] voxelKeys = null;

    /**
     * For each possible voxel sample byte value, this is an index into the vectorColorIndexes array (if >= 0)
     * which in turn is the index into the colorscale vector.
     * If the index is -1, the corresponding sample byte value is outside the keyRanges.
     */
    transient int[] voxelKeyIndexes;

    /**
     * Voxels in a voxelKey range can have either the top color or the bottom color of the range.
     * We create a set of vertexArrays for each possible voxel color.  The number of such sets is
     * the number of voxelKeys times 2.
     */
    transient int[] voxelColorIndexes;

    /** experimental drawmode switch */
    transient boolean drawQuads = false;
    static protected final float[] nullNormalVector = new float[]{nullValue, nullValue, nullValue};

    static public final int DATA = 3;
    static public final int CROP_CHANGED = 0;

    static public final byte nullByte = StsParameters.nullByte;
    static public final int nullUnsignedInt = 255;

    static boolean hasVBO = false;
    static boolean vboChecked = false;
    static boolean hasPP = false;
    static boolean ppChecked = false;

    static protected StsObjectPanel objectPanel = null;

    static StsComboBoxFieldBean displayAttributeBean = new StsComboBoxFieldBean(StsSeismicVolume.class, "displayAttribute", "Basemap Attribute");

    static StsDateFieldBean bornField = new StsDateFieldBean(StsSeismicVolume.class, "bornDate", "Born Date:");
    static StsDateFieldBean deathField = new StsDateFieldBean(StsSeismicVolume.class, "deathDate", "Death Date:");
    static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsSeismicVolume.class, "colorscale");

    static public final StsFieldBean[] seismicDisplayFields =
        {
            new StsBooleanFieldBean(StsSeismicVolume.class, "isVisible", "Enable"),
            bornField, deathField,
            new StsBooleanFieldBean(StsSeismicVolume.class, "readoutEnabled",
                "Mouse Readout"),
            displayAttributeBean,
             colorscaleBean

        };

    static public final StsFieldBean[] seismicPropertyFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsSeismicVolume.class, "name", true, "Name"),
            new StsStringFieldBean(StsSeismicVolume.class, "zDomainString", false, "Z Domain"),
            new StsStringFieldBean(StsSeismicVolume.class, "segyFilename", false, "SEGY Filename"),
            new StsStringFieldBean(StsSeismicVolume.class, "segyFileDate", false, "SEGY creation date"),
            new StsIntFieldBean(StsSeismicVolume.class, "nRows", false, "Number of Lines"),
            new StsIntFieldBean(StsSeismicVolume.class, "nCols", false, "Number of Crosslines"),
            new StsIntFieldBean(StsSeismicVolume.class, "nSlices", false, "Number of Samples"),
            new StsDoubleFieldBean(StsSeismicVolume.class, "xOrigin", false, "X Origin"),
            new StsDoubleFieldBean(StsSeismicVolume.class, "yOrigin", false, "Y Origin"),
            new StsFloatFieldBean(StsSeismicVolume.class, "xInc", false, "X Inc"),
            new StsFloatFieldBean(StsSeismicVolume.class, "yInc", false, "Y Inc"),
            new StsFloatFieldBean(StsSeismicVolume.class, "zInc", false, "Z Inc"),
            new StsFloatFieldBean(StsSeismicVolume.class, "xMin", false, "X Loc Min"),
            new StsFloatFieldBean(StsSeismicVolume.class, "yMin", false, "Y Loc Min"),
            new StsFloatFieldBean(StsSeismicVolume.class, "zMin", false, "Min Z or T"),
            new StsFloatFieldBean(StsSeismicVolume.class, "zMax", false, "Max Z or T"),
            new StsFloatFieldBean(StsSeismicVolume.class, "angle", false, "Angle to Line Direction"),
            new StsFloatFieldBean(StsSeismicVolume.class, "rowNumMin", false, "Min Line"),
            new StsFloatFieldBean(StsSeismicVolume.class, "rowNumMax", false, "Max Line"),
            new StsFloatFieldBean(StsSeismicVolume.class, "colNumMin", false, "Min Crossline"),
            new StsFloatFieldBean(StsSeismicVolume.class, "colNumMax", false, "Max Crossline"),
            new StsFloatFieldBean(StsSeismicVolume.class, "dataMin", false, "Data Min"),
            new StsFloatFieldBean(StsSeismicVolume.class, "dataMax", false, "Data Max"),
            new StsFloatFieldBean(StsSeismicVolume.class, "dataAvg", false, "Data Avg"),
        };

    static StsTimer timer = null;
    static boolean runTimer = false;

    static final boolean debug = false;
    static final boolean voxelDebug = false;

    public String getGroupname()
    {
        return group3d;
    }

    public StsSeismicVolume()
    {
    }

    public StsSeismicVolume(boolean persistent)
    {
        super(persistent);
    }

    public StsSeismicVolume(StsModel model, StsSeismicBoundingBox seismicBoundingBox, StsRotatedGridBoundingBox cropBox, String mode)
    {
        this(false);
        StsToolkit.copySubToSuperclass(seismicBoundingBox, this, StsSeismicBoundingBox.class, StsSeismicBoundingBox.class, true);
        StsToolkit.copySubToSuperclass(cropBox, this, StsRotatedGridBoundingBox.class, StsBoundingBox.class, true);
        isDataFloat = true;
        setName(stemname);
        allocateVolumes(mode);
    }

    /*
        public StsSeismicVolume(StsFile file, StsModel model) throws FileNotFoundException, StsException
        {
            this(file.getDirectory(), file.getFilename(), model);
        }
    */
    /** construct this seismicVolume from files.  Find/read the header file. initialize(model) will find/read the volume files. */
    protected StsSeismicVolume(StsModel model, StsFile file) throws FileNotFoundException, StsException
    {
        super(false);
        if (!file.exists())
        {
            throw new FileNotFoundException();
        }
        loadFromFile(file.getDirectory(), file.getFilename(), model);
    }

    protected void loadFromFile(String directory, String filename, StsModel model) throws FileNotFoundException, StsException
    {
        String pathname = directory + filename;
        StsParameterFile.initialReadObjectFields(pathname, this, getLoadSubClass(), getLoadSuperClass());
        // restore stsDirectory to directory since it is overwritten by read above
        stsDirectory = directory;
        setName(getStemname());
        if (!initialize(model))
        {
            throw new FileNotFoundException(pathname);
        }
        isVisible = true;
        getSeismicVolumeClass().setIsVisibleOnCursor(true);
    }

    static public StsSeismicVolume constructor(StsFile file, StsModel model, boolean addToModel) throws FileNotFoundException, StsException
    {
        return StsSeismicVolume.checkLoadFromFilename(model, file, addToModel);
    }

    static public StsSeismicVolume constructTestVolume(String filePathname)throws FileNotFoundException, StsException
    {
        StsModel model = StsModel.constructor();
        StsFile file = StsFile.constructor(filePathname);
        return constructor(file, model, false);
    }

    static public StsSeismicVolume checkLoadFromStemname(StsModel model, String directory, String stemname, boolean addToModel)
    {
        String filename = createHeaderFilename(StsSeismicBoundingBox.group3d, stemname);
        try
        {
            StsFile file = StsFile.constructor(directory, filename);
            return checkLoadFromFilename(model, file, addToModel);
        }
        catch (FileNotFoundException fnfe)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find file " + filename);
            return null;
        }
        catch (StsException e)
        {
            StsException.outputException("StsSeismicVolume.checkloadFromStemname failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public StsSeismicVolume checkLoadFromFilename(StsModel model, StsFile file, boolean addToModel) throws FileNotFoundException, StsException
    {
        if (!file.exists()) return null;
        StsSeismicVolume seismicVolume = new StsSeismicVolume(model, file);
        if (addToModel)
        {
            if (!seismicVolume.addToProject(true))
            {
//                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add volume to project: " + file.getFilename());
                return null;
            }
            seismicVolume.addToModel();
        }
        return seismicVolume;
    }

    public boolean initialize(StsFile file, StsModel model)
    {
        try
        {
            String pathname = file.getDirectory() + file.getFilename();
            StsParameterFile.initialReadObjectFields(pathname, this, StsSeismicBoundingBox.class, StsBoundingBox.class);
            setName(getStemname());
            initializeColorscale();
            stsDirectory = file.getDirectory();
            if (!initialize(model))
            {
                return false;
            }
            isVisible = true;
            getSeismicVolumeClass().setIsVisibleOnCursor(true);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.loadFile() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public boolean initialize()
    {
        return initialize(currentModel);
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            return initialize(model, "r");
        }
        catch (Exception e)
        {
            new StsMessage(null, StsMessage.ERROR, "Failed to open Seismic line file: " + e.getMessage());
            StsMessage.printMessage("Failed to find file. Error: " + e.getMessage());
            return false;
        }
    }

    public boolean initialize(StsModel model, String mode)
    {
        super.initialize(model);
        blocksMemory = model.getProject().getBlocksMemoryManager();
        //		dataCubeMemory = model.project.getDataCubeMemory();
        clearCache();
        if (!allocateVolumes(mode))
        {
            return false;
        }
        initializeColorscale();
        setDataHistogram();
        byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
        zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
        //		if (!model.project.checkSetZDomain(zDomainByte, zDomainByte))
        //		{
        //			return false;
        //		}
        initDisplayAttributes();
        return true;
    }

	public void initializeColorscale()
	{
		try
		{
			if (colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSeismicClass seismicClass = (StsSeismicClass)getCreateStsClass();
				// jbw
				if (name != null && (name.contains("velocity") || this.name.contains("Vint") || this.name.contains("Vave")) )
					colorscale = new StsColorscale("Velocity", spectrumClass.getSpectrum(spectrumClass.SPECTRUM_RAINBOW), getVelocityScalar()*dataMin, getVelocityScalar()*dataMax);
				else
					colorscale = new StsColorscale("Seismic", spectrumClass.getSpectrum(seismicClass.getSeismicSpectrumName()), getVelocityScalar()*dataMin, getVelocityScalar()*dataMax);
				colorscale.setEditRange(getVelocityScalar()*dataMin, getVelocityScalar()*dataMax);
			}
			seismicColorList = new StsColorList(colorscale);
			colorscale.addActionListener(this);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

    public Class getLoadSubClass() { return StsSeismicBoundingBox.class; }

    public Class getLoadSuperClass() { return StsMainObject.class; }

    static public StsSeismicVolume initializeAttributeVolume(StsModel model, StsSeismicBoundingBox volume, float dataMin, float dataMax, boolean isDataFloat, boolean writeHeaderFile, String stemname, String attributeName, String mode)
    {
        stemname = stemname + "." + attributeName;
        return initializeAttributeVolume(model, volume, dataMin, dataMax, isDataFloat, writeHeaderFile, stemname, mode);
    }

    static public StsSeismicVolume initializeAttributeVolume(StsModel model, StsSeismicBoundingBox volume, float dataMin, float dataMax, boolean isDataFloat, boolean writeHeaderFile, String stemname, String mode)
    {
        StsSeismicVolume attributeVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, stemname);
        if (attributeVolume != null) attributeVolume.delete();
        attributeVolume = new StsSeismicVolume(false);
        byte zDomain = volume.getZDomain();
        attributeVolume.initializeVolume(model, volume, dataMin, dataMax, isDataFloat, writeHeaderFile, stemname, zDomain, false, mode);
        return attributeVolume;
    }

    /*
         public void initializeVolume(StsModel model, StsSeismicBoundingBox volume, boolean createFiles, String stemname, byte zDomain, boolean addToModel)
         {
            initializeVolume(model, volume, 0.0f, 0.0f, createFiles, stemname, zDomain, addToModel);
         }
    */
    public void initializeVolume(StsModel model, StsSeismicBoundingBox volume, float dataMin, float dataMax, boolean isDataFloat, boolean writeHeaderFile, String stemname, byte zDomain, boolean addToModel, String mode)
    {
        // Initialize this new volume (StsSeismicBoundingBox) by copying from "this" volume;
        // this has the side effect of overwriting some members in StsSeismicBoundingBox which need to be recomputed
        // such as row/col/slice filenames (from stemname) and dataMin, dataMax, and zDomain.
        StsToolkit.copySubToSuperclass(volume, this, StsSeismicBoundingBox.class, StsBoundingBox.class, false);

        stsDirectory = volume.stsDirectory;
        setIsDataFloat(isDataFloat);
        setName(stemname);
        setFilenames(stemname);
        deleteExistingFiles();
        setDataMin(dataMin);
        setDataMax(dataMax);
        setZDomain(zDomain);
        if (writeHeaderFile) writeHeaderFile();
        initializeColorscale();
        if (addToModel) addToModel();

        //         if (isDataFloat) createMappedRowBuffer("rw");
        allocateVolumes(mode);
    }

    public StsSeismicVolumeClass getSeismicVolumeClass()
    {
        return (StsSeismicVolumeClass) currentModel.getCreateStsClass(getClass());
    }

    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return getSeismicVolumeClass().getWiggleDisplayProperties();
    }

    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            int modifiers = e.getModifiers();
            if (modifiers == StsColorscale.VOXELS_CHANGED)
            {
                voxelVertexArraySets = null;
            }
            else
            {
                seismicColorList.setColorListChanged(true);
                currentModel.clearDisplayTextured3dCursors(this);
                fireActionPerformed(e);
            }
            currentModel.viewObjectRepaint(this, this);
        }
        else
        {
            fireActionPerformed(e);
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
        return;
    }

    public void addToModel()
    {
        super.addToModel();
        currentModel.checkAddToCursor3d(this);
    }

    public boolean getIsPixelMode()
    {
        return getSeismicVolumeClass().getIsPixelMode();
    }

    public void setDataHistogram()
    {
        if (dataHist == null || colorscaleBean == null) return;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    colorscaleBean.setHistogram(dataHist);
                    colorscaleBean.revalidate();
                }
            }
        );
    }

    public void setDataHistogram(float[] histogram)
    {
        dataHist = histogram;
        setDataHistogram();
    }

    public float[] getHistogram()
    {
        return dataHist;
    }

    public boolean allocateVolumes(String mode)
    {

        try
        {
            StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
            if(Main.isJarDB)
            {
               // Set file names
                String rootDir = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
                File rfile = Main.jar.uncompressTo(rootDir, rowCubeFilename);
                File cfile = Main.jar.uncompressTo(rootDir, colCubeFilename);
                File sfile = Main.jar.uncompressTo(rootDir, sliceCubeFilename);
                if((rfile == null) || (cfile == null) || (sfile == null))
                {
                    new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find all required volume files in archive for volume: " + this.getName());
                    return false;
                }
                if (filesMapBlocks == null)
                {
                    filesMapBlocks = new StsCubeFileBlocks[3];
                    filesMapBlocks[XDIR] = new StsCubeFileBlocks(XDIR, nRows, nCols, nSlices, 1, cfile, 1, memoryManager);
                    filesMapBlocks[YDIR] = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, 1, rfile, 1, memoryManager);;
                    filesMapBlocks[ZDIR] = new StsCubeFileBlocks(ZDIR, nRows, nCols, nSlices, 1, sfile, 1, memoryManager);
                }
                isDataFloat = rowFloatFilename != null;
                if (rowFloatFilename != null && fileMapRowFloatBlocks == null)
                try
                {
                    File fltFile = Main.jar.uncompressTo(rootDir,rowFloatFilename);
                    fileMapRowFloatBlocks = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, 1, fltFile, 4, memoryManager);
                }

                // if the float volume is not available, we will assume application can use the row byte volume (XDIR volume) allocated above
                // and compute floats as needed from the scaled bytes
                catch (Exception e)
                {
                    rowFloatFilename = null;
                    isDataFloat = false;
                }

            }
            else
            {
                if (filesMapBlocks == null)
                {
                    filesMapBlocks = new StsCubeFileBlocks[3];
                    filesMapBlocks[XDIR] = new StsCubeFileBlocks(XDIR, nRows, nCols, nSlices, 1, this, colCubeFilename, mode, 1, memoryManager);
                    filesMapBlocks[YDIR] = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, 1, this, rowCubeFilename, mode, 1, memoryManager);;
                    filesMapBlocks[ZDIR] = new StsCubeFileBlocks(ZDIR, nRows, nCols, nSlices, 1, this, sliceCubeFilename, mode, 1, memoryManager);
                }
                isDataFloat = rowFloatFilename != null;
                if (rowFloatFilename != null && fileMapRowFloatBlocks == null)
                try
                {
                    fileMapRowFloatBlocks = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, 1, this, rowFloatFilename, mode, 4, memoryManager);
                }

                // if the float volume is not available, we will assume application can use the row byte volume (XDIR volume) allocated above
                // and compute floats as needed from the scaled bytes
                catch (Exception e)
                {
                    rowFloatFilename = null;
                    isDataFloat = false;
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.allocateVolumes) failed.", e, StsException.WARNING);
            return false;
        }
    }

    public StsCubeFileBlocks[] getFilesMapBlocks() { return this.filesMapBlocks; }

    public boolean getSliceBlockTraceValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
    {
        return filesMapBlocks[ZDIR].getSliceBlockTraceValues(row, col, sliceMin, sliceMax, values);
    }

    public boolean getRowBlockTraceValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
    {
        return filesMapBlocks[YDIR].getRowBlockTraceValues(row, col, sliceMin, sliceMax, values);
    }

    public boolean getSlice(int nSlice, byte[] values)
    {
        return filesMapBlocks[ZDIR].getPlane(nSlice, values);
    }

    public ByteBuffer getSliceBuffer(int nSlice)
    {
        return filesMapBlocks[ZDIR].getByteBufferPlane(nSlice, FileChannel.MapMode.READ_ONLY).slice();
    }

    public boolean setSliceBufferPosition(ByteBuffer sliceBuffer, int row, int col)
    {
        try
        {
            sliceBuffer.position(row * nCols + col);
            return true;
        }
        catch (Exception e)
        {
            return false;
        }
    }

    public int getSlicePosition(int nSlice)
    {
        return filesMapBlocks[ZDIR].getPlanePosition(nSlice);
    }

    public void printMemorySummary()
    {
        long planeMemoryUsed = 0, sliceBlockMemoryUsed = 0;
        /*
                if (dataCubeMemory != null)
                {
                    dataCubeMemory.printMemorySummary();
                }
        */
        for (int n = 0; n < 3; n++)
            filesMapBlocks[n].printMemorySummary();
    }

    public boolean setupReadRowFloatBlocks()
    {
        if (fileMapRowFloatBlocks != null)
        {
            return true;
        }
        try
        {
            return constructReadFloatBlocks();
        }
        catch (Exception e)
        {
            StsMessage.printMessage("Couldn't find/read float volume cube: " + rowFloatFilename);
            return false;
        }
    }

    public boolean setupWriteRowFloatBlocks()
    {
        if (fileMapRowFloatBlocks != null)
        {
            return true;
        }
        try
        {
            return constructReadFloatBlocks();
        }
        catch (Exception e)
        {
            StsMessage.printMessage("Couldn't find/read float volume cube: " + rowFloatFilename);
            return false;
        }
    }

    public boolean setupRowByteBlocks()
    {
        try
        {
            StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
            filesMapBlocks[YDIR] = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, this, rowCubeFilename, "r", 1, memoryManager);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.setupRowByteBlocks() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public void deleteRowFloatBlocks()
    {
        if (fileMapRowFloatBlocks != null)
        {
            StsBlocksMemoryManager blocksMemoryManager = currentModel.getProject().getBlocksMemoryManager();
            blocksMemoryManager.clearBlocks(fileMapRowFloatBlocks);
            //			fileMapRowFloatBlocks.deleteAllBlocks();
            fileMapRowFloatBlocks = null;
        }
        return;
    }

    private boolean constructReadFloatBlocks() throws FileNotFoundException
    {
        return constructFloatBlocks("r");
    }

    private boolean constructWriteFloatBlocks() throws FileNotFoundException
    {
        return constructFloatBlocks("rw");
    }

    private boolean constructFloatBlocks(String mode) throws FileNotFoundException
    {
        try
        {
            StsBlocksMemoryManager blocksMemoryManager = currentModel.getProject().getBlocksMemoryManager();
            fileMapRowFloatBlocks = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, this, rowFloatFilename, mode, 4, blocksMemoryManager);
            return true;
        }
        catch (Exception e)
        {
            StsException.systemError("StsSeismicVolume.constructReadFloatBlocks failed.");
            return false;
        }
    }

    public int getNumberOfFilledCells(StsProgressPanel progressPanel)
    {
        int count = 0;
        int slices = 0;
        if(nSlices <= 0)
            return -1;
        for(float i=zMin; i<=zMax; i = i+zInc)
        {
            progressPanel.incrementCount();
            progressPanel.setDescription("Completed slice #" + (slices+1));
            byte[] plane = this.readBytePlaneData(ZDIR, i);
            if(plane == null)
                continue;
            for(int j=0; j<plane.length; j++)
            {
                if(plane[j] != -1)
                   count++;
            }
            slices++;
        }
        return count;
    }
    public float getYFromRowNum(float rowNum)
    {
        return (float) ((((rowNum - rowNumMin) / rowNumInc) * yInc) + (double) yMin);
    }

    public float getXFromColNum(float colNum)
    {
        return (float) ((((colNum - colNumMin) / colNumInc) * xInc) + (double) xMin);
    }

    public String getSegyFileDate()
    {
        if (segyLastModified == 0)
        {
            File segyFile = new File(segyDirectory + segyFilename);
            if (segyFile != null)
            {
                segyLastModified = segyFile.lastModified();
            }
        }
        DateFormat dateFormat = DateFormat.getDateTimeInstance();
        return dateFormat.format(new Date(segyLastModified));
    }

    public String getDate()
    {
        return null;
    }
	public void treeObjectSelected()
	{
        currentModel.setCurrentObject(this);
		currentModel.getGlPanel3d().checkAddView(StsView3d.class);
		currentModel.win3dDisplayAll();
	}

	public boolean anyDependencies()
	{
		StsCrossplot[] cp = (StsCrossplot[])currentModel.getCastObjectList(
			StsCrossplot.class);
		for (int n = 0; n < cp.length; n++)
		{
			StsSeismicBoundingBox[] volumes = cp[n].getVolumes();
			for (int j = 0; j < cp[n].volumes.getSize(); j++)
			{
				if (this == cp[n].volumes.getElement(j))
				{
					StsMessageFiles.infoMessage("Seismic PostStack3d " + getName() +
												" used by Crossplot " +
												cp[n].getName());
					return true;
				}
			}
		}
		return false;
	}

	public boolean delete()
	{
		close();
        super.delete(); // sets this.index to -1
        return true;
	}

    public String getLabel()
    {
        return stemname;
    }

    public int getIntValue(int row, int col, int slice)
    {
        if (row < 0 || row >= nRows || col < 0 || col >= nCols || slice < 0 || slice >= nSlices)
        {
            StsException.systemError("StsSeismicVolume.getValue() failed for row: " + row + " col: " + col + " plane: " + slice);
            return 0;
        }

        byte[] rowData = readBytePlaneData(StsCursor3d.YDIR, getYCoor(row));
        if (rowData == null)
        {
            return 0;
        }
        byte signedByteValue = rowData[col * nSlices + slice];
        return StsMath.signedByteToUnsignedInt(signedByteValue);
    }

    public int getIntValue(float x, float y, float z)
    {
        int row = getNearestBoundedRowCoor(y);
        int col = getNearestBoundedColCoor(x);
        int plane = getNearestBoundedSliceCoor(z);
        return getIntValue(row, col, plane);
    }

    public int getIntValue(float[] xyz)
    {
        int row = getNearestBoundedRowCoor(xyz[1]);
        int col = getNearestBoundedColCoor(xyz[0]);
        int plane = getNearestBoundedSliceCoor(xyz[2]);
        return getIntValue(row, col, plane);
    }

    public int getIntValue(int row, int col, float z)
    {
        int v0, v1, v;

        float kF = (float) ((z - zMin) / zInc);
        int k = (int) kF;
        float dk = kF - k;
        if (dk < 0.5f)
        {
            return getIntValue(row, col, k);
        }
        else
        {
            return getIntValue(row, col, k + 1);
        }
    }

    public float getScaledByteValue(int row, int col, float z)
    {
        float kF = (z - zMin) / zInc;
        int k = (int) kF;
        float dk = kF - k;
        int v0 = getIntValue(row, col, k);
        int v1 = getIntValue(row, col, k + 1);
        float v = v0 + dk * (v1 - v0);
        return dataMin + (v / 254) * (dataMax - dataMin);
    }

    public int getPlaneValue(int dir, float[] xyz)
    {
        int row = -1, col = -1, slice = -1;
        byte signedByteValue;
        byte[] planeData;
        try
        {
            row = getNearestBoundedRowCoor(xyz[1]);
            col = getNearestBoundedColCoor(xyz[0]);
            slice = getNearestBoundedSliceCoor(xyz[2]);
            signedByteValue = filesMapBlocks[YDIR].readByteValue(row, col, slice);
            return StsMath.signedByteToUnsignedInt(signedByteValue);
        }
        catch (Exception e)
        {
            StsException.systemError(
                "StsSeismicVolume.getPlaneValue() failed for row " + row +
                    " col " + col +
                    " slice " +
                    slice);
            return 0;
        }
    }

    public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, byte[] values)
    {
        return filesMapBlocks[YDIR].getRowCubeByteValues(row, col, sliceMin, sliceMax, values);
    }

    public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax,
                                  int dir, boolean useByteCubes,
                                  double[] values)
    {
        int nValues = sliceMax - sliceMin + 1;
        if (nValues <= 0) return false;

        try
        {
            if (!useByteCubes && fileMapRowFloatBlocks != null)
            {
                float[] floatData = new float[nValues];
                if (!fileMapRowFloatBlocks.getRowCubeFloatValues(row, col, sliceMin, sliceMax, floatData))
                    return false;
                for (int n = 0; n < nValues; n++)
                    values[n] = (double) floatData[n];
            }
            else
            {
                byte[] byteData = new byte[nValues];
                if (!filesMapBlocks[YDIR].getRowCubeByteValues(row, col, sliceMin, sliceMax, byteData))
                    return false;
                for (int n = 0; n < nValues; n++)
                    values[n] = getScaledValue(byteData[n]);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.systemError("StsSeismicVolume.getValue() failed for row " + row + " col " + col + " sliceMin " + sliceMin +
                " sliceMax " + sliceMax + "values.length " + values.length);
            return false;
        }
    }

    public float[] getTraceValues(int row, int col)
    {
        float[] floatData = new float[nSlices];
        if (!getTraceValues(row, col, 0, nSlices - 1, StsCursor3d.XDIR, true, floatData)) return null;
        return floatData;
    }

    public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, boolean useByteCubes, float[] floatData)
    {
        int index, n;
        int nValues = sliceMax - sliceMin + 1;
        if (nValues <= 0)
        {
            StsException.systemError("StsSeismicVolume.getValues(....floats) called with nValues negative: " + nValues);
            return false;
        }
        if (floatData == null)
        {
            StsException.systemError("StsSeismicVolume.getValues(....floats) called with floatData = null.");
            return false;
        }
        if (nValues > floatData.length)
        {
            StsException.systemError("StsSeismicVolume.getValues(....floats) called with nValues (" + nValues + ") > floatData.length (" + floatData.length + ").");
            return false;
        }

        try
        {
            if (!useByteCubes && fileMapRowFloatBlocks != null)
            {
                return fileMapRowFloatBlocks.getRowCubeFloatValues(row, col, sliceMin, sliceMax, floatData);
            }
            else
            {
                byte[] planeData;
                switch (dir)
                {
                    case XDIR:
                        /*
                            if (mainDebug && !planeInMemory(XDIR, col))
                            {
                                System.out.println("Loaded 1 xdir planes at row: " +
                                                   row + " col: " + col);
                            }
                        */
                        planeData = readBytePlaneData(XDIR, getXCoor(col));
                        index = row * nSlices + sliceMin;
                        n = 0;
                        for (int slice = sliceMin; slice <= sliceMax; slice++,
                            index++)
                        {
                            floatData[n++] = getScaledValue(planeData[index]);
                        }
                        break;
                    case YDIR:
                        /*
                            if (mainDebug && !planeInMemory(YDIR, row))
                            {
                                System.out.println("Loaded 1 ydir planes at row: " +
                                                   row + " col: " + col);
                            }
                        */
                        planeData = readBytePlaneData(YDIR, getYCoor(row));
                        index = col * nSlices + sliceMin;
                        n = 0;
                        for (int slice = sliceMin; slice <= sliceMax; slice++,
                            index++)
                        {
                            floatData[n++] = getScaledValue(planeData[index]);
                        }
                        break;
                    case ZDIR:
                        n = 0;
                        index = row * nCols + col;
                        int nPlanesToLoad = 0;
                        for (int slice = sliceMin; slice <= sliceMax; slice++)
                        {
                            /*
                                if (mainDebug && !planeInMemory(ZDIR, slice))
                                {
                                    nPlanesToLoad++;
                                    System.out.println("Loaded slice plane " +
                                        slice + " at row: " + row + " col: " + col);
                                }
                            */
                            planeData = readBytePlaneData(ZDIR, getZCoor(slice));
                            floatData[n++] = getScaledValue(planeData[index]);
                        }
                        if (debug && nPlanesToLoad > 0)
                        {
                            System.out.println("Loaded " + nPlanesToLoad + " slice planes at row: " + row + " col: " + col);
                        }
                        break;
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getTraceValues", "row " + row + " col " + col + " sliceMin " + sliceMin +
                " sliceMax " + sliceMax + "values.length " + floatData.length, e);
            return false;
        }
    }

    /*
        public boolean getTraceValues(int dirNo, int row, int col, byte[] values)
        {
            try
            {
                int index;
                ByteBuffer planeData;
                switch (dirNo)
                {
                    case XDIR:
                        planeData = filesMapBlocks[XDIR].readByteBufferPlane(col);
                        index = row * nCroppedSlices;
                        for (int slice = 0; slice < nCroppedSlices; slice++, index++)
                        {
                            values[slice] = planeData.get(index);
                        }
                        return true;
                    case YDIR:
                        planeData = filesMapBlocks[YDIR].readByteBufferPlane(row);
                        index = col * nCroppedSlices;
                        for (int slice = 0; slice < nCroppedSlices; slice++, index++)
                        {
                            values[slice] = planeData.get(index);
                        }
                        return true;
                    default:
                        return false;
                }
            }
            catch (Exception e)
            {
                StsException.systemError(
                    "StsSeismicVolume.getValue() failed for row " + row + " col " +
                    col +
                    "values.length " +
                    values.length);
                return false;
            }
        }
    */
    /*
        public void initializeVelocityTrace(float rowF, float colF)
        {
            if (interpolatedTrace == null)interpolatedTrace = new InterpolatedTrace();
            interpolatedTrace.initialize(rowF, colF);
        }

        public float getVelocityTraceValue(int nSlice)
        {
            return interpolatedTrace.getValue(nSlice);
        }

        public class InterpolatedTrace
        {
            int row = -1, col = -1;
            float dx, dy;
            int index;
            float[] trace00, trace01, trace10, trace11;
            int type = 0;

            InterpolatedTrace()
            {
                setupReadRowFloatBlocks();
                trace00 = new float[nCroppedSlices];
                trace01 = new float[nCroppedSlices];
                trace10 = new float[nCroppedSlices];
                trace11 = new float[nCroppedSlices];
            }

            void initialize(float rowF, float colF)
            {
                row = (int)rowF;
                if (row == nRows - 1)
                {
                    row--;
                }
                dy = rowF - row;

                col = (int)colF;
                if (col == nCols - 1)
                {
                    col--;
                }
                dx = colF - col;

                if (dy == 0.0f)
                {
                    if (dx == 0.0f)
                    {
                        getTraceValues(row, col, 0, nCroppedSlices - 1, YDIR, false, trace00);
                        type = 0;
                    }
                    else
                    { // dy == 0.0f, point is along bottom row
                        getTraceValues(row, col, 0, nCroppedSlices - 1, YDIR, false, trace00);
                        getTraceValues(row, col + 1, 0, nCroppedSlices - 1, YDIR, false,
                                  trace01);
                        type = 1;
                    }
                }
                else
                { // dy != 0.0f
                    if (dx == 0.0f)
                    { // point is along left column
                        getTraceValues(row, col, 0, nCroppedSlices - 1, YDIR, false, trace00);
                        getTraceValues(row + 1, col, 0, nCroppedSlices - 1, YDIR, false,
                                  trace10);

                        type = 3;
                    }
                    else
                    { // point is in the middle
                        getTraceValues(row, col, 0, nCroppedSlices - 1, YDIR, false, trace00);
                        getTraceValues(row, col + 1, 0, nCroppedSlices - 1, YDIR, false,
                                  trace01);
                        getTraceValues(row + 1, col, 0, nCroppedSlices - 1, YDIR, false,
                                  trace10);
                        getTraceValues(row + 1, col + 1, 0, nCroppedSlices - 1, YDIR, false,
                                  trace11);
                        type = 4;
                    }
                }
            }

            float getValue(int nSlice)
            {
                nSlice = StsMath.minMax(nSlice, 0, nCroppedSlices - 1);
                if (type == 0)
                { // point is in lower-left corner
                    return trace00[nSlice];
                }
                else if (type == 1)
                { // point is along bottom row
                    return (1.0f - dx) * trace00[nSlice] + dx * trace01[nSlice];
                }
                else if (type == 3)
                { // point is along left column
                    return (1.0f - dy) * trace00[nSlice] + dy * trace10[nSlice];
                }
                else
                { // type == 4, point is in the middle
                    float value1 = (1.0f - dx) * trace00[nSlice] +
                        dx * trace01[nSlice];
                    float value2 = (1.0f - dx) * trace10[nSlice] +
                        dx * trace11[nSlice];
                    return (1.0f - dy) * value1 + dy * value2;
                }
            }
        }
    */
    public float getNearestFloatValue(float x, float y, float z)
    {
        int row = getNearestBoundedRowCoor(y);
        int col = getNearestBoundedColCoor(x);
        int plane = getNearestBoundedSliceCoor(z);
        return fileMapRowFloatBlocks.getRowCubeFloatValue(row, col, plane);
    }

    public float getTrilinearFloatValue(float x, float y, float z)
    {
        return getTrilinearFloatValue(x, y, z, false);
    }

    public float getTrilinearFloatValue(float x, float y, float z,
                                        boolean debug)
    {
        if (fileMapRowFloatBlocks == null)
        {
            return nullValue;
        }
        float rowF = getBoundedRowCoor(y);
        float colF = getBoundedColCoor(x);
        float sliceF = getBoundedSliceCoor(z);
        int row = (int) rowF;
        int col = (int) colF;
        int slice = (int) sliceF;
        float dRow = rowF - row;
        float dCol = colF - col;
        float dSlice = sliceF - slice;
        float value;
        if (dSlice == 0.0f)
        {
            value = getBilinearFloatValue(row, dRow, col, dCol, slice);
        }
        else
        {
            float topValue = getBilinearFloatValue(row, dRow, col, dCol, slice);
            float botValue = getBilinearFloatValue(row, dRow, col, dCol,
                slice + 1);
            value = topValue * (1.0f - dSlice) + botValue * dSlice;
        }
        if (debug)
        {
            System.out.println("getTrilinearFloatValue mainDebug. row " + row +
                " col " + col + " slice " + slice +
                " value " + value);
        }
        return value;
    }

    public float getBilinearFloatValue(int row, float dRow, int col, float dCol,
                                       int nSlice)
    {
        if (dRow == 0.0f)
        {
            if (dCol == 0.0f)
            {
                return fileMapRowFloatBlocks.getRowCubeFloatValue(row, col, nSlice);
            }
            else
            {
                float v00 = fileMapRowFloatBlocks.getRowCubeFloatValue(row, col,
                    nSlice);
                float v01 = fileMapRowFloatBlocks.getRowCubeFloatValue(row, col + 1,
                    nSlice);
                return v00 * (1.0f - dCol) + v01 * dCol;
            }
        }
        else
        {
            if (dCol == 0.0f)
            {
                float v00 = fileMapRowFloatBlocks.getRowCubeFloatValue(row, col,
                    nSlice);
                float v10 = fileMapRowFloatBlocks.getRowCubeFloatValue(row + 1, col,
                    nSlice);
                return v00 * (1.0f - dRow) + v10 * dRow;
            }
            else
            {
                float v00 = fileMapRowFloatBlocks.getRowCubeFloatValue(row, col,
                    nSlice);
                float v01 = fileMapRowFloatBlocks.getRowCubeFloatValue(row, col + 1,
                    nSlice);
                float v0 = v00 * (1.0f - dCol) + v01 * dCol;
                float v10 = fileMapRowFloatBlocks.getRowCubeFloatValue(row + 1, col,
                    nSlice);
                float v11 = fileMapRowFloatBlocks.getRowCubeFloatValue(row + 1,
                    col + 1, nSlice);
                float v1 = v10 * (1.0f - dCol) + v11 * dCol;
                return v0 * (1.0f - dRow) + v1 * dRow;
            }
        }
    }

    public byte getSliceBlockValue(int row, int col, float z)
    {
        int slice = Math.round(getSliceCoor(z));
        if (!isInsideRange(row, col, slice)) return nullByte;
        return filesMapBlocks[ZDIR].readByteValue(row, col, slice);
    }

    public byte getSliceBlockValue(int row, int col, int slice)
    {
        if (!isInsideRange(row, col, slice)) return nullByte;
        return filesMapBlocks[ZDIR].readByteValue(row, col, slice);
    }

    public byte getSliceBlockValue(float rowF, float colF, float z)
    {
        int row = Math.round(rowF);
        int col = Math.round(colF);
        int slice = Math.round(getSliceCoor(z));
        if (!isInsideRange(row, col, slice)) return nullByte;
        return filesMapBlocks[ZDIR].readByteValue(row, col, slice);
    }

    public float getRowFloatBlockValue(float rowF, float colF, float z)
    {
        int row = Math.round(rowF);
        int col = Math.round(colF);
        int slice = Math.round(getSliceCoor(z));
        return getRowFloatBlockValue(row, col, slice);
    }

    public float getRowFloatBlockValue(int row, int col, int slice)
    {
        if (!isInsideRange(row, col, slice)) return nullValue;
        return fileMapRowFloatBlocks.readFloatValue(row, col, slice);
    }

    /*
        public float getSliceByteBlockFloatValue(int row, int col, int slice)
        {
            if (!isInsideRange(row, col, slice))return nullValue;
            byte byteValue = filesMapBlocks[ZDIR].readByteValue(row, col, slice);
            if(byteValue == nullByte)
                return nullValue;
            else
                return StsMath.signedByteToFloatWithScale(byteValue, scale, scaleOffset);
        }
    */
    public void clearCache()
    {
        //		if(dataCubeMemory == null || dataCubes == null)return;
        //		dataCubeMemory.clearCache(dataCubes);
    }

    static public void close(StsSeismicVolume volume)
    {
        if (volume == null) return;
        volume.close();
    }

    public String getSegyFilename()
    {
        return segyFilename;
    }

    public void setSegyFilename(String segyFilename)
    {
        this.segyFilename = segyFilename;
    }

    public String getSegyDirectory()
    {
        return segyDirectory;
    }

    public void setSegyDirectory(String segyDirectory)
    {
        this.segyDirectory = segyDirectory;
    }

    public String getStsDirectory()
    {
        return stsDirectory;
    }

    public void setStsDirectory(String stsDirectory)
    {
        this.stsDirectory = stsDirectory;
    }

    public void setStemname(String stemname)
    {
        this.stemname = stemname;
    }

    public boolean getIsXLineCCW()
    {
        return isXLineCCW;
    }

    public void setIsXLineCCW(boolean isXLineCCW)
    {
        this.isXLineCCW = isXLineCCW;
    }

    public boolean getIsRegular()
    {
        return isRegular;
    }

    public void setIsRegular(boolean isRegular)
    {
        this.isRegular = isRegular;
    }

    //    public float getSampleSpacing() { return getZInc(); }
    //    public void setSampleSpacing(float sampleSpacing) { setZInc(sampleSpacing); }

    //    public void setLineAngle(float lineAngle) { setAngle(lineAngle); }
    //    public float getLineAngle() { return getAngle(); }

    public boolean getIsVisibleOnCursor()
    {
        return isVisible && getSeismicVolumeClass().getIsVisibleOnCursor();
    }

    //    public StsRotatedGridBoundingBox getBoundingBox() { return new StsRotatedGridBoundingBox(this); }

    public byte[] getCurrentCursorPlaneData()
    {
        return null;
    }

    public StsFieldBean[] getDisplayFields()
    {
        displayAttributeBean.setListItems(displayAttributes);
        return seismicDisplayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return seismicPropertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    /** volume has been removed from display; delete allocated planes, etc */
    public void deleteTransients()
    {
        clearCache();
    }

    public StsSurface constructSurface(String name, StsColor stsColor, byte type)
    {
        // note that angle is 0 because surface coordinates will be in rotated coordinate system of volume
        byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
        byte vUnits = currentModel.getProject().getDepthUnits();
        byte hUnits = currentModel.getProject().getXyUnits();
        return StsSurface.constructSurface(name, stsColor, type, this.getNCols(), getNRows(), xOrigin, yOrigin,
            getXInc(), getYInc(), getXMin(), getYMin(), angle, null, false, StsParameters.nullValue, zDomainByte, vUnits, hUnits, null);
    }

    public boolean canDisplayZDomain()
    {
        StsProject project = currentModel.getProject();
        return project.canDisplayZDomain(getZDomain()) || project.hasVelocityModel();
    }

    /** Delete the seismic colorlist. Can be called ONLY when GL context is active (during display operation). */
    /*
    private void deleteColorList(GL gl)
    {
     gl.deleteLists(colorListNum, 1);
     colorListNum = 0;
    }
    */
    // StsSurfaceDisplayable interface:

    /** Gets the display type (Cell or Grid); */
    //	public byte getDisplayType() { return StsSurfaceDisplayable.ORTHO_GRID_CENTERED; }
    /** This property can be isVisible on this surace. */
    public boolean isDisplayable()
    {
        return true;
    }

    /** Call for getting a cell-centered color */
    public Color getCellColor(int row, int col, int layer)
    {
        return Color.RED;
    }

    /** Call for getting a grid-centered color */
    public Color getGridColor(int row, int col, float z)
    {
        return getColor(row, col, z);
    }

    public Color[][] get2dColorArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset)
    {
        int nSurfaceRows = surfaceBoundingBox.nRows;
        int nSurfaceCols = surfaceBoundingBox.nCols;

        Color[][] colors = new Color[nSurfaceRows][nSurfaceCols];
        Color[] colorscaleColors = colorscale.getNewColorsInclTransparency();
        int value;
        int currentK = -1;
        byte[] planeData = null;

        if (timer == null)
        {
            timer = new StsTimer();
        }
        timer.start();

        int surfaceRowMin = surfaceBoundingBox.getNearestBoundedRowCoor(yMin);
        int surfaceRowMax = surfaceBoundingBox.getNearestBoundedRowCoor(yMax);
        int surfaceColMin = surfaceBoundingBox.getNearestBoundedColCoor(xMin);
        int surfaceColMax = surfaceBoundingBox.getNearestBoundedColCoor(xMax);

        Color nullColor = colorscaleColors[0];

        for (int surfaceRow = 0; surfaceRow < nSurfaceRows; surfaceRow++)
        {
            if (surfaceRow < surfaceRowMin || surfaceRow > surfaceRowMax)
            {
                for (int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
                {
                    colors[surfaceRow][surfaceCol] = nullColor;
                }
            }
            else
            {
                int row = surfaceRow - surfaceRowMin;
                for (int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
                {
                    if (surfaceCol < surfaceColMin || surfaceCol > surfaceColMax)
                    {
                        colors[surfaceRow][surfaceCol] = nullColor;
                    }
                    else
                    {
                        int col = surfaceCol - surfaceColMin;

                        float kF = (float) ((z[surfaceRow][surfaceCol] + zOffset - zMin) / zInc);
                        int k = (int) kF;
                        if (kF - k > 0.5f)
                        {
                            k++;
                        }
                        if (k != currentK)
                        {
                            planeData = readBytePlaneData(ZDIR, getZCoor(k));
                            currentK = k;
                        }
                        if (planeData == null)
                        {
                            value = 0;
                        }
                        else
                        {
                            byte signedByteValue = planeData[row * nCols + col];
                            value = StsMath.signedByteToUnsignedInt(signedByteValue);
                        }
                        colors[surfaceRow][surfaceCol] = colorscaleColors[value];
                    }
                }
            }
        }
        timer.stopPrint("get2dColorArray");
        return colors;
    }

    public byte[] getByteArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset)
    {
        int nSurfaceRows = surfaceBoundingBox.nRows;
        int nSurfaceCols = surfaceBoundingBox.nCols;

        byte[] bytes = new byte[nSurfaceRows * nSurfaceCols];

        int value;
        int currentK = -1;
        byte[] planeData = null;

        if (timer == null)
        {
            timer = new StsTimer();
        }
        timer.start();

        int surfaceRowMin = surfaceBoundingBox.getNearestBoundedRowCoor(yMin);
        int surfaceRowMax = surfaceBoundingBox.getNearestBoundedRowCoor(yMax);
        int surfaceColMin = surfaceBoundingBox.getNearestBoundedColCoor(xMin);
        int surfaceColMax = surfaceBoundingBox.getNearestBoundedColCoor(xMax);

        int n = 0;
        for (int surfaceRow = 0; surfaceRow < nSurfaceRows; surfaceRow++)
        {
            if (surfaceRow < surfaceRowMin || surfaceRow > surfaceRowMax)
            {
                continue;
            }

            int row = surfaceRow - surfaceRowMin;
            for (int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
            {
                if (surfaceCol < surfaceColMin || surfaceCol > surfaceColMax)
                {
                    continue;
                }

                int col = surfaceCol - surfaceColMin;

                float kF = (float) ((z[surfaceRow][surfaceCol] + zOffset - zMin) / zInc);
                int k = (int) kF;
                if (kF - k > 0.5f)
                {
                    k++;
                }
                if (k != currentK)
                {
                    planeData = readBytePlaneData(ZDIR, getZCoor(k));
                    currentK = k;
                }
                if (planeData != null)
                {
                    bytes[surfaceRow * nSurfaceCols + surfaceCol] = planeData[row * nCols + col];
                }
            }
        }
        timer.stopPrint("getByteArray");
        return bytes;
    }

    /**
     * Return number of rows and cols on isVisible 2D plane from 3D cube of nRows, nCols, nCroppedSlices.
     * Return nCroppedSlices,nRows for X-plane, nCroppedSlices,nCols for Y-plane, and nRows,nCols for Z-plane.
     */
    public int[] getCursorDisplayNRowCols(int dir)
    {
        switch (dir)
        {
            case XDIR:
                return new int[]
                    {
                        nSlices, nRows};
            case YDIR:
                return new int[]
                    {
                        nSlices, nCols};
            case ZDIR:
                return new int[]
                    {
                        nCols, nRows};
            default:
                return null;
        }
    }

    /**
     * Return number of rows and cols on data plane from 3D cube of nRows, nCols, nCroppedSlices.
     * Texture goes down first trace on vertical planes and across first row on horizontal plane.
     * Return nCols, nCroppedSlices for X-plane, nRows,nCroppedSlices for Y-plane, and nRows,nCols for Z-plane.
     */
    public int[] getCursorDataNRowCols(int dir)
    {
        switch (dir)
        {
            case XDIR:
                return new int[]
                    {
                        nRows, nSlices};
            case YDIR:
                return new int[]
                    {
                        nCols, nSlices};
            case ZDIR:
                return new int[]
                    {
                        nRows, nCols};
            default:
                return null;
        }
    }

    public byte[] readBytePlaneData(int dir, float dirCoordinate)
    {
        try
        {
            int nPlane = this.getCursorPlaneIndex(dir, dirCoordinate);
            if (nPlane == -1) return null;
            return filesMapBlocks[dir].readBytePlane(nPlane);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean canExport() { return true; }

    public boolean export()
    {
        try
        {

            StsSeismicExportPanel.createDialog(currentModel, this, "Export Seismic Volume", false);
            return true;
        }
        catch (Exception e)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to export " + getName());
            return false;
        }
    }

    public float[] readRowPlaneFloatData(int nPlane)
    {
        return fileMapRowFloatBlocks.readFloatPlane(nPlane);
    }

    public float[] readRowPlaneFloatData(float y)
    {
        int nPlane = getCursorPlaneIndex(YDIR, y);
        return fileMapRowFloatBlocks.readFloatPlane(nPlane);
    }

    public boolean readRowPlaneFloatData(float y, float[] floats)
    {
        int nPlane = getCursorPlaneIndex(YDIR, y);
        return fileMapRowFloatBlocks.readFloatPlane(nPlane, floats);
    }

    public byte[] readRowPlaneByteData(int nPlane)
    {
        return readBytePlaneData(YDIR, nPlane);
    }

    /*
        public ByteBuffer readRowPlaneByteBufferData(int nPlane)
        {
            return filesMapBlocks[YDIR].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY);
        }

        public FloatBuffer readRowPlaneFloatBufferData(int nPlane)
        {
            return fileMapRowFloatBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY).asFloatBuffer();
        }
    */
    public ByteBuffer readByteBufferPlane(int dir, float dirCoordinate)
    {
        try
        {
            if (dir == ZDIR && currentAttribute != null && currentAttribute != nullAttribute)
                return currentAttribute.getByteBuffer();
            else
            {
                int nPlane = this.getCursorPlaneIndex(dir, dirCoordinate);
                ByteBuffer byteBuffer = readByteBufferPlane(dir, nPlane);
                return checkInterpolateUserNull(dir, nPlane, byteBuffer);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public ByteBuffer readByteBufferPlane(int dir, int nPlane)
    {
        try
        {
            if (nPlane == -1) return null;
            if (filesMapBlocks == null || filesMapBlocks[dir] == null) return null;
            return filesMapBlocks[dir].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
            return null;
        }
    }

    private ByteBuffer checkInterpolateUserNull(int dir, int nPlane, ByteBuffer byteBuffer)
    {
        boolean fillPlaneNulls = getSeismicVolumeClass().getFillPlaneNulls();
        if (!fillPlaneNulls) return byteBuffer;
        StsVolumeInterpolation volumeInterpolator = StsVolumeInterpolation.getInstance(this);
        return volumeInterpolator.interpolatePlane(dir, nPlane, byteBuffer);
    }

    /*
        public ByteBuffer readByteBufferPlane(int dirNo, float dirCoordinate)
        {
            try
            {
                int nPlane = this.getCursorPlaneIndex(dirNo, dirCoordinate);
                if(nPlane == -1) return null;
                if(filesMapBlocks == null || filesMapBlocks[dirNo] == null) return null;
                MappedByteBuffer blockBuffer = filesMapBlocks[dirNo].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY);
                return blockBuffer.slice();
            }
            catch (Exception e)
            {
                StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
                return null;
            }
        }
    */

    // remove planes not between rowKeepMin and rowKeepMax as long as they are not between excludeRowMin and excludeRowMax.
    // then read planes between rowKeepMin and rowKeepMax to fill in any that are missing
    /*
        public void checkAllocateRowPlanes(int rowKeepMin, int rowKeepMax, int excludeRowMin, int excludeRowMax)
        {
            if (floatRowFilename != null)
            {
                try
                {
                    if (rowFloatDataCube == null)
                    {
                        rowFloatDataCube = new StsFloatDataCube(stsDirectory, floatRowFilename, "r", YDIR, dataCubeMemory, nRows, nCols, nCroppedSlices);
                    }
                    rowFloatDataCube.removeRowPlanes(rowKeepMin, rowKeepMax, excludeRowMin, excludeRowMax);
                    for (int row = rowKeepMin; row <= rowKeepMax; row++)
                    {
                        rowFloatDataCube.readPlaneData(row);
                    }
                    return;
                }
                catch (FileNotFoundException fnfe)
                {
                    floatRowFilename = null;
                }
            }
            dataCubes[YDIR].removeRowPlanes(rowKeepMin, rowKeepMax, excludeRowMin, excludeRowMax);
            for (int row = rowKeepMin; row <= rowKeepMax; row++)
            {
                readPlaneData(YDIR, row);
            }
            if (mainDebug)
            {
                System.out.println("Planes in memory: " + dataCubeMemory.nPlanesInMemory);
            }
        }
    */

    public void setCheckWriteCursorPlaneFlag(boolean calculateDataPlane)
    {
        this.calculateDataPlane = calculateDataPlane;
        //        System.out.println("Set calculateDataPlane: " + calculateDataPlane);
    }

    public byte[] readBytePlaneData(int dir, int nPlane)
    {
        return filesMapBlocks[dir].readBytePlane(nPlane);
    }

    private Object editPlaneData(int dir, int nPlane, byte[] data)
    {
        switch (type)
        {
            // convert from avg velocity to instantaneous velocity
            case StsParameters.SAMPLE_TYPE_VEL_AVG:
            {
                float scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
                float scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, dataMin);
                if (dir == ZDIR) return (Object) data;
                int nTraces = this.getNCursorRows(dir);
                byte[] outData = new byte[data.length];
                int n = 0;
                for (int t = 0; t < nTraces; t++)
                {
                    outData[n] = data[n];
                    float v1 = getScaledValue(data[n]);
                    n++;
                    for (int s = 1; s < nSlices; s++, n++)
                    {
                        float v0 = v1;
                        v1 = getScaledValue(data[n]);
                        float iv = s * v1 - (s - 1) * v0;
                        outData[n] = StsMath.floatToSignedByte254WithScale(iv, scale, scaleOffset);
                    }
                }
                return (Object) outData;
            }
            default:
                return (Object) data;
        }
    }

    public void setCalculateDataMethod(StsMethod method)
    {
        calculateDataPlaneMethod = method;
    }
    /*
        final public boolean planeInMemory(int dirNo, int nPlane)
        {
            return dataCubes[dirNo].getPlaneData(nPlane) != null;
        }

        final public boolean floatPlaneInMemory(int nPlane)
        {
            return rowFloatDataCube.planes[nPlane] != null;
        }

        public boolean traceInMemory(int sliceMin, int sliceMax)
        {
            Object[] planes = dataCubes[ZDIR].planes;
            for (int n = sliceMin; n <= sliceMax; n++)
            {
                if (planes[n] == null)
                {
                    return false;
                }
            }
            return true;
        }
    */
    /*
    public void deleteCache()
    {
     dataCubeMemory.deleteCache();
    }
    */

    public final boolean isByteValueNull(byte byteValue)
    {
        return byteValue == -1;
    }

    /** in byte range 0 to 254, the zero crossing is at this value */
    public float getUnsignedByteAverage()
    {
        return 127;
    }
    /*
        public synchronized void addItemListener(ItemListener listener)
        {
            if(itemListeners == null) itemListeners = new StsItemListeners();
            itemListeners.add(listener);
        }

        public synchronized void removeItemListener(ItemListener listener)
        {
            if(itemListeners == null) return;
            itemListeners.remove(listener);
        }

        protected void fireItemStateChanged(ItemEvent e)
        {
            if (itemListeners == null) return;
            itemListeners.fireItemStateChanged(e);
        }
    */

    //    public String getSuperVolumeName() { return superVolumeName; }
    //    public void setSuperVolumeName(String superVolumeName) { this.superVolumeName = superVolumeName; }

    /* get or create surfaceTexture for this surface. */

    public StsSurfaceTexture getSurfaceTexture(StsSurface surface)
    {
        if (surfaceTextures != null)
        {
            for (StsSurfaceTexture surfaceTexture : surfaceTextures)
            {
                if (surfaceTexture.surface == surface)
                    return surfaceTexture;
            }
        }
        SeismicSurfaceTexture surfaceTexture = new SeismicSurfaceTexture(surface, this);
        surfaceTextures = (StsSurfaceTexture[]) StsMath.arrayAddElement(surfaceTextures, surfaceTexture);
        return surfaceTexture;
    }

    class SeismicSurfaceTexture extends StsSurfaceTexture
    {
        StsSeismicVolume volume;

        public SeismicSurfaceTexture(StsSurface surface, StsSeismicVolume volume)
        {
            super(surface);
            this.volume = volume;
            initializeColorscaleActionListener();
        }

        public byte[] getTextureData()
        {
            byte[] textureData = null;
            float[][] pointsZ = surface.getPointsZ();
            byte[][] pointsNull = surface.getPointsNull();
            if (surface.sameAs(volume))
            {
                textureData = new byte[nRows * nCols];
                int n = 0;
                for (int row = 0; row < nRows; row++)
                {
                    for (int col = 0; col < nCols; col++)
                    {
                        if (pointsNull[row][col] != StsParameters.SURF_PNT)
                        {
                            textureData[n++] = -1;
                        }
                        else
                        {
                            float z = pointsZ[row][col] + surface.getOffset();
                            textureData[n++] = getSliceBlockValue(row, col, z);
                            //						textureData[n++] = getSliceBlockValue(row, col, pointsZ[row][col]);
                        }
                    }
                }
            }
            else
            {
                int nRows = surface.getNRows();
                int nCols = surface.getNCols();

                textureData = new byte[nRows * nCols];

                float seismicColStartF = getColCoor(surface.getXMin());
                float seismicRowStartF = getRowCoor(surface.getYMin());
                float seismicRowInc = surface.getYInc() / yInc;
                float seismicColInc = surface.getXInc() / xInc;

                int n = 0;
                float seismicRowF = seismicRowStartF;
                for (int row = 0; row < nRows; row++)
                {
                    float seismicColF = seismicColStartF;
                    for (int col = 0; col < nCols; col++)
                    {
                        if (pointsNull[row][col] != StsParameters.SURF_PNT)
                        {
                            textureData[n++] = -1;
                        }
                        else
                        {
                            float z = pointsZ[row][col] + surface.getOffset();
                            textureData[n++] = getSliceBlockValue(seismicRowF, seismicColF, z);
                            //						textureData[n++] = getSliceBlockValue(seismicRowF, seismicColF, pointsZ[row][col]);

                        }
                        seismicColF += seismicColInc;
                    }
                    seismicRowF += seismicRowInc;
                }
            }
            return textureData;
        }

        public float[] getHistogram()
        {
            return volume.getHistogram();
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            return volume.getColorDisplayListNum(gl, nullsFilled);
        }

        public boolean isDisplayable() { return true; }

        public String getName() { return volume.getName(); }

        public StsColorscale getColorscale() { return volume.getColorscale(); }

        public float getDataMin()
        {
            return dataMin;
        }
        public float getDataMax()
        {
            return dataMax;
        }
    }

    public byte[] getSeismicCurtainData(StsGridPoint[] gridPoints)
    {
        if (gridPoints == null)
        {
            return null;
        }
        int nPoints = gridPoints.length;
        int nSamples = getNSlices();
        StsGridPoint gridPoint;
        byte[] data = new byte[nPoints * nSamples];
        byte[] trace = new byte[nSamples];
        byte[] transparentTrace = new byte[nSamples];
        /*
                try
                {
                    if (!constructByteBlocks())
                    {
                        return null;
                    }
                }
                catch (Exception e)
                {
                    return null;
                }
        */
        for (int n = 0; n < nSamples; n++)
        {
            transparentTrace[n] = -1;
        }

        int row1 = gridPoints[0].row;
        int col1 = gridPoints[0].col;
        boolean traceOK;
        int i = 0;
        for (int n = 1; n < nPoints; n++)
        {
            int row0 = row1;
            int col0 = col1;
            row1 = gridPoints[n].row;
            col1 = gridPoints[n].col;
            int row = Math.min(row0, row1);
            int col = Math.min(col0, col1);
            traceOK = filesMapBlocks[YDIR].getRowCubeByteValues(row, col, 0, nSlices - 1, trace);
            if (traceOK)
            {
                System.arraycopy(trace, 0, data, i, nSamples);
            }
            else
            {
                System.arraycopy(transparentTrace, 0, data, i, nSamples);
            }
            i += nSamples;
        }
        System.arraycopy(transparentTrace, 0, data, i, nSamples);

        //		deleteRowByteBlocks();
        return data;
    }

    public float computeFloatFromByte(byte byteValue)
    {
        return StsMath.signedByteToFloat(byteValue, dataMin, dataMax);
    }

    public void display(StsGLPanel glPanel)
    {
        GL gl = glPanel.getGL();

        if (!isVisible) return;

        if (!currentModel.getProject().canDisplayZDomain(getZDomain())) return;

        byte projectZDomain = currentModel.getProject().getZDomain();
        if (projectZDomain != zDomainDisplayed)
            zDomainDisplayed = projectZDomain;

        if (colorscale != null && colorscale.canDrawVoxels()) drawIsoSurfaces(gl);
    }

    public void deleteDisplayList(GL gl)
    {
        if (displayListNum > 0)
        {
            gl.glDeleteLists(displayListNum, 1);
            displayListNum = 0;
        }
    }

    void drawIsoSurfaces(GL gl)
    {
        if (!checkSetupByteOrFloatRowBlocks()) return;
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        if (fileMapRowFloatBlocks != null)
            drawFloatIsoSurfaces(gl);
        else
            drawByteIsoSurfaces(gl);
    }

    public boolean checkSetupByteOrFloatRowBlocks()
    {
        isDataFloat = rowFloatFilename != null;
        if (isDataFloat)
        {
            if (fileMapRowFloatBlocks != null) return true;
            return setupReadRowFloatBlocks();
        }
        return true;  // using byteBlocks
        /*
            else
            {
                if (fileMapRowByteBlocks != null)return true;
                return setupRowByteBlocks();
            }
        */
    }

    private void drawFloatIsoSurfaces(GL gl)
    {
    }

    public Object getDisplayAttribute() { return currentAttribute; }

    public void setDisplayAttribute(Object attribute)
    {
        if (currentAttribute == attribute) return;
        currentAttribute = (Attribute) attribute;
        basemapDisplayChanged = true;
        currentModel.getCursor3d().clearTextureDisplays(ZDIR);
        currentModel.viewObjectChanged(this, this);
	}

    public boolean displayingBasemap()
    {
        return currentAttribute != nullAttribute;
    }

    public boolean hasAttributeChanged()
    {
        return basemapDisplayChanged;
    }

    public boolean displayingAttribute()
    {
        return currentAttribute != nullAttribute;
    }

    public boolean basemapTextureChanged()
    {
        if (basemapDisplayChanged)
        {
            basemapDisplayChanged = false;
            return true;
        }
        else
            return false;
    }

    public void initDisplayAttributes()
    {
        currentAttribute = nullAttribute;
        displayAttributes = (Attribute[]) StsMath.arrayAddElement(displayAttributes, nullAttribute);
        if (attributeNames == null) return;
        for (int n = 0; n < attributeNames.length; n++)
            displayAttributes = (Attribute[]) StsMath.arrayAddElement(displayAttributes, new Attribute(attributeNames[n]));
    }

    /*
// test version:  draws points along trace which are > isoValue
    private void drawByteIsoSurfaces(GL gl)
    {
     try
     {
      gl.begin(GL.GL_POINTS);
      gl.pointSize(1.0f);
      //			gl.glDisable(GL.GL_DEPTH_TEST);
      StsColor.setGLJavaColor(gl, getColor(isoValue));
      int[] rowPlane = new int[nCols * nCroppedSlices];

      float y = yMin;
      for (int row = 0; row < nRows; row++, y += yInc)
      {
    fileMapRowByteBlocks.getPlaneValues(row, rowPlane);
    int n = 0;
    float x = xMin;
    for (int col = 0; col < nCols; col++, x += xInc)
    {
     float z = zMin;
     int v1 = rowPlane[n++];
     int nPoints = 0;
     for (int slice = 1; slice < nCroppedSlices; slice++, z += zInc)
     {
      int v0 = v1;
      v1 = rowPlane[n++];
      if (v0 >= isoValue && v1 <= isoValue || v0 <= isoValue && v1 >= isoValue)
      {
    float f;
    if(v0 == v1)
     f = 0.5f;
    else
     f = (isoValue - v0) / (v1 - v0);
    gl.glVertex3f(x, y, z + f * zInc);
    nPoints++;
      }
     }
     if(row < 5) System.out.println("row " + row + " col " + col + " nPoints " + nPoints);
    }
      }
     }
     catch (Exception e)
     {
      StsException.outputException("StsSeismicVolume.drawByteIsoSurface() failed.", e, StsException.WARNING);
     }
     finally
     {
      gl.glEnd();
      //			gl.glEnable(GL.GL_DEPTH_TEST);
     }
    }
    */
    /*
// test version: draws point for top or bottom face if face average value is > isoValue
    private void drawByteIsoSurfaces(GL gl)
    {
    try
     {
      gl.begin(GL.GL_POINTS);
      gl.pointSize(1.0f);
      //			gl.glDisable(GL.GL_DEPTH_TEST);
      StsColor.setGLJavaColor(gl, getColor(isoValue));
      int[] rowPlane0 = new int[nCols * nCroppedSlices];
      int[] rowPlane1 = new int[nCols * nCroppedSlices];
      int[] savePlane;

      fileMapRowByteBlocks.getPlaneValues(0, rowPlane1);
      float yCenter = yMin;
      for (int row = 0; row < nRows - 1; row++, yCenter += yInc)
      {
    savePlane = rowPlane0;
    rowPlane0 = rowPlane1;
    rowPlane1 = savePlane;
    fileMapRowByteBlocks.getPlaneValues(row, rowPlane1);
    int n = 0;
    float xCenter = xMin;
    for (int col = 0; col < nCols - 1; col++, xCenter += xInc)
    {
     float zCenter = zMin;
     for (int slice = 0; slice < nCroppedSlices - 1; slice++, n++, zCenter += zInc)
     {
      float[] gradient = null;
      int v0 = rowPlane0[n];
      int vr = rowPlane1[n];
      int vc = rowPlane0[n + nCols];
      int vs = rowPlane0[n + 1];
      int v = (v0 + vr + vs + vc)/4;
      if (v >= isoValue)
      {
    gradient = new float[] { vc - v0, vr - v0, vs - v0 };
    gl.glNormal3fv(gradient);
    gl.glVertex3f(xCenter, yCenter, zCenter);
      }
     }
     n++;
    }
      }
     }
     catch (Exception e)
     {
      StsException.outputException("StsSeismicVolume.drawByteIsoSurface() failed.", e, StsException.WARNING);
     }
     finally
     {
      gl.glEnd();
      //			gl.glEnable(GL.GL_DEPTH_TEST);
     }
    }
    */
    /*
// test version: draws point along each  edge where isoValue occurs
    private void drawByteIsoSurfaces(GL gl)
    {
    try
     {
      gl.begin(GL.GL_POINTS);
      gl.pointSize(1.0f);
      //			gl.glDisable(GL.GL_DEPTH_TEST);
      StsColor.setGLJavaColor(gl, getColor(isoValue));
      int[] rowPlane0 = new int[nCols * nCroppedSlices];
      int[] rowPlane1 = new int[nCols * nCroppedSlices];
      int[] savePlane;

      fileMapRowByteBlocks.getPlaneValues(0, rowPlane1);
      float yCenter = yMin;
      for (int row = 0; row < nRows - 1; row++, yCenter += yInc)
      {
    savePlane = rowPlane0;
    rowPlane0 = rowPlane1;
    rowPlane1 = savePlane;
    fileMapRowByteBlocks.getPlaneValues(row, rowPlane1);
    int n = 0;
    float xCenter = xMin;
    for (int col = 0; col < nCols - 1; col++, xCenter += xInc)
    {
     float zCenter = zMin;
     for (int slice = 0; slice < nCroppedSlices - 1; slice++, n++, zCenter += zInc)
     {
      float[] gradient = null;
      int v0 = rowPlane0[n];
      int vr = rowPlane1[n];
      int vc = rowPlane0[n + nCols];
      int vs = rowPlane0[n + 1];
      if (v0 < isoValue)
      {
    if (vr > isoValue)
    {
     if (gradient == null)
      gradient = new float[]
    {
    vc - v0, vr - v0, vs - v0};
     gl.glNormal3fv(gradient);
     float f = (isoValue - v0) / (vr - v0);
     gl.glVertex3f(xCenter, yCenter + yInc / 2, zCenter);
    }
    if (vc > isoValue)
    {
     if (gradient == null)
      gradient = new float[]
    {
    vc - v0, vr - v0, vs - v0};
     gl.glNormal3fv(gradient);
     float f = (isoValue - v0) / (vc - v0);
     gl.glVertex3f(xCenter + xInc / 2, yCenter, zCenter);
    }
    if (vs > isoValue)
    {
     if (gradient == null)
      gradient = new float[]
    {
    vc - v0, vr - v0, vs - v0};
     gl.glNormal3fv(gradient);
     float f = (isoValue - v0) / (vs - v0);
     gl.glVertex3f(xCenter, yCenter, zCenter + zInc / 2);
    }
      }
     }
     n++;
    }
      }
     }
     catch (Exception e)
     {
      StsException.outputException("StsSeismicVolume.drawByteIsoSurface() failed.", e, StsException.WARNING);
     }
     finally
     {
      gl.glEnd();
      //			gl.glEnable(GL.GL_DEPTH_TEST);
     }
    }
    */

    void drawByteIsoSurfaces(GL gl)
    {
        if (voxelVertexArraySets == null)
        {
            try
            {
                constructByteIsoSurfaces();
            }
            catch (java.lang.OutOfMemoryError E)
            {
                System.out.println("Out of memory creating voxels!");
            }
        }
        if (voxelVertexArraySets == null) return;
        int nVertices = 0;
        for (int n = 0; n < voxelVertexArraySets.length; n++)
            nVertices += voxelVertexArraySets[n].getNVertices();
        if (timer == null) timer = new StsTimer();
        timer.start();
        for (int n = 0; n < voxelVertexArraySets.length; n++)
        //		for (int n = 1; n < voxelVertexArraySets.length; n += 2)
        {
            StsColor.setGLJavaColor(gl, getColor(voxelColorIndexes[n]));
            voxelVertexArraySets[n].draw(gl);
        }
        timer.stopPrint("nVertices drawn " + nVertices);
    }

    static int nGridsInSet = 0;
    static int nCreated = 0;
    static int nTotal = 0;
    static int[] indices = null;
    static int[] count = null;
    static int nIndicesAllocated = 0;

    private void constructByteIsoSurfaces()
    {
        createVoxelIndexVector();
        //        boolean[] isVoxelInsideVector = getVoxelIsInsideVector();
        //        if(isVoxelInsideVector == null) return;
        isoSurfaceID = 0;

        int[] rowPlane0 = new int[nCols * nSlices];
        int[] rowPlane1 = new int[nCols * nSlices];
        int[] savePlane;

        VoxelTracePoints[] voxelTracePointsRow0 = new VoxelTracePoints[nCols];
        VoxelTracePoints[] voxelTracePointsRow1 = new VoxelTracePoints[nCols];
        VoxelTracePoints[] saveVoxelTracePointsRow;

        isoValueGridsRow0 = new IsoValueGridSet(100);
        isoValueGridsRow1 = new IsoValueGridSet(100);
        IsoValueGridSet saveIsoValueGridsRow;

        isoValueGrids = new IsoValueGridSet(100);

        // progress

        StsStatusArea statusArea = currentModel.win3d.statusArea;
        statusArea.addProgress();
        statusArea.setText("Creating Voxels");
        statusArea.setMinimum(0);
        statusArea.setMaximum(100);
        statusArea.addAll();

        // construct pairs (blob entry & exit) along each trace
        float y = getYMin();
        getRowPlaneValues(y, rowPlane1);
        int n = 0;
        for (int col = 0; col < nCols; col++, n += nSlices)
            voxelTracePointsRow1[col] = new VoxelTracePoints(0, col, rowPlane1, n, nSlices);

        statusArea.setProgress(10);
        checkAddTraceBlobs(voxelTracePointsRow1[0], null, null);
        for (int col = 1; col < nCols; col++)
            checkAddTraceBlobs(voxelTracePointsRow1[col], null, voxelTracePointsRow1[col - 1]);
        if (debug) System.out.println("Completed row 0. ");
        y += yInc;
        statusArea.setProgress(20);
        for (int row = 1; row < nRows; row++, y += yInc)
        {
            savePlane = rowPlane0;
            rowPlane0 = rowPlane1;
            rowPlane1 = savePlane;
            statusArea.setProgress(40.f / nRows * row);
            getRowPlaneValues(y, rowPlane1);
            //			checkFreeRowBlocks(0, row);
            saveIsoValueGridsRow = isoValueGridsRow0;
            isoValueGridsRow0 = isoValueGridsRow1;
            isoValueGridsRow1 = saveIsoValueGridsRow;

            saveVoxelTracePointsRow = voxelTracePointsRow0;
            voxelTracePointsRow0 = voxelTracePointsRow1;
            voxelTracePointsRow1 = saveVoxelTracePointsRow;

            n = 0;
            for (int col = 0; col < nCols; col++, n += nSlices)
                voxelTracePointsRow1[col] = new VoxelTracePoints(row, col, rowPlane1, n, nSlices);

            checkAddTraceBlobs(voxelTracePointsRow1[0], voxelTracePointsRow0[0], null);
            for (int col = 1; col < nCols; col++)
            {
                checkAddTraceBlobs(voxelTracePointsRow1[col], voxelTracePointsRow0[col], voxelTracePointsRow1[col - 1]);
            }
            isoValueGridsRow0.processRowGrids(row);
            int rowM1 = row - 1;
            if (debug) System.out.println("Processed row " + rowM1 + " memory used " + getMemoryUsed());
            if (row % 100 == 0) isoValueGrids.processDisconnectedGrids(row);
        }
        isoValueGridsRow1.addToGrids(isoValueGrids); // add last row to volumeGrids
        isoValueGrids.initialize(); // classInitialize each isoValueGrid and remove volumeGridSet
        for (n = 0; n < voxelVertexArraySets.length; n++)
            voxelVertexArraySets[n].trim();

        clearVoxelAllocations();
        statusArea.removeProgress();
    }

    private void clearVoxelAllocations()
    {
        indices = null;
        count = null;
        nIndicesAllocated = 0;
    }

    private boolean getRowPlaneValues(float y, int[] data)
    {
        byte[] byteData;

        int nPlane = this.getCursorPlaneIndex(YDIR, y);
        if (nPlane == -1) return false;
        return filesMapBlocks[YDIR].getRowCubePlaneValues(nPlane, data);
    }

    public FloatBuffer getRowPlaneFloatBuffer(int row)
    {
        return fileMapRowFloatBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY).asFloatBuffer();
    }

   public FloatBuffer getRowPlaneFloatBuffer(int row, int colMin)
    {
        try
        {
            FloatBuffer rowFloatBuffer =  fileMapRowFloatBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY).asFloatBuffer();
            rowFloatBuffer.position(colMin*nSlices);
            return rowFloatBuffer;
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getRowPlaneFloatBuffer", "Failed to get floatBuffer at row " + row + " col " + colMin);
            return null;
        }
    }

    public ByteBuffer getRowPlaneByteBuffer(int row)
    {
        return fileMapRowFloatBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY);
    }

    /*
        private void checkFreeRowBlocks(int minRow, int maxRow)
        {
            if(fileMapRowByteBlocks != null) fileMapRowByteBlocks.checkFreeBlocks(minRow, maxRow);
        }
    */
    private long getMemoryUsed()
    {
        long memoryUsed = 0;
        for (int n = 0; n < voxelVertexArraySets.length; n++)
            memoryUsed += voxelVertexArraySets[n].getMemoryUsed();
        return memoryUsed;
    }

    /*
    private boolean[] getVoxelIsInsideVector()
    {
     StsColorscale.VoxelKey[] voxelKeys = colorscale.getVoxelKeys();
     if(voxelKeys == null || voxelKeys.length == 0) return null;
     boolean[] isInside = new boolean[256];
     for(int i = 0; i < 256; i++)
      isInside[i] = false;
     for(int n = 0; n < voxelKeys.length; n++)
     {
      int min = voxelKeys[n].getMin();
      int max = voxelKeys[n].getMax();
      for(int i = min; i <= max; i++)
    isInside[i] = true;
     }
     return isInside;
    }
    */
    /**
     * For each voxelKey range, create 3 blob index numbers and a corresponding color:
     * top, middle, and bottom. So if we have 2 ranges, we will have 6 indices and colors.
     */

    private void createVoxelIndexVector()
    {
        int nKeys;
        voxelKeys = colorscale.getVoxelKeys();
        if (voxelKeys == null || (nKeys = voxelKeys.length) == 0) return;
        voxelKeyIndexes = new int[256];
        for (int i = 0; i < 256; i++)
            voxelKeyIndexes[i] = -1;
        voxelVertexArraySets = new VoxelVertexArrays[nKeys * 2];
        voxelColorIndexes = new int[nKeys * 2];
        int nn = 0;
        for (int n = 0; n < nKeys; n++)
        {
            int min = voxelKeys[n].getMin();
            int max = voxelKeys[n].getMax();
            for (int i = min; i <= max; i++)
                voxelKeyIndexes[i] = n;
            voxelVertexArraySets[nn] = new VoxelVertexArrays();
            voxelColorIndexes[nn++] = min;
            voxelVertexArraySets[nn] = new VoxelVertexArrays();
            voxelColorIndexes[nn++] = max;
        }
    }

    /** add top and bot points for the newTrace, getting IDs from prevTrace points */
    private void checkAddTraceBlobs(VoxelTracePoints newVoxelTracePoints, VoxelTracePoints prevRowVoxelTracePoints,
                                    VoxelTracePoints prevColVoxelTracePoints)
    {
        if (prevRowVoxelTracePoints != null) this.addCorrelatedTraceBlobs(newVoxelTracePoints.traceBlobs,
            prevRowVoxelTracePoints.traceBlobs);
        if (prevColVoxelTracePoints != null) addCorrelatedTraceBlobs(newVoxelTracePoints.traceBlobs,
            prevColVoxelTracePoints.traceBlobs);
        addUncorrelatedTraceBlobs(newVoxelTracePoints.traceBlobs);
    }

    private void addCorrelatedTraceBlobs(TraceBlob[] newBlobs, TraceBlob[] prevBlobs)
    {
        int nNewBlobs = newBlobs.length;
        if (nNewBlobs == 0) return;

        int nPrevBlobs = prevBlobs.length;
        if (nPrevBlobs == 0) return;

        int no = 0;
        TraceBlob prevBlob = prevBlobs[no++];
        int prevBlobBotIndex = prevBlob.botIndex;

        for (int nn = 0; nn < nNewBlobs; nn++)
        {
            TraceBlob newBlob = newBlobs[nn];
            int newBlobTopIndex = newBlob.topIndex;
            int newBlobBotIndex = newBlob.botIndex;

            while (no < nPrevBlobs && prevBlobBotIndex < newBlobTopIndex)
            {
                prevBlob = prevBlobs[no++];
                prevBlobBotIndex = prevBlob.botIndex;
            }
            int prevBlobTopIndex = prevBlob.topIndex;
            if (prevBlobBotIndex >= newBlobTopIndex - 1 && prevBlobTopIndex <= newBlobBotIndex + 1)
            {
                if (prevBlob.voxelKeyIndex == newBlob.voxelKeyIndex && prevBlob.isTopMax == newBlob.isTopMax && prevBlob.isBotMax == newBlob.isBotMax)
                    addTraceBlob(newBlob, prevBlob);
            }
        }
    }

    private void addUncorrelatedTraceBlobs(TraceBlob[] newBlobs)
    {
        int nNewPoints = newBlobs.length;
        if (nNewPoints == 0) return;
        for (int nn = 0; nn < nNewPoints; nn++)
        {
            TraceBlob newBlob = newBlobs[nn];
            if (newBlob.id == -1) addTraceBlob(newBlob, null);
        }
    }

    private void addTraceBlob(TraceBlob newBlob, TraceBlob prevBlob)
    {
        VoxelIsoValueGrid voxelIsoValueGrid;
        int id;

        if (prevBlob == null)
        {
            id = isoSurfaceID++;
            voxelIsoValueGrid = this.isoValueGridsRow1.addGrid(id, newBlob.voxelKeyIndex);
        }
        else
        {
            id = prevBlob.id;
            voxelIsoValueGrid = getVoxelGrid(id);
            if (voxelIsoValueGrid == null)
            {
                StsException.systemError("StsSeismicVolume.addTracePoint() failed. Couldn't get voxelIsoValueGrid for id " +
                    id + " from prevPoint at row " + prevBlob.row + " col " + prevBlob.col);
                return;
            }
        }
        newBlob.id = id;
        voxelIsoValueGrid.addTraceBlob(newBlob);
    }

    VoxelIsoValueGrid getVoxelGrid(int id)
    {
        VoxelIsoValueGrid voxelIsoValueGrid;

        voxelIsoValueGrid = isoValueGridsRow1.getGrid(id);
        if (voxelIsoValueGrid != null) return voxelIsoValueGrid;
        voxelIsoValueGrid = isoValueGridsRow0.getGrid(id);
        if (voxelIsoValueGrid != null) return voxelIsoValueGrid;
        voxelIsoValueGrid = isoValueGrids.getGrid(id);
        if (voxelIsoValueGrid != null) return voxelIsoValueGrid;
        return null;
    }

    class IsoValueGridSet
    {
        int number = 0;
        int length;
        VoxelIsoValueGrid[] voxelGrids = new VoxelIsoValueGrid[100];

        IsoValueGridSet(int length)
        {
            this.length = length;
            voxelGrids = new VoxelIsoValueGrid[length];
        }

        VoxelIsoValueGrid addGrid(int id, int voxelKeyIndex)
        {
            VoxelIsoValueGrid newVoxelGrid = new VoxelIsoValueGrid(id, voxelKeyIndex);
            addGrid(newVoxelGrid);
            return newVoxelGrid;
        }

        void addGrid(VoxelIsoValueGrid voxelIsoValueGrid)
        {
            if (number == length)
            {
                int newLength = length * 2;
                VoxelIsoValueGrid[] newVoxelGrids = new VoxelIsoValueGrid[newLength];
                System.arraycopy(voxelGrids, 0, newVoxelGrids, 0, length);
                voxelGrids = newVoxelGrids;
                length = newLength;
            }
            voxelGrids[number++] = voxelIsoValueGrid;
        }

        VoxelIsoValueGrid getGrid(int target)
        {
            int high = number, low = -1, probe;
            while (high - low > 1)
            {
                probe = (high + low) / 2;
                int id = voxelGrids[probe].id;
                if (id > target)
                    high = probe;
                else
                    low = probe;
            }
            if (low == -1 || voxelGrids[low].id != target)
                return null;
            else
                return voxelGrids[low];
        }

        /**
         * This IsoValueGridSet is for the row before row just finished.
         * If a grid in this row is disconnected (doesn't have entry in row just finished),
         * then either delete it if it is a small point, or write it out (disconnected so we won't refer to it again).
         * If not disconnected, add it to volume set.
         */
        void processRowGrids(int row)
        {
            for (int n = 0; n < number; n++)
            {
                VoxelIsoValueGrid voxelIsoValueGrid = voxelGrids[n];
                if (voxelIsoValueGrid.id == 4428)
                {
                    //System.out.println("Found it!!");
                }
                boolean disconnected = voxelIsoValueGrid.isDisconnected(row);
                if (disconnected)
                {
                    if (voxelIsoValueGrid.isTooSmall()) // delete: disconnected and too small
                        voxelGrids[n] = null;
                    else // put in storage: disconnected and size ok
                    {
                        voxelIsoValueGrid.initialize();
                        nCreated++;
                        nTotal++;
                    }
                }
                else // !disconnected: add to current patchGrids as we will need it again
                {
                    isoValueGrids.addGrid(voxelIsoValueGrid);
                    nGridsInSet++;
                    nTotal++;
                }
            }
            //			System.out.println("    nonNulls in volumeSet " + nGridsInSet + " wrote out " + nCreated + ". Total " + nTotal);
            number = 0;
        }

        void processRowGridsXX(int row)
        {
            for (int n = 0; n < number; n++)
            {
                VoxelIsoValueGrid voxelIsoValueGrid = voxelGrids[n];
                if (voxelIsoValueGrid.id == 4428)
                {
                    //System.out.println("Found it!!");
                }
                boolean disconnected = voxelIsoValueGrid.isDisconnected(row);
                if (!voxelIsoValueGrid.isTooSmall())
                {
                    if (!disconnected)
                    {
                        isoValueGrids.addGrid(voxelIsoValueGrid);
                        nGridsInSet++;
                    }
                    else // disconnected && more than one point, classInitialize and write out
                    {
                        voxelIsoValueGrid.initialize();
                        nCreated++;
                    }
                    nTotal++;
                }
                else if (disconnected)
                    voxelGrids[n] = null;
            }
            //			System.out.println("    nonNulls in volumeSet " + nGridsInSet + " wrote out " + nCreated + ". Total " + nTotal);
            number = 0;
        }

        /**
         * This IsoValueGridSet is for the row before row just finished.
         * If a grid in this row is disconnected (doesn't have entry in row just finished),
         * then either delete it if it is a single point, or add it to the volume set.
         * If not disconnected, add it to volume set.
         */
        void processDisconnectedGrids(int row)
        {
            int numberRemoved = 0;
            boolean move = false;
            boolean connectedSequence = false;
            int nn = 0;
            int firstConnected = 0; // first connected grid after one or more disconnected
            int lastConnected = 0; // last connectedGrid in connectedGrid sequence beginning with firstConnected
            int n = 0;
            VoxelIsoValueGrid voxelIsoValueGrid = null;
            try
            {
                for (n = 0; n < number; n++)
                {
                    voxelIsoValueGrid = voxelGrids[n];
                    if (voxelIsoValueGrid.isDisconnected(row))
                    {
                        voxelIsoValueGrid.initialize();
                        numberRemoved++;
                        nCreated++;
                        nGridsInSet--;
                        if (move && connectedSequence)
                        {
                            int nConnected = lastConnected - firstConnected + 1;
                            System.arraycopy(voxelGrids, firstConnected, voxelGrids, nn, nConnected);
                            nn += nConnected;
                            connectedSequence = false;
                        }
                        move = true;
                    }
                    else // found a connected grid
                    {
                        if (!move)
                        {
                            nn++;
                        }
                        else
                        {
                            if (connectedSequence)
                            {
                                lastConnected = n;
                            }
                            else
                            {
                                firstConnected = n;
                                lastConnected = n;
                                connectedSequence = true;
                            }
                        }
                    }
                }
                if (move && connectedSequence)
                {
                    int nConnected = lastConnected - firstConnected + 1;
                    if (nConnected < number)
                    {
                        System.arraycopy(voxelGrids, firstConnected, voxelGrids, nn, nConnected);
                        nn += nConnected;
                        connectedSequence = false;
                    }
                }
                number = nn;

                int checkTotal = nCreated + nGridsInSet;
                System.out.println("    removed " + numberRemoved + " nWritten " + nCreated + " nonNulls in set " +
                    nGridsInSet + " total " + nTotal + " check Total " + checkTotal + " nn " + nn);
            }
            catch (Exception e)
            {
                StsException.outputException("StsSeismicVolume.processDisConnectedGrids() failed at n: " + n, e,
                    StsException.WARNING);
            }
        }

        /** Called for the last row only; unconditionally add these points to volume unless they are single points. */
        void addToGrids(IsoValueGridSet volumeGridSet)
        {
            int numberAdded = 0;
            for (int n = 0; n < number; n++)
            {
                VoxelIsoValueGrid voxelIsoValueGrid = voxelGrids[n];
                if (!voxelIsoValueGrid.isTooSmall())
                {
                    volumeGridSet.addGrid(voxelIsoValueGrid);
                    nGridsInSet++;
                    nTotal++;
                    numberAdded++;
                }
                voxelGrids[n] = null;
            }
            number = 0;
            System.out.println("    added " + numberAdded + " in last row. Total " + nTotal);
        }

        /*
        void addGrids(StsBinaryFile isoValueGridsFile)
        {
         for(int n = 0; n < number; n++)
         {
          PatchGrid isoValueGrid = voxelGrids[n];
          if(!isoValueGrid.isSmaller(4))
          {
        isoValueGrid.classInitialize(isoValueGridsFile);
        nWritten++;
        nTotal++;
          }
          voxelGrids[n] = null;
         }
        System.out.println("    nonNulls in volumeSet" + nonNullsInSet + " wrote out " + nWritten + ". Total " + nTotal);
         number = 0;
        }
        */

        void initialize()
        {
            for (int n = 0; n < number; n++)
                if (voxelGrids[n] != null) voxelGrids[n].initialize();
            voxelGrids = null;
        }
    }

    class VoxelTracePoints
    {
        int row;
        int col;
        int nTraceBlobs = 0;
        TraceBlob[] traceBlobs = new TraceBlob[10];

        VoxelTracePoints(int row, int col, int[] values, int n0, int nValues)
        {
            this.row = row;
            this.col = col;
            int nn = n0;
            int index0;
            int index1 = -1;
            int value0, value1 = -1;
            float traceIndexF = -1.0f;
            boolean isMax;
            boolean isInsideBlob = false;
            TraceBlob currentBlob = null;
            int currentBlobVoxelKeyIndex = -1;
            try
            {
                for (int n = 0; n < nValues; n++)
                {

                    value0 = value1;
                    value1 = values[nn++];
                    index0 = index1;
                    index1 = voxelKeyIndexes[value1];
                    if (index1 >= 0 && currentBlob == null) // hit top of blob
                    {
                        currentBlobVoxelKeyIndex = index1;
                        currentBlob = addTraceBlob(row, col, currentBlobVoxelKeyIndex);
                        currentBlob.setTop(n, value0, value1);
                    }
					else if (currentBlob != null && index1 != currentBlobVoxelKeyIndex) // we have exited the currentBlob
                    {
                        if (index1 >= 0) // if this point belongs to a different blob, close the current and start a new one
                        {
                            currentBlob.setBot(n, nValues, value0, value1); // close out current blob: add bottom
                            currentBlobVoxelKeyIndex = index1;
                            currentBlob = addTraceBlob(row, col, currentBlobVoxelKeyIndex); // add new blob
                            currentBlob.setTop(n, value0, value1);
                        }
                        else // we have exited bottom of current blob
                        {
                            currentBlob.setBot(n, nValues, value0, value1);
                            currentBlobVoxelKeyIndex = index1; // should be -1
                            currentBlob = null;
                        }
                    }
                }
                traceBlobs = trimBlobs(traceBlobs, nTraceBlobs);
            }
            catch (Exception e)
            {
                StsException.outputException("StsSeismicVolume.VoxelTracePoints.constructor() failed. value = " + value1 +
                    " traceIndexF " + traceIndexF, e, StsException.WARNING);
            }
        }

        TraceBlob addTraceBlob(int row, int col, int voxelKeyIndex)
        {
            int traceBlobsLength = traceBlobs.length;
            if (nTraceBlobs >= traceBlobsLength)
            {
                TraceBlob[] newTraceBlobs = new TraceBlob[traceBlobsLength + 10];
                System.arraycopy(traceBlobs, 0, newTraceBlobs, 0, nTraceBlobs);
                traceBlobs = newTraceBlobs;
            }
            TraceBlob blob = new TraceBlob(row, col, voxelKeyIndex);
            traceBlobs[nTraceBlobs++] = blob;
            return blob;
        }

        TraceBlob[] trimBlobs(TraceBlob[] blobs, int nBlobs)
        {
            int blobsLength = blobs.length;
            if (nBlobs == blobsLength) return blobs;
            TraceBlob[] newBlobs = new TraceBlob[nBlobs];
            System.arraycopy(blobs, 0, newBlobs, 0, nBlobs);
            return newBlobs;
        }
    }

    class TraceBlob
    {
        int row;
        int col;
        int id = -1;
        int topIndex;
        int botIndex;
        float topIndexF;
        float botIndexF;
        int voxelKeyIndex;
        boolean isTopMax;
        boolean isBotMax;

        TraceBlob(int row, int col, int voxelKeyIndex)
        {
            this.row = row;
            this.col = col;
            this.voxelKeyIndex = voxelKeyIndex;
        }

        /**
         * For this trace which passes thru a blob, set the information defining
         * the top intersection.  traceIndex is the top of the blob.  The min or
         * max value lies between value0 and value1.  Get this interpolated value
         * and set topIndexF as the interpolated index value of the top of blob.
         *
         * @param traceIndex sequential index of this point in the trace
         * @param value0     attribute value of previous trace point
         * @param value1     attribute value of this point
         */
        void setTop(int traceIndex, int value0, int value1)
        {
            boolean null_0 = (value0 == nullUnsignedInt);
            boolean null_1 = (value1 == nullUnsignedInt);
            if (null_0) value0 = 0;
            if (null_1) value1 = 0;

            this.topIndex = traceIndex;

            // if we at the first index, we are starting inside the blob, so set the top to 0
            if (traceIndex == 0)
            {
                int mid = (voxelKeys[voxelKeyIndex].max + voxelKeys[voxelKeyIndex].min) / 2;
                isTopMax = value1 >= mid;
                topIndexF = 0.0f;
            }
            else
            {
                int max = voxelKeys[voxelKeyIndex].max;
                if (value0 > max)
                {
                    // previous value (value0) is greater than max so we are crossing the max value as we enter the blob (value1 must be <= max)
                    if (value1 == value0)
                        topIndexF = traceIndex - 0.5f;
                    else
                        topIndexF = traceIndex - (float) (value1 - max) / (value1 - value0);
                    isTopMax = true;
                }
                else
                {
                    // since we know we have entered blob and value0 must be < max, we also know value0 < min
                    int min = voxelKeys[voxelKeyIndex].min;
                    if (value1 == value0)
                        topIndexF = traceIndex - 0.5f;
                    else
                    {
                        topIndexF = traceIndex - (float) (value1 - min) / (value1 - value0);
                    }
                    isTopMax = false;
                }
            }
        }

        void setBot(int traceIndex, int nTraceValues, int value0, int value1)
        {
            boolean null_0 = (value0 == nullUnsignedInt);
            boolean null_1 = (value1 == nullUnsignedInt);
            if (null_0) value0 = 0;
            if (null_1) value1 = 0;

            this.botIndex = traceIndex;
            if (traceIndex == nTraceValues - 1)
            {
                int mid = (voxelKeys[voxelKeyIndex].max + voxelKeys[voxelKeyIndex].min) / 2;
                isBotMax = !null_1 && value1 >= mid;
                botIndexF = (float) traceIndex;
            }
            else
            {
                int max = voxelKeys[voxelKeyIndex].max;
                if (value1 > max)
                {
                    // we have exited the blob.  if value1 > max, then we must have crossed max on exit
                    if (value1 == value0)
                        botIndexF = traceIndex - 0.5f;
                    else
                    {
                        float f = (float) (value1 - max) / (value1 - value0);
                        if (f < 0.0f || f > 1.0f)
                            System.out.println("Bummer");
                        botIndexF = traceIndex - (float) (value1 - max) / (value1 - value0);
                    }
                    isBotMax = true;
                }
                else
                {
                    // we have exited the blob.  value1 < min, so we crossed min value on exit
                    int min = voxelKeys[voxelKeyIndex].min;
                    if (value1 == value0)
                        botIndexF = traceIndex - 0.5f;
                    else
                    {
                        float f = (float) (value1 - min) / (value1 - value0);
                        if (f < 0.0f || f > 1.0f)
                            System.out.println("Bummer");
                        botIndexF = traceIndex - (float) (value1 - min) / (value1 - value0);
                    }
                    isBotMax = false;
                }
            }
        }
    }

    class VoxelTracePoint // implements Comparable
    {
        int row;
        int col;
        int id = -1;
        float traceIndexF;
        int voxelKeyIndex;
        boolean isMax;

        VoxelTracePoint(int row, int col, float traceIndexF, int voxelKeyIndex, boolean isMax)
        {
            this.row = row;
            this.col = col;
            this.traceIndexF = traceIndexF;
            this.voxelKeyIndex = voxelKeyIndex;
            this.isMax = isMax;
        }
        /*
         void setID(int id)
         {
          this.id = id;
         }

         public int compareTo(Object obj)
         {
          return colorIndex - ( (VoxelTracePoint) obj).colorIndex;
         }
        */
    }

    class VoxelIsoValueGrid
    {
        int id;
        int voxelKeyIndex;
        int rowMin = StsParameters.largeInt;
        int rowMax = -StsParameters.largeInt;
        int colMin = StsParameters.largeInt;
        int colMax = -StsParameters.largeInt;
        int nGridRows, nGridCols;
        TraceBlob[] traceBlobs = new TraceBlob[2];
        int nTraceBlobs = 0;
        float[][] topPointsZ = null;
        float[][] botPointsZ = null;

        VoxelIsoValueGrid(int id, int voxelKeyIndex)
        {
            this.id = id;
            this.voxelKeyIndex = voxelKeyIndex;
            if (indices == null)
            {
                nIndicesAllocated = 4;
                indices = new int[nIndicesAllocated];
                count = new int[nIndicesAllocated];
            }
        }

        void addTraceBlob(TraceBlob blob)
        {
            addToList(blob);
            int row = blob.row;
            int col = blob.col;

            rowMin = Math.min(rowMin, row);
            rowMax = Math.max(rowMax, row);
            colMin = Math.min(colMin, col);
            colMax = Math.max(colMax, col);
        }

        void addToList(TraceBlob blob)
        {
            int length = traceBlobs.length;
            if (nTraceBlobs >= length)
            {
                TraceBlob[] newTraceBlobs = new TraceBlob[2 * length];
                System.arraycopy(traceBlobs, 0, newTraceBlobs, 0, length);
                traceBlobs = newTraceBlobs;
            }
            traceBlobs[nTraceBlobs++] = blob;
        }

        /*
         void classInitialize()
         {
          if (topPointsZ != null)return;

          nRows = rowMax - rowMin + 1;
          nCols = colMax - colMin + 1;
          if (nRows <= 1 && nCols <= 1)return;

//	        System.out.println("Initialize isoValueGrid. rowMin " + rowMin + " rowMax " + rowMax + " colMin " + colMin +
//							   " colMax " + colMax + " nTracePoints " + nTracePoints);

          topPointsZ = new float[nRows][nCols];
          botPointsZ = new float[nRows][nCols];
          for (int row = 0; row < nRows; row++)
          {
        for (int col = 0; col < nCols; col++)
        {
         topPointsZ[row][col] = nullValue;
         botPointsZ[row][col] = nullValue;
        }
          }
          for (int n = 0; n < nTraceBlobs; n++)
          {
        int row = patchPoints[n].row;
        int col = patchPoints[n].col;
        float traceIndexF = patchPoints[n].topIndexF;
        topPointsZ[row - rowMin][col - colMin] = zMin + traceIndexF * zInc;
        traceIndexF = patchPoints[n].botIndexF;
        botPointsZ[row - rowMin][col - colMin] = zMin + traceIndexF * zInc;
          }

          float[][] topXyzVector = new float[nRows*nCols][3];
          float[][] topNormalVector = new float[nRows*nCols][];
          float[][] botXyzVector = new float[nRows*nCols][3];
          float[][] botNormalVector = new float[nRows*nCols][];
          int n = 0;
          float y = yMin + rowMin * yInc;
          for (int row = 0; row < nRows; row++, y += yInc)
          {
        float x = xMin + colMin * xInc;
        for (int col = 0; col < nCols; col++, x += xInc, n++)
        {
         if(topPointsZ[row][col] == nullValue) continue;
         topXyzVector[n][0] = x;
         topXyzVector[n][1] = y;
         topXyzVector[n][2] = topPointsZ[row][col];
         topNormalVector[n] = computeNormal(topPointsZ, row, col, true);
         botXyzVector[n][0] = x;
         botXyzVector[n][1] = y;
         botXyzVector[n][2] = botPointsZ[row][col];
         botNormalVector[n] = computeNormal(botPointsZ, row, col, false);
        }
          }

          int index = 0;
          // count number of sequences
          int nSeq = 0;
          boolean isSequence = false;
          for (int row = 0; row < nRows-1; row++)
          {
        for (int col = 0; col < nCols; col++)
        {
         if (topPointsZ[row][col] != nullValue && topPointsZ[row + 1][col] != nullValue)
         {
          if (!isSequence)
          {
           nSeq++;
           isSequence = true;
          }
         }
         else
          isSequence = false;
        }
        if (isSequence)
        {
         isSequence = false;
        }
          }

          int[][] indexArray = new int[nSeq][];
          int[] indices = new int[2 * (nCols + 1)];
          nSeq = 0;
          int nIndices = 0;
          isSequence = false;
          for (int row = 0; row < nRows-1; row++)
          {
        for (int col = 0; col < nCols; col++)
        {
         if (topPointsZ[row][col] != nullValue && topPointsZ[row + 1][col] != nullValue)
         {
          if (!isSequence)
          {
           nIndices = 0;
           isSequence = true;
          }
          indices[nIndices++] = row * nCols + col;
          indices[nIndices++] = (row + 1) * nCols + col;
         }
         else
         {
          if (isSequence)
          {
           indexArray[nSeq] = new int[nIndices];
           System.arraycopy(indices, 0, indexArray[nSeq], 0, nIndices);
           isSequence = false;
           nSeq++;
          }
         }
        }
        if (isSequence)
        {
         indexArray[nSeq] = new int[nIndices];
         System.arraycopy(indices, 0, indexArray[nSeq], 0, nIndices);
         isSequence = false;
         nSeq++;
        }
          }


          boolean isTopMax = patchPoints[0].isTopMax;
          boolean isBotMax = patchPoints[0].isBotMax;
          if (patchPoints[0].isTopMax)
        voxelVertexArraySets[voxelKeyIndex * 2 + 1].3topXyzVector, topNormalVector, indexArray);
          else
        voxelVertexArraySets[voxelKeyIndex * 2].addVertices(topXyzVector, topNormalVector, indexArray);

          if (patchPoints[0].isBotMax)
        voxelVertexArraySets[voxelKeyIndex * 2 + 1].addVertices(botXyzVector, botNormalVector, indexArray);
          else
        voxelVertexArraySets[voxelKeyIndex * 2].addVertices(botXyzVector, botNormalVector, indexArray);

          // done with tracePoints and gridZ, so free them
          patchPoints = null;
          topPointsZ = null;
          botPointsZ = null;
          topNormalVector = null;
          botNormalVector = null;
         }
        */

        /**
         * This IsoValueGrid contains an rectangular array of traceBlobs, nRows by
         * nCols. We will take the top and bottom points of each blob and form a
         * surface which is one-half cell wider on each side, i.e., we have square
         * cells centered around each original grid point.
         */

        void initialize()
        {
            //			if (topPointsZ != null)return;

            nGridRows = rowMax - rowMin + 1;
            nGridCols = colMax - colMin + 1;
            if (nGridRows < 1 || nGridCols < 1) return;

            //	        System.out.println("Initialize isoValueGrid. rowMin " + rowMin + " rowMax " + rowMax + " colMin " + colMin +
            //							   " colMax " + colMax + " nTracePoints " + nTracePoints);

            float[][] topGridZ = new float[nGridRows][nGridCols];
            float[][] botGridZ = new float[nGridRows][nGridCols];
            for (int row = 0; row < nGridRows; row++)
            {
                for (int col = 0; col < nGridCols; col++)
                {
                    topGridZ[row][col] = nullValue;
                    botGridZ[row][col] = nullValue;
                }
            }
            for (int n = 0; n < nTraceBlobs; n++)
            {
                int row = traceBlobs[n].row;
                int col = traceBlobs[n].col;
                float topIndexF = traceBlobs[n].topIndexF;
                topGridZ[row - rowMin][col - colMin] = zMin + topIndexF * zInc;
                float botIndexF = traceBlobs[n].botIndexF;
                if (botIndexF < topIndexF) botIndexF = topIndexF;
                botGridZ[row - rowMin][col - colMin] = zMin + botIndexF * zInc;
            }

            float[][] topPointsZ = this.computeCellCenteredPoints(topGridZ, nGridRows, nGridCols);
            float[][] botPointsZ = this.computeCellCenteredPoints(botGridZ, nGridRows, nGridCols);

            boolean isTopMax = traceBlobs[0].isTopMax;
            boolean isBotMax = traceBlobs[0].isBotMax;

            // done with these
            if (!voxelDebug)
            {
                topGridZ = null;
                botGridZ = null;
                traceBlobs = null;
            }

            float[][] topXyzVector = new float[(nGridRows + 1) * (nGridCols + 1)][3];
            float[][] botXyzVector = new float[(nGridRows + 1) * (nGridCols + 1)][3];
            float[][] topNormalVector = new float[(nGridRows + 1) * (nGridCols + 1)][];
            float[][] botNormalVector = new float[(nGridRows + 1) * (nGridCols + 1)][];
            //			float minZ = largeFloat;
            //			float maxZ = -largeFloat;

            // construct row y vector
            float[] rowY = new float[nGridRows + 1];
            if (rowMin == 0)
                rowY[0] = yMin;
            else
                rowY[0] = yMin + (rowMin - 0.5f) * yInc;
            for (int row = 1; row < nGridRows; row++)
                rowY[row] = rowY[row - 1] + yInc;
            if (rowMax == nRows - 1)
                rowY[nGridRows] = yMax;
            else
                rowY[nGridRows] = rowY[nGridRows - 1] + yInc;

            // construct col x vector
            float[] colX = new float[nGridCols + 1];
            if (colMin == 0)
                colX[0] = xMin;
            else
                colX[0] = xMin + (colMin - 0.5f) * xInc;

            for (int col = 1; col < nGridCols; col++)
                colX[col] = colX[col - 1] + yInc;
            if (colMax == nCols - 1)
                colX[nGridCols] = xMax;
            else
                colX[nGridCols] = colX[nGridCols - 1] + xInc;

            int n = 0;
            for (int row = 0; row < nGridRows + 1; row++)
            {
                for (int col = 0; col < nGridCols + 1; col++, n++)
                {
                    // if this top point is null, bottom is also; supply null vector for both
                    // as they will be copied to array passed to glNormalPointer
                    if (topPointsZ[row][col] == nullValue)
                    {
                        topNormalVector[n] = nullNormalVector;
                        botNormalVector[n] = nullNormalVector;
                        continue;
                    }
                    topXyzVector[n][0] = colX[col];
                    topXyzVector[n][1] = rowY[row];
                    topXyzVector[n][2] = topPointsZ[row][col];
                    topNormalVector[n] = computeNormal(topPointsZ, row, col, true);
                    botXyzVector[n][0] = colX[col];
                    botXyzVector[n][1] = rowY[row];
                    botXyzVector[n][2] = botPointsZ[row][col];
                    botNormalVector[n] = computeNormal(botPointsZ, row, col, true);
                }
            }
            int nPoints = 2 * nGridRows * (nGridCols + 1);
            if (nPoints > nIndicesAllocated)
            {
                nIndicesAllocated = nPoints * 2;
                indices = new int[nIndicesAllocated];
                count = new int[nIndicesAllocated];
            }
            int nIndices = 0;
            int nTotalIndices = 0;
            int nCount = 0;
            for (int row = 0; row < nGridRows; row++)
            {
                for (int col = 0; col < nGridCols + 1; col++)
                {
                    if (topPointsZ[row][col] != nullValue && topPointsZ[row + 1][col] != nullValue)
                    {
                        indices[nTotalIndices++] = row * (nGridCols + 1) + col;
                        indices[nTotalIndices++] = (row + 1) * (nGridCols + 1) + col;
                        nIndices += 2;
                    }
                    else if (nIndices > 0)
                    {
                        count[nCount++] = nIndices;
                        nIndices = 0;
                    }
                }
                if (nIndices > 0)
                {
                    count[nCount++] = nIndices;
                    nIndices = 0;
                }
            }

            if (isTopMax)
                voxelVertexArraySets[voxelKeyIndex * 2 + 1].addVertices(topXyzVector, topNormalVector, indices, nTotalIndices, count, nCount);
            else
                voxelVertexArraySets[voxelKeyIndex * 2].addVertices(topXyzVector, topNormalVector, indices, nTotalIndices, count, nCount);

            if (isBotMax)
                voxelVertexArraySets[voxelKeyIndex * 2 + 1].addVertices(botXyzVector, botNormalVector, indices, nTotalIndices, count, nCount);
            else
                voxelVertexArraySets[voxelKeyIndex * 2].addVertices(botXyzVector, botNormalVector, indices, nTotalIndices, count, nCount);

            if (voxelDebug)
            {
                topGridZ = null;
                botGridZ = null;
                traceBlobs = null;
            }

            topPointsZ = null;
            botPointsZ = null;

            topXyzVector = null;
            topNormalVector = null;
            botXyzVector = null;
            botNormalVector = null;
        }

        private float[][] computeCellCenteredPoints(float[][] gridZ, int nRows, int nCols)
        {
            float[][] pointsZ = new float[nRows + 1][nCols + 1];

            float p00, p01, p10, p11;
            // compute first cell centered row
            p01 = gridZ[0][0];
            pointsZ[0][0] = p01;
            // computer intermediate points in first row
            for (int col = 1; col < nCols; col++)
            {
                p00 = p01;
                p01 = gridZ[0][col];
                if (p00 != nullValue)
                {
                    if (p01 != nullValue)
                        pointsZ[0][col] = (p00 + p01) / 2;
                    else
                        pointsZ[0][col] = p00;
                }
                else // p00 == nullValue
                {
                    if (p01 != nullValue)
                        pointsZ[0][col] = p01;
                    else
                        pointsZ[0][col] = nullValue;
                }
            }
            // compute last col point in first row
            pointsZ[0][nCols] = p01;
            // compute intermediate rows
            for (int row = 1; row < nRows; row++)
            {
                // compute first col in row
                p01 = gridZ[row - 1][0];
                p11 = gridZ[row][0];
                if (p01 != nullValue)
                {
                    if (p11 != nullValue)
                        pointsZ[row][0] = (p01 + p11) / 2;
                    else
                        pointsZ[row][0] = p01;
                }
                else // p01 == nullValue
                {
                    if (p11 != nullValue)
                        pointsZ[row][0] = p11;
                    else
                        pointsZ[row][0] = nullValue;
                }
                // compute intermediate cols
                for (int col = 1; col < nCols; col++)
                {
                    p00 = p01;
                    p10 = p11;
                    p01 = gridZ[row - 1][col];
                    p11 = gridZ[row][col];
                    int nValues = 0;
                    float sum = 0.0f;
                    if (p00 != nullValue)
                    {
                        sum += p00;
                        nValues++;
                    }
                    if (p01 != nullValue)
                    {
                        sum += p01;
                        nValues++;
                    }
                    if (p10 != nullValue)
                    {
                        sum += p10;
                        nValues++;
                    }
                    if (p11 != nullValue)
                    {
                        sum += p11;
                        nValues++;
                    }
                    if (nValues > 0)
                        pointsZ[row][col] = sum / nValues;
                    else
                        pointsZ[row][col] = nullValue;
                }
                // compute last col
                p00 = p01;
                p10 = p11;
                if (p00 != nullValue)
                {
                    if (p10 != nullValue)
                        pointsZ[row][nCols] = (p00 + p10) / 2;
                    else
                        pointsZ[row][nCols] = p00;
                }
                else // p00 == nullValue
                {
                    if (p10 != nullValue)
                        pointsZ[row][nCols] = p10;
                    else
                        pointsZ[row][nCols] = nullValue;
                }
            }
            // compute first col in last cell centered row
            p01 = gridZ[nRows - 1][0];
            pointsZ[nRows][0] = p01;
            // compute intermediate col points in last row
            for (int col = 1; col < nCols; col++)
            {
                p00 = p01;
                p01 = gridZ[nRows - 1][col];
                if (p00 != nullValue)
                {
                    if (p01 != nullValue)
                        pointsZ[nRows][col] = (p00 + p01) / 2;
                    else
                        pointsZ[nRows][col] = p00;
                }
                else // p00 == nullValue
                {
                    if (p01 != nullValue)
                        pointsZ[nRows][col] = p01;
                    else
                        pointsZ[nRows][col] = nullValue;
                }
            }
            // compute last col point in last row
            pointsZ[nRows][nCols] = p01;

            return pointsZ;
        }

        boolean isDisconnected(int row)
        {
            return rowMax < row;
        }

        boolean isOnePoint()
        {
            return rowMax - rowMin <= 0 && colMax - colMin <= 0;
        }

        boolean isTooSmall()
        {
            if (smallestVoxel == 0) return false;
            return rowMax - rowMin < smallestVoxel && colMax - colMin < smallestVoxel;
        }

        /*
         void write(StsBinaryFile isoValueGridsFile)
         {
          int[] ints = new int[]
        {
        id, rowMin, rowMax, colMin, colMax, nRows, nCols};
          isoValueGridsFile.setIntegerValues(ints);
          isoValueGridsFile.setBooleanValues(new boolean[]
            {isTop});
          for (int row = 0; row <= nRows; row++)
        isoValueGridsFile.setFloatValues(pointsZ[row]);

          nIsoValueGrids++;
         }

         boolean read(StsBinaryFile isoValueGridsFile)
         {
          int[] ints = isoValueGridsFile.getIntegerValues();
          if (ints == null)return false;
          id = ints[0];
          rowMin = ints[1];
          rowMax = ints[2];
          colMin = ints[3];
          colMax = ints[4];
          nRows = ints[5];
          nCols = ints[6];

          boolean[] bools = isoValueGridsFile.getBooleanValues();
          if (bools == null)return false;
          isTop = bools[0];

          pointsZ = new float[nRows + 1][];
          for (int row = 0; row <= nRows; row++)
          {
        pointsZ[row] = isoValueGridsFile.getFloatValues();
        if (pointsZ[row] == null)return false;
          }
          return true;
         }
        */
        /*
         void draw(GL gl)
         {
          if(nRows <= 1 && nCols <= 1) return;
          gl.glBegin(GL.GL_POINTS);
          float y = yMin + rowMin * yInc;
          for (int row = 0; row < nRows; row++, y += yInc)
          {
        float x = xMin + colMin * xInc;
        for (int col = 0; col <= nCols; col++, x += xInc)
        {
         float z = pointsZ[row][col];
         if (z == nullValue) continue;
         gl.glVertex3f(x, y, z);
        }
          }
          gl.glEnd();
         }
        */
        /*
        void draw(GL gl)
        {
         float zRowM1, zRowP1, zColM1, zColP1;
         if (nRows <= 1 && nCols <= 1)return;
         boolean drawing = false;
         float zRowp1 = nullValue;
         float zRowm1 = nullValue;
         float y1 = yMin + rowMin * yInc;
         for (int row = 0; row < nRows; row++)
         {
          float y0 = y1;
          y1 += yInc;
          float x = xMin + colMin * xInc;
          for (int col = 0; col <= nCols; col++, x += xInc)
          {
        float z0 = pointsZ[row][col];
        float z1 = pointsZ[row + 1][col];
        if (z0 != nullValue && z1 != nullValue)
        {
         if (!drawing)
         {
          gl.glBegin(GL.GL_QUAD_STRIP);
          drawing = true;
         }

         float[] normal = computeNormal(row, col);
         gl.glNormal3fv(normal);
         gl.glVertex3f(x, y0, z0);
         normal = computeNormal(row + 1, col);
         gl.glNormal3fv(normal);
         gl.glVertex3f(x, y1, z1);
        }
        else
        {
         if (drawing)
         {
          gl.glEnd();
          drawing = false;
         }
        }
          }
          if (drawing)
          {
        gl.glEnd();
        drawing = false;
          }
         }
        }
        */
        /*
         void debugPrint()
         {
          System.out.println(" PatchGrid mainDebug: id " + id + "rowMin " + rowMin + " rowMax " + rowMax + " colMin " +
           colMin + " colMax " + colMax);
          System.out.println("                     pointsZ[rowMin][colMin] " + pointsZ[0][0] +
           "pointsZ[rowMax][colMax] " + pointsZ[nRows][nCols]);
         }
        */
        float[] computeNormal(float[][] pointsZ, int row, int col, boolean isTop)
        {
            float dZdY, dZdX;
            float zCenter = pointsZ[row][col];
            if (row < nGridRows - 1 && pointsZ[row + 1][col] != nullValue)
                dZdY = (pointsZ[row + 1][col] - zCenter) / yInc;
            else if (row > 0 && pointsZ[row - 1][col] != nullValue)
                dZdY = (zCenter - pointsZ[row - 1][col]) / yInc;
            else
                dZdY = 0.0f;

            if (col < nGridCols - 1 && pointsZ[row][col + 1] != nullValue)
                dZdX = (pointsZ[row][col + 1] - zCenter) / xInc;
            else if (col > 0 && pointsZ[row][col - 1] != nullValue)
                dZdX = (zCenter - pointsZ[row][col - 1]) / xInc;
            else
                dZdX = 0.0f;
            if (isTop)
                return new float[]
                    {
                        dZdX, dZdY, -1.0f};
            else
                return new float[]
                    {
                        -dZdX, -dZdY, 1.0f};
        }
    }

    class VoxelVertexArrays
    {
        VoxelVertexArray[] voxelVertexArrays = new VoxelVertexArray[10];
        VoxelVertexArray currenVoxelArray = null;
        int nVertexArrays = 0;
        long memoryUsed = 0;

        public VoxelVertexArrays()
        {
            addArray();
        }

        void addArray()
        {
            if (currenVoxelArray != null)
            {
                memoryUsed += 8 * currenVoxelArray.nCoordinates + 4 * currenVoxelArray.nTotalIndices;
            }
            int length = voxelVertexArrays.length;
            if (nVertexArrays == length)
            {
                VoxelVertexArray[] newVoxelVertexArrays = new VoxelVertexArray[length + 10];
                System.arraycopy(voxelVertexArrays, 0, newVoxelVertexArrays, 0, length);
                voxelVertexArrays = newVoxelVertexArrays;
            }
            currenVoxelArray = new VoxelVertexArray();
            voxelVertexArrays[nVertexArrays++] = currenVoxelArray;
        }

        void addVertices(float[][] vertexVector, float[][] normalVector, int[] indices, int nIndices, int[] count, int nCount)
        {
            if (currenVoxelArray.addVertices(vertexVector, normalVector, indices, nIndices, count, nCount)) return;
            //			currenVoxelArray.trim();
            addArray();
            if (!currenVoxelArray.addVertices(vertexVector, normalVector, indices, nIndices, count, nCount))
            {
                StsException.systemError("StsSeismicVolume.VoxelVertexArrays.addVertices() failed.");
            }
        }

        int getNVertices()
        {
            int nVertices = 0;
            for (int n = 0; n < voxelVertexArrays.length; n++)
            {
                if (voxelVertexArrays[n] == null) continue;
                nVertices += voxelVertexArrays[n].nTotalVertices;
            }
            return nVertices;
        }

        long getMemoryUsed()
        {
            if (currenVoxelArray == null) return memoryUsed;
            return memoryUsed + 8 * currenVoxelArray.nCoordinates + 4 * currenVoxelArray.nTotalIndices;
        }

        void trim()
        {
            //			currenVoxelArray.trim();
            VoxelVertexArray[] newVoxelVertexArrays;
            if (nVertexArrays < voxelVertexArrays.length)
            {
                newVoxelVertexArrays = new VoxelVertexArray[nVertexArrays];
                System.arraycopy(voxelVertexArrays, 0, newVoxelVertexArrays, 0, nVertexArrays);
                voxelVertexArrays = newVoxelVertexArrays;
            }
        }


        boolean hasVBO(GL gl)
        {
            if (vboChecked) return hasVBO;

            if (gl.isExtensionAvailable("GL_ARB_vertex_buffer_object"))
                hasVBO = true;

            vboChecked = true;

            return hasVBO;
        }

        boolean hasPointParams(GL gl)
        {
            if (ppChecked) return hasPP;

            if (gl.isExtensionAvailable("GL_ARB_point_parameters"))
                hasPP = true;

            if (gl.isExtensionAvailable("GL_EXT_point_parameters"))
                hasPP = true;

            ppChecked = true;

            return hasPP;
        }

        double sqr(double a) { return (a * a); }
	public float computeFloatFromByte(byte byteValue)
	{
		return StsMath.signedByteToFloat(byteValue, dataMin, dataMax);
	}

	public void display(StsGLPanel glPanel)
	{
		GL gl = glPanel.getGL();

		if (!isVisible)return;

		if (!currentModel.getProject().canDisplayZDomain(getZDomain()))return;

		byte projectZDomain = currentModel.getProject().getZDomain();
		if (projectZDomain != zDomainDisplayed)
			zDomainDisplayed = projectZDomain;

        if (colorscale != null  && colorscale.canDrawVoxels()) drawIsoSurfaces(gl);
	}

	public void deleteDisplayList(GL gl)
	{
		if (displayListNum > 0)
		{
			gl.glDeleteLists(displayListNum, 1);
			displayListNum = 0;
		}
	}

	void drawIsoSurfaces(GL gl)
	{
		if (!checkSetupByteOrFloatRowBlocks())return;
		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		if (fileMapRowFloatBlocks != null)
			drawFloatIsoSurfaces(gl);
		else
			drawByteIsoSurfaces(gl);
	}

	public boolean checkSetupByteOrFloatRowBlocks()
	{
		isDataFloat = rowFloatFilename != null;
		if (isDataFloat)
		{
			if (fileMapRowFloatBlocks != null)return true;
			return setupReadRowFloatBlocks();
		}
		return true;  // using byteBlocks
	/*
		else
		{
			if (fileMapRowByteBlocks != null)return true;
			return setupRowByteBlocks();
		}
	*/
	}

	private void drawFloatIsoSurfaces(GL gl)
	{
	}

	public Object getDisplayAttribute() { return currentAttribute; }

	public void setDisplayAttribute(Object attribute)
	{
		if(currentAttribute == attribute) return;
		currentAttribute = (Attribute)attribute;
		basemapDisplayChanged = true;
		currentModel.getCursor3d().clearTextureDisplays(ZDIR);
		currentModel.viewObjectChanged(this, this);
	}

    public boolean displayingBasemap()
    {
        return currentAttribute != nullAttribute;
    }

    public boolean hasBasemapChanged()
	{
		return basemapDisplayChanged;
	}

	public boolean displayingAttribute()
	{
		return currentAttribute != nullAttribute;
	}

	public boolean basemapTextureChanged()
	{
		if (basemapDisplayChanged)
		{
			basemapDisplayChanged = false;
			return true;
		}
		else
			return false;
	}

	public void initDisplayAttributes()
	{
		currentAttribute = nullAttribute;
		displayAttributes = (Attribute[])StsMath.arrayAddElement(displayAttributes, nullAttribute);
		if(attributeNames == null) return;
		for (int n = 0; n < attributeNames.length; n++)
			displayAttributes = (Attribute[])StsMath.arrayAddElement(displayAttributes, new Attribute(attributeNames[n]));
	}

    private void setPointParams(GL gl)
    {
        if (!hasPointParams(gl)) return;
        int error = gl.glGetError();
        float maxSize = 50.f;

        gl.glPointParameterf(GL.GL_POINT_SIZE_MIN, 2.0f);
        double[] projectionMatrix = new double[16];
        int[] viewport = new int[4];
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, projectionMatrix, 0);
        gl.glGetIntegerv(GL.GL_VIEWPORT, viewport, 0);
        double H = viewport[2];
        double h = 2.0 / projectionMatrix[0];
        double D0 = Math.sqrt(2.0 * H / h);
        double k = 1.0 / (1.0 + 2 * sqr(1. / projectionMatrix[0]));
        float[] atten = new float[3];
        atten[0] = 1.f;
        atten[1] = 0.0f;
        k /= 500.f;
        atten[2] = (float) sqr(1 / D0) * (float) k;
        //System.out.println(" atten "+atten[2]);
        //atten[2] = 0.000001f;

        gl.glPointParameterfv(GL.GL_DISTANCE_ATTENUATION_EXT, atten, 0);
        error = gl.glGetError();
        if (error != 0)
        {
            GLU glu = new GLU();
            System.out.println("pointParams err code " + error + " " + glu.gluErrorString(error));
        }
        System.out.println("atten " + atten[2]);
        float[] ft = new float[1];
        gl.glGetFloatv(GL.GL_POINT_SIZE_MAX, ft, 0);
        //System.out.println("point max "+ft[0]);
        gl.glPointSize(ft[0] > maxSize ? maxSize : ft[0]);
        gl.glPointParameterf(GL.GL_POINT_SIZE_MAX, ft[0] > maxSize ? maxSize : ft[0]);
    }

    private int genObj(GL gl)
    {
        int[] tmp = new int[1];
        gl.glGenBuffers(1, tmp, 0);
        return tmp[0];
    }

    void draw(GL gl)
    {
        int a = 0, n = 0;
        try
        {
            gl.glPointSize(5.0f);
            setPointParams(gl);

            //gl.glEnable(GL.GL_POINT_SMOOTH);

            gl.glEnableClientState(GL.GL_VERTEX_ARRAY);
            gl.glEnableClientState(GL.GL_NORMAL_ARRAY);

            if (hasVBO(gl))
            {
                for (a = 0; a < nVertexArrays; a++)
                {
                    if (voxelVertexArrays[a].vertexBufferObject == 0)
                    {
                        voxelVertexArrays[a].vertexBufferObject = genObj(gl);
                        voxelVertexArrays[a].vertexBuffer.rewind();
                        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, voxelVertexArrays[a].vertexBufferObject);
                        gl.glBufferData(GL.GL_ARRAY_BUFFER_ARB, voxelVertexArrays[a].vertexBuffer.capacity() * 4,
                            voxelVertexArrays[a].vertexBuffer, GL.GL_STATIC_DRAW);
                    }
                    if (voxelVertexArrays[a].normalBufferObject == 0)
                    {
                        voxelVertexArrays[a].normalBufferObject = genObj(gl);
                        voxelVertexArrays[a].normalBuffer.rewind();
                        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, voxelVertexArrays[a].normalBufferObject);
                        gl.glBufferData(GL.GL_ARRAY_BUFFER_ARB, voxelVertexArrays[a].normalBuffer.capacity() * 4,
                            voxelVertexArrays[a].normalBuffer, GL.GL_STATIC_DRAW);
                    }
                }
            }

            for (a = 0; a < nVertexArrays; a++)
            {
                FloatBuffer vertexBuffer = voxelVertexArrays[a].vertexBuffer;
                FloatBuffer normalBuffer = voxelVertexArrays[a].normalBuffer;
                vertexBuffer.rewind();
                normalBuffer.rewind();

                if (hasVBO(gl))
                {
                    gl.glEnableClientState(GL.GL_VERTEX_ARRAY);
                    gl.glBindBuffer(GL.GL_ARRAY_BUFFER, voxelVertexArrays[a].vertexBufferObject);
                    gl.glVertexPointer(3, GL.GL_FLOAT, 0, 0);


                    {
                        gl.glEnableClientState(GL.GL_NORMAL_ARRAY);
                        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, voxelVertexArrays[a].normalBufferObject);
                        gl.glNormalPointer(GL.GL_FLOAT, 0, 0);
                    }

                }
                else
                {
                    gl.glVertexPointer(3, GL.GL_FLOAT, 0, vertexBuffer);
                    gl.glNormalPointer(GL.GL_FLOAT, 0, normalBuffer);
                }

                int nTotalCount = voxelVertexArrays[a].nTotalCount;
                ShortBuffer indexBuffer = voxelVertexArrays[a].indexBuffer;
                indexBuffer.rewind();
                int[] countVector = voxelVertexArrays[a].countVector;
                //					IntBuffer countBuffer = voxelVertexArrays[a].countBuffer;
                //					countBuffer.rewind();
                int offset = 0;
                //					nTotalCount = 2;

                if (drawQuads)
                {
                    for (n = 0; n < nTotalCount; n++)
                    {
                        int count = countVector[n];
                        //						int count = countBuffer.get();
                        indexBuffer.position(offset);
                        if (debug)
                        {
                            short[] indices = new short[count];
                            indexBuffer.get(indices);

                            System.out.print(" offset " + offset + " count " + count);
                            for (int nn = 0; nn < 2; nn++)
                            {
                                int i = 3 * indices[nn];
                                float x = vertexBuffer.get(i);
                                float y = vertexBuffer.get(i + 1);
                                float z = vertexBuffer.get(i + 2);
                                System.out.print(" i,x,y,z " + i + " " + x + " " + y + " " + z);
                            }
                            System.out.println();
                            indexBuffer.position(offset);
                        }
                        gl.glDrawElements(GL.GL_QUAD_STRIP, count, GL.GL_UNSIGNED_SHORT, indexBuffer);
                        offset += count;
                    }
                }
                else
                {
                    gl.glDrawArrays(GL.GL_POINTS, 0, voxelVertexArrays[a].nTotalVertices);
                }
            }
            gl.glDisableClientState(GL.GL_VERTEX_ARRAY);
            gl.glDisableClientState(GL.GL_NORMAL_ARRAY);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.draw() failed for vertexArray " + a + " indexArray " + n, e, StsException.WARNING);
        }
    }
    /*
     private void validateArrays(VoxelVertexArray[] voxelVertexArrays, int vertexArrayIndex, int indexArrayIndex)
     {
         VoxelVertexArray vertexArray = voxelVertexArrays[vertexArrayIndex];
         if (vertexArray == null)
         {
             StsException.systemError("ValidateArrays failed. VoxelVertexArray " + vertexArrayIndex + " is null.");
             return;
         }
         int[] indexArray = vertexArray.indexArrays[indexArrayIndex];
         if (indexArray == null)
         {
             StsException.systemError("ValidateArrays failed. IndexArray " + indexArrayIndex + " for VoxelVertexArray " + vertexArrayIndex + " is null.");
         }
         FloatBuffer vertexXYZs = vertexArray.vertexXYZs;
//			int nVertexXYZs = vertexXYZs.length;
         FloatBuffer vertexNormals = vertexArray.vertexNormals;
//			int nVertexNormals = vertexNormals.length;
         float minZ = largeFloat;
         float maxZ = -largeFloat;
         for (int n = 0; n < indexArray.length; n++)
         {
             int vertexIndex = indexArray[n];
             if (vertexIndex < 0 || vertexIndex > nVertexXYZs - 1)
             {
                 StsException.systemError("ValidateArrays failed. IndexArray " + n + " for VoxelVertexArray " + vertexArrayIndex +
                                          " has index of " + vertexIndex + " which is out of range 0 to " + nVertexXYZs);
                 return;

             }
             int coordinateIndex = vertexIndex * 3;
             float x = vertexXYZs[coordinateIndex];
             float y = vertexXYZs[coordinateIndex + 1];
             float z = vertexXYZs[coordinateIndex + 2];
             if (x < xMin - xInc || x > xMax + xInc || y < yMin - yInc || y > yMax + yInc || z < 0.0f || z > 2500.0f)
             {
                 StsException.systemError("ValidateArrays failed. Vertex index " + vertexIndex + " for VoxelVertexArray " + vertexArrayIndex +
                                          " has bad coordinates " + x + " " + y + " " + z);
                 return;
             }
             minZ = Math.min(minZ, z);
             maxZ = Math.max(maxZ, z);
         }
         if (Math.abs(minZ - maxZ) > 100.0f)
         {
             System.out.println("minZ " + minZ + " maxZ " + maxZ + " is large difference.");
         }
     }
    */
    }

    class VoxelVertexArray
    {
        FloatBuffer vertexBuffer;
        FloatBuffer normalBuffer;
        ShortBuffer indexBuffer;
        int[] countVector;
        int nTotalVertices = 0;
        int nTotalIndices = 0;
        int nTotalCount = 0;
        //		int nVerticesMax = 0;
        int nCoordinates = 0;
        int vertexMax; // maximum number of vertices
        int indexMax; // maximum number of indices
        int countMax; // maximum number of quadStrips
        int vertexBufferObject;
        int normalBufferObject;

        public VoxelVertexArray()
        {
            try
            {
                if (debug) System.out.println("Allocating Buffer for vertex arrays.  indexBufferMax: " + indexMax);
                indexMax = Short.MAX_VALUE - 1;
                vertexMax = indexMax;
                countMax = indexMax / 4;
                vertexBuffer = BufferUtil.newFloatBuffer(3 * vertexMax);
                normalBuffer = BufferUtil.newFloatBuffer(3 * vertexMax);
                indexBuffer = BufferUtil.newShortBuffer(indexMax);
                countVector = new int[countMax];

                int nCoordinates = 0;
            }
            catch (Exception e)
            {
                StsException.outputException("StsSeismicVolume.VoxelVertexArray.constructor() failed. nVerticesMax: " + indexMax,
                    e, StsException.WARNING);
            }
        }

        /**
         * This set of vertex and normal pairs describes a series of quadStrips.
         * The number of quadStrips is indicesVector.length.  For each quadStrip
         * there is a series of indices of the vector&normals for that strip.
         */

        boolean addVertices(float[][] newXYZs, float[][] newNormals, int[] indices, int nIndices, int[] count, int nCount)
        {
            /*
            if(nIndexArrays == 5768)
            {
             System.out.println("Found 5768.");
            }
            */
            int nNewVertices = newXYZs.length;
            boolean allocationOk = true;
            {
                if (nNewVertices + nTotalVertices > vertexMax) allocationOk = false;
                if (nCount + nTotalCount > countMax) allocationOk = false;
                if (nIndices + nTotalIndices > indexMax) allocationOk = false;
            }
            if (!allocationOk)
            {
                float floatBufferFill = 100.0f * nTotalVertices / vertexMax;
                float indexBufferFill = 100.0f * nTotalIndices / indexMax;
                float countBufferFill = 100.0f * nTotalCount / countMax;
                System.out.println("Buffers full.  floatBufferFill: " + floatBufferFill + " indexBufferFill: " + indexBufferFill + " countBufferFill: " + countBufferFill);
                return false;
            }
            int firstVertexIndex = nTotalVertices;
            for (int n = 0; n < nNewVertices; n++)
            {
                vertexBuffer.put(newXYZs[n]);
                newXYZs[n] = null;
                normalBuffer.put(newNormals[n]);
                newNormals[n] = null;
                nCoordinates += 3;
            }
            nTotalVertices += nNewVertices;

            for (int n = 0; n < nIndices; n++)
            {
                int index = indices[n] + firstVertexIndex;
                indexBuffer.put((short) index);
            }
            nTotalIndices += nIndices;

            for (int n = 0; n < nCount; n++)
                countVector[nTotalCount++] = count[n];
            return true;
        }
    }

    public void setBornDate(String born)
    {
        if (!StsDateFieldBean.validateDateInput(born))
        {
            bornField.setValue(StsDateFieldBean.convertToString(bornDate));
            return;
        }
        super.setBornDate(born);
    }

    public void setDeathDate(String death)
    {
        if (!StsDateFieldBean.validateDateInput(death))
        {
            deathField.setValue(StsDateFieldBean.convertToString(deathDate));
            return;
        }
        super.setDeathDate(death);
    }

    private float getVelocityScalar()
    {
        if (this.oneOrTwoWayVelocity.equals(StsParameters.TWO_WAY_VELOCITY))
          return 2.0f;
        return 1.0f;
    }
    public float getScaledValue(byte byteValue)
    {
        if (byteValue == StsParameters.nullByte) return nullValue;
        float f = (float) StsMath.signedByteToUnsignedInt(byteValue);
        float value = dataMin + (f / 254) * (dataMax - dataMin);

        value *=getVelocityScalar();
        return value;
    }
    public String getUnits()
    {
        if (velocityUnits != null)
             return velocityUnits;
        else
             return "";
    }

    public float getDataMin() { return this.dataMin * getVelocityScalar(); }

    public float getDataMax() { return this.dataMax * getVelocityScalar(); }

    public void setDataRange(float dataMin, float dataMax)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        if(colorscale != null)
            colorscale.setRange(dataMin * getVelocityScalar(), dataMax* getVelocityScalar());
        else
            initializeColorscale();
    }

    public float getMinDepthAtTime(float time)
    {
        StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
	    return velocityModel.getMinDepthAtTime(time);
    }

    public float getMaxDepthAtTime(float time)
    {
        StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
	    return velocityModel.getMaxDepthAtTime(time);
    }

    public String getTypeAsString() { return "Amplitude"; }

    public String getUnitsString() { return "none"; }

    static public StsSeismicVolume[] getZDomainVolumes(StsSeismicVolume[] volumes, byte currentZDomain)
    {
        StsSeismicVolume[] domainVolumes = (StsSeismicVolume[])StsMath.arraycopy(volumes);
        int nVolumes = domainVolumes.length;
        int nDomainVolumes = 0;
        for(int n = 0; n < nVolumes; n++)
        {
            byte volumeZDomain = volumes[n].getZDomain();
            if(volumeZDomain == currentZDomain)
                domainVolumes[nDomainVolumes++] = volumes[n];
        }
        domainVolumes = (StsSeismicVolume[])StsMath.trimArray(domainVolumes, nDomainVolumes);
        return domainVolumes;
    }
}

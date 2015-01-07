package com.Sts.Types;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.nio.*;

public class StsSeismicBoundingBox extends StsRotatedGridBoundingBox implements StsXYGridable, Serializable, Cloneable
{
    /** Histogram of the data distribution */
    public float[] dataHist = new float[255];

    /** directory containing this seismic volume */
    public String segyDirectory;

    /** directory where S2S Volumes are written */
    public String stsDirectory;

    /** filename for this segy volume */
    public String segyFilename;

    /** Original file is archived */
    public boolean isArchived = false;

    /** filename of cube (not complete path) */
    public String stemname;

    /** filename for output inline planes cube */
    public String rowCubeFilename = null;

    /** filename for output crossline planes cube */
    public String colCubeFilename = null;

    /** filename for output trace planes cube */
    public String sliceCubeFilename = null;

    /** filename for output trace planes cube */
    public String attributeFilename = null;

    /** filename for output of native resolution inline planes cube */
    public String rowFloatFilename = null;

    /** SegY last modified date:time */
    public long segyLastModified;

    /** time the s2s files written (header, row, col, slice, float files */
//    public long outputTime;

    /** true if crosslineAngle is 90 degrees CCW from line angle */
    public boolean isXLineCCW = true;

    /** Indicates cube is regular (rectangular) */
    public boolean isRegular = true;

    /** Indicates whether cube is time or depth */
    public String zDomain = StsParameters.TD_NONE_STRING;

    /** String indicating the horizontal units for this volume */
    public String horzUnitsString = StsParameters.D_NONE_STRING;

    /** String indicating the vertical units for this volume */
    public String vertUnitsString = StsParameters.T_NONE_STRING;

    /** Attribute Names */
    public String[] attributeNames;

    /** number of attributes saved per traceHeader */
    public int nAttributes;

    /** volume type assigned by IO processing routine: prestack, poststack, VSP-prestack, etc */
    public byte volumeType = StsParameters.NONE;
    //TODO this could/should be moved to StsSeismicLine perhaps since it is only used for 2d and has subclasses StsSegyLine2d and StsSeismicLine2d.
    public float cdpInterval = 10.0f;

    /** Indicates if values are to be clipped and if so, clipped to null or min/max */
    public byte clipType = CLIP_NONE;
    /** indicates this volume has nulls specified, either an S2Snull for clipped nulls or padded traces,
     *  or a null value defined by the user (see userNullValue).
     */
    public boolean hasNulls = false;
    /** Indicates user has provided a null value.  Samples in the volume with this value will be displayed as nulls,
     *  but unlike S2Snulls, the interpolation option is applied, they will replace by an interpolated value.
     */
    public boolean hasUserNull;
    /** Null value defined by the user.  These values will be displayed as nulls, but unlike S2Snulls, the interpolation
     *  option will replace them with an interpolated value.
     */
    public float userNull = nullValue;

    /** for velocity volumes **/
	public String velocityUnits = StsParameters.VEL_UNITS_NONE;
    public String oneOrTwoWayVelocity = StsParameters.NONE_STRING;

    /** has segy data information; used in SEGY processing by prestack/poststack 2d/3d and vsp processing */
    transient public StsSegyData segyData;
    /** status string used in SEGY I/O processing */
    transient public int status = STATUS_OK;
    transient public String statusString = STATUS_OK_STR;
    public String assocLineName = null;
    /**
     * zDomain currently being displayed. Changing domains requires building new display lists and textures;
     * in which case zDomainDisplayed is set to none, display() method deletes display lists, rebuilds displays
     * for current project zDomain and sets ZDomainDisplayed to this zDomain.
     */
    transient public byte zDomainDisplayed = StsParameters.TD_NONE;
    /** filename of ascii header file for this volume */
    transient public String headerFilename = null;
    /** prefix of files that constitute one volume; file names are ...groupName.stemName... */
//	transient public String groupname = groupNone;

    transient protected float verticalScalar = 1.0f;
    transient protected float horizontalScalar = 1.0f;

    transient protected KnownPoint[] knownPoints;
    transient protected boolean useKnownPoints = false;

    /** Total samples in each of 255 steps */
    transient protected int dataCnt[] = new int[255];
    transient protected int ttlHistogramSamples = 0;

    transient public boolean isDataFloat = false; // if true, floatRowFilename will be defined/persisted

    /**
     * data is scaled to byte range: b = 254*(f - dataMin)/(dataMax - dataMin) which can be written more efficiently as:
     * b = scaleOffset + f*scale where scale = 254/(dataMax - dataMin) and scaleOffset = -dataMin*scale
     */
    transient protected float scale = 0.0f;
    transient protected float scaleOffset = 0.0f;
    //    transient public byte scaleType = SCALE_TRACE;
    transient String dataRangeType;

    transient static final boolean debug = false;

    static final public double doubleNullValue = StsParameters.doubleNullValue;
    static final public float nullValue = StsParameters.nullValue;

    static final public String suffixFilter = "sgy";
    static final public String headerFormat = "txt";
    static final public String group3d = "seis3d";
    static final public String group2d = "seis2d";
    static final public String group3dPrestack = "prestack3d";
    static final public String group2dPrestack = "prestack2d";
    static final public String groupVsp = "seisVsp";
    static final public String groupNone = "none";
    static final public String attributeFormat = "attributes";
    static final public String binFormat = "bin";
    static final public String byteFormat = "bytes";
    static final public String floatFormat = "floats";
    static final public String inline = "inline";
    static final public String xline = "xline";
    static final public String trace = "trace";

    static final public String STATUS_GRID_BAD_STR     = "Lines bad";
    static final public String STATUS_GEOMETRY_BAD_STR = "Coors bad";
    static final public String STATUS_TRACES_BAD_STR   = "Traces bad";
    static final public String STATUS_HEADER_BAD_STR   = "Header bad";
    static final public String STATUS_FILE_BAD_STR     = "File bad";
    static final public String STATUS_UNKNOWN_STR      = "Unknown";
    static final public String STATUS_OK_STR           = "OK";
    static final public String STATUS_FILE_OK_STR      = "File OK";
    static final public String STATUS_HEADER_OK_STR    = "Header OK";
    static final public String STATUS_TRACES_OK_STR    = "Traces OK";
    static final public String STATUS_GEOMETRY_OK_STR  = "Coors OK";
    static final public String STATUS_GRID_OK_STR      = "Lines OK";

    static final public byte STATUS_GEOMETRY_BAD     = 0;
    static final public byte STATUS_GRID_BAD = 1;
    static final public byte STATUS_TRACES_BAD   = 2;
    static final public byte STATUS_FILE_BAD   = 3;
    static final public byte STATUS_HEADER_BAD     = 4;
    static final public byte STATUS_OK      = 5;
    static final public byte STATUS_HEADER_OK      = 6;
    static final public byte STATUS_FILE_OK    = 7;
    static final public byte STATUS_TRACES_OK    = 8;
    static final public byte STATUS_GRID_OK  = 9;
    static final public byte STATUS_GEOMETRY_OK      = 10;

    static final public byte CLIP_NONE = 0;
    static final public byte CLIP_NULL = 1;
    static final public byte CLIP_MINMAX = 2;

    static final public String CLIP_NONE_STRING = "No clip: scale to full range";
    static final public String CLIP_NULL_STRING = "Null clip: null outside clip range";
    static final public String CLIP_MINMAX_STRING = "Min/max clip: min or max outside clip range";

    static final public String[] CLIP_STRINGS = new String[] { CLIP_NONE_STRING , CLIP_NULL_STRING , CLIP_MINMAX_STRING };

    static final public String[] statusText =
    { STATUS_GEOMETRY_BAD_STR, STATUS_GRID_BAD_STR, STATUS_TRACES_BAD_STR, STATUS_FILE_BAD_STR, STATUS_HEADER_BAD_STR,
      STATUS_OK_STR, STATUS_HEADER_OK_STR, STATUS_FILE_OK_STR, STATUS_TRACES_OK_STR, STATUS_GRID_OK_STR, STATUS_GEOMETRY_OK_STR
    };

    static final public byte SCALE_TRACE = 0;
    static final public byte SCALE_NORMALIZED = 1;
    static final public byte SCALE_VOLUME = 2;
    static final public byte SCALE_TRACE_NORMALIZED = 3;

    static public final String clipDataRange = "Clip data range.";
    static public final String fullDataRange = "Full data range.";
    static public final String symmetricDataRange = "Full symmetric data range.";
    static public final String[] dataRangeTypes = new String[]{clipDataRange, fullDataRange, symmetricDataRange};

    final static public Class packSuperClass = StsBoundingBox.class;
    //	transient long segyFileSize;
    transient public StsSegyIO segyIO = null;

    public StsSeismicBoundingBox()
    {
    }

    public StsSeismicBoundingBox(boolean persistent)
    {
        super(persistent);
    }

    public StsSeismicBoundingBox(StsFile file, String stsDirectory)
    {
        super(false);
        segyDirectory = file.getDirectory();
        segyFilename = file.getFilename();
        File sfile = new File(segyDirectory + segyFilename);
        segyLastModified = sfile.lastModified();
        this.stsDirectory = stsDirectory;
    }

    public StsSeismicBoundingBox(StsRotatedGridBoundingBox boundingBox, boolean persistent)
    {
        super(boundingBox, persistent);
    }


    public StsSeismicBoundingBox(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax)
    {
        super(xMin, xMax, yMin, yMax, zMin, zMax);
    }

    public void setZDomain(byte zDomain)
    {
        this.zDomain = StsParameters.TD_ALL_STRINGS[zDomain];
        if (isPersistent()) dbFieldChanged("zDomain", this.zDomain);
    }

    public String getZDomainString()
    {
        return zDomain;
    }

    public void setZDomainString(String zDomainString)
    {
        this.zDomain = zDomainString;
        dbFieldChanged("zDomain", this.zDomain);
    }

    public void setHorzUnitsString(String horzUnitsString)
    {
        this.horzUnitsString = horzUnitsString;
        horizontalScalar = getScalar(horzUnitsString);
    }

    private float getScalar(String unitsString)
    {
        if (zDomain == StsParameters.TD_TIME_STRING)
        {
            byte unitsType = StsParameters.getDistanceUnitsFromString(unitsString);
            return currentModel.getProject().getTimeScalar(unitsType);
        }
        else
        {
            byte unitsType = StsParameters.getDistanceUnitsFromString(unitsString);
            return currentModel.getProject().getDepthScalar(unitsType);
        }
    }

    public String getHorzUnitsString() { return horzUnitsString; }


    public void setVertUnitsString(String vertUnitsString)
    {
        this.vertUnitsString = vertUnitsString;
        verticalScalar = getScalar(vertUnitsString);
    }

    public String getVertUnitsString() { return vertUnitsString; }

    public byte getZDomain() { return StsParameters.getZDomainFromString(zDomain); }

    public void setSegyDirectory(String segyDirectory)
    {
        this.segyDirectory = segyDirectory;
    }

    public String getSegyDirectory()
    {
        return segyDirectory;
    }

    public void setStsDirectory(String stsDirectory)
    {
        this.stsDirectory = stsDirectory;
    }

    public String getStsDirectory()
    {
        return stsDirectory;
    }

    public void setSegyFilename(String segyFilename)
    {
        this.segyFilename = segyFilename;
    }

    public String getSegyFilename()
    {
        return segyFilename;
    }

    public void setStemname(String stemname)
    {
        this.stemname = stemname;
    }

    public String getStemname()
    {
        return stemname;
    }

    public void setRowCubeFilename(String rowCubeFilename)
    {
        this.rowCubeFilename = rowCubeFilename;
    }

    public String getRowCubeFilename()
    {
        return rowCubeFilename;
    }

    public File getRowCubeFile()
    {
        return new File(stsDirectory, rowCubeFilename);
    }

    public File getColCubeFile()
    {
        return new File(stsDirectory, colCubeFilename);
    }

    public File getSliceCubeFile()
    {
        return new File(stsDirectory, sliceCubeFilename);
    }

    public void setColCubeFilename(String colCubeFilename)
    {
        this.colCubeFilename = colCubeFilename;
    }

    public String getColCubeFilename()
    {
        return colCubeFilename;
    }

    public void setSliceCubeFilename(String sliceCubeFilename)
    {
        this.sliceCubeFilename = sliceCubeFilename;
    }

    public String getSliceCubeFilename()
    {
        return sliceCubeFilename;
    }

    public String getCreateHeaderFilename(String stemname)
    {
        if (headerFilename == null)
            headerFilename = createFilename(headerFormat, stemname);
        return headerFilename;
    }

    /**
     * returns a header filename with this stemname
     *
     * @param stemname
     * @return
     */

    public String createHeaderFilename(String stemname)
    {
        return createFilename(headerFormat, stemname);
    }

    static public String createHeaderFilename(String groupname, String stemname)
    {
        return createFilename(groupname, headerFormat, stemname);
    }

    public String getCreateHeaderFilename()
    {
        return getCreateHeaderFilename(stemname);
    }

    public String getCreateAttributeFilename(String stemname)
    {
        if (attributeFilename == null)
            attributeFilename = createFilename(attributeFormat, stemname);
        return attributeFilename;
    }

    public String getCreateAttributeFilename()
    {
        return getCreateAttributeFilename(stemname);
    }

    public String createFilename(String format)
    {
        return getGroupname() + "." + format + "." + stemname;
    }

    public String createFilename(String format, String name)
    {
        return getGroupname() + "." + format + "." + name;
    }

    static public String createRowByteFilename(String group, String name)
    {
        return createFilename(group, byteFormat, inline, name);
    }

    static public String createColByteFilename(String group, String name)
    {
        return createFilename(group, byteFormat, xline, name);
    }

    static public String createSliceByteFilename(String group, String name)
    {
        return createFilename(group, byteFormat, trace, name);
    }

    static public String createRowFloatFilename(String group, String name)
    {
        return createFilename(group, floatFormat, inline, name);
    }

    public String getCreateRowFloatFilename(String stemname)
    {
        rowFloatFilename = createFilename(floatFormat + "-" + inline, stemname);
        return rowFloatFilename;
    }


    static public String createFilename(String group, String format, String stemname)
    {
        return group + "." + format + "." + stemname;
    }

    static public String createFilename(String group, String format, String planeType, String stemname)
    {
        return group + "." + format + "-" + planeType + "." + stemname;
    }

    public String getGroupname()
    {
        new StsMessage(currentModel.win3d, StsMessage.WARNING, "Developer: need to implement getGroup() in calling class " + StsToolkit.getSimpleClassname(this));
        return "none";
    }

    public File getHeaderFile(String stemname)
    {
        return new File(stsDirectory, getCreateHeaderFilename(stemname));
    }

    public boolean isVolumeRegular()
    {
        return isRegular;
    }

    public String getName()
    {
        return stemname;
    }

    public void setName(String name)
    {
        this.stemname = name;
        super.setName(name);
        if (isPersistent())
        {
            dbFieldChanged("stemname", name);
            dbFieldChanged("name", name);
        }
        if (currentModel != null) currentModel.refreshObjectPanel(this);
    }

    public int getCursorPlaneIndex(int dirNo, float dirCoordinate)
    {
        switch (dirNo)
        {
            case XDIR: // crossline direction
                return getNearestColCoor(dirCoordinate);
            case YDIR: // inline direction
                return getNearestRowCoor(dirCoordinate);
            case ZDIR:
                return getNearestSliceCoor(dirCoordinate);

            default:
                return -1;
        }
    }
/*
	public void setFilenames()
	{
		setFilenames(stemname);
	}
*/


    static final boolean debugDirectoryAlias = false;

       /** returns true if pathname exists and file is found ok; otherwise false */
    public boolean checkStsDirectoryForFilename(String filename) throws FileNotFoundException
    {
        if(debugDirectoryAlias) System.out.print("old pathname: " + stsDirectory + filename);
        if(filename == null) return false;
        if(StsFile.constructor(stsDirectory, filename).exists())
        {
            if(debugDirectoryAlias) System.out.println("  new pathname: IDENTICAL");
            return true;
        }
        String aliasDirectory = currentModel.getProject().getDirectoryAlias(stsDirectory);
        while(aliasDirectory != null)
        {
        	if(debugDirectoryAlias) System.out.println("  found alias directory for: " + stsDirectory + " : " + aliasDirectory);
        	stsDirectory = aliasDirectory;
        	if(StsFile.constructor(stsDirectory, filename).exists())
        	{
        		if(debugDirectoryAlias) System.out.println("  found filename " + filename + " in alias directory: " + aliasDirectory);
        		return true;
        	}
        	else
        	{
        		aliasDirectory = currentModel.getProject().getDirectoryAlias(stsDirectory);
        		continue;
        	}
        }
        // either an alias directory doesn't exist or file not found in alias directory: ask user for new alias directory
        while(aliasDirectory == null)
        {
            Frame frame = new Frame(); // win3d may not be realized yet so the dialog is thrown in the top left instead of in the center of screen
            aliasDirectory = StsToolkit.findDirectoryForFile(frame, stsDirectory, filename);
            if(aliasDirectory == null)
            {
                if(debugDirectoryAlias) System.out.println("  User canceled.");
                new StsMessage(frame, StsMessage.ERROR, "User canceled directory selection. Unable to load " + filename);
                break;
                //throw new FileNotFoundException("User canceled directory selection. Unable to load " + filename);
            }
            else if(!(new File(aliasDirectory)).isDirectory())
            {
                if(debugDirectoryAlias) System.out.println("  User supplied directory does not exist.");
                new StsMessage(frame,StsMessage.ERROR, "User supplied directory ("+aliasDirectory+") does not exist for file: " + filename);
                aliasDirectory = null;
                //throw new FileNotFoundException("User supplied directory ("+aliasDirectory+") does not exist for file: " + filename);
            }
            else if(!StsFile.constructor(aliasDirectory, filename).exists())
            {
                if(debugDirectoryAlias) System.out.println("  File does not exist in user supplied directory (" + aliasDirectory + ").");
                new StsMessage(frame,StsMessage.ERROR, "File does not exist in user supplied directory (" + aliasDirectory + ").");
                aliasDirectory = null;
                //throw new FileNotFoundException("File does not exist in user supplied directory (" + aliasDirectory + ") does not exist.");
            }
            else
            {
                currentModel.getProject().addDirectoryAlias(stsDirectory, aliasDirectory);
                stsDirectory = aliasDirectory;
                StsMessageFiles.infoMessage("New supplied project path is " + stsDirectory);
                if(debugDirectoryAlias) System.out.println("  New pathname: " + stsDirectory + filename);
            }
        }
        // waste of time to change field since the old directory is forever a part of this db and
        // the alias redefinition will have to take place everytime
        // dbFieldChanged("stsDirectory", stsDirectory);
        return false;
    }

    public void createVolumeFilenames()
    {
        rowCubeFilename = getVolumeFilename(inline);
        colCubeFilename = getVolumeFilename(xline);
        sliceCubeFilename = getVolumeFilename(trace);
    }

    public void createVolumeFilenames(String stemname)
    {
        rowCubeFilename = getVolumeFilename(inline, stemname);
        colCubeFilename = getVolumeFilename(xline, stemname);
        sliceCubeFilename = getVolumeFilename(trace, stemname);
    }

    public String createFloatRowVolumeFilename(String stemname)
    {
//        if (!isDataFloat) return null;
        rowFloatFilename = createFilename(floatFormat + "-" + inline, stemname);
        return rowFloatFilename;
    }

    public String getVolumeFilename(String planeType)
    {
        return createFilename(byteFormat + "-" + planeType);
    }

    public String getVolumeFilename(String planeType, String stemname)
    {
        return createFilename(byteFormat + "-" + planeType, stemname);
    }

    public void setFilenames(String name)
    {
        createVolumeFilenames(name);
        getCreateAttributeFilename(name);
        if (isDataFloat) getCreateRowFloatFilename(name);
    }

    /*
        public boolean isDataFloat()
        {
            if (floatRowFilename == null)
            {
                isDataFloat = false;
            }
            else
            {
                File file = new File(stsDirectory, floatRowFilename);
                isDataFloat = file.exists();
            }
            return isDataFloat;
        }
    */
    public void setIsDataFloat(boolean isFloat)
    { isDataFloat = isFloat; }

    public boolean getIsDataFloat() { return isDataFloat; }

    public void deleteExistingFiles()
    {
        if (rowCubeFilename != null) deleteExistingFile(stsDirectory, rowCubeFilename);
        if (colCubeFilename != null) deleteExistingFile(stsDirectory, colCubeFilename);
        if (sliceCubeFilename != null) deleteExistingFile(stsDirectory, sliceCubeFilename);
        if (rowFloatFilename != null) deleteExistingFile(stsDirectory, rowFloatFilename);
        if (attributeFilename != null) deleteExistingFile(stsDirectory, attributeFilename);
    }

    public void deleteExistingFile(String directory, String filename)
    {
        File file = new File(directory, filename);
        if (file.exists())
        {
            file.delete();
        }
    }

    public boolean writeHeaderFile()
    {
        return writeHeaderFile("");
    }

    public boolean writeHeaderFile(String extension)
    {
        try
        {
            String filename = getCreateHeaderFilename(stemname) + extension;
            StsParameterFile.writeObjectFields(stsDirectory + filename, this, StsSeismicBoundingBox.class, StsMainObject.class);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "writeHeaderFile()", e);
            return false;
        }
    }

    public File getHeaderFile()
    {
        return new File(stsDirectory + headerFilename);
    }

    public boolean derivedFileOK(String derivedFileStemname)
    {
        String filename;
        File file;
        try
        {
            // check that timestamp on derived file is later than source file
            file = getHeaderFile();
            long volumeTime = file.lastModified();
            filename = getCreateHeaderFilename(derivedFileStemname);
            file = new File(stsDirectory + filename);
            if (!file.exists()) return false;
            long derivedTime = file.lastModified();
            if (derivedTime < volumeTime) return false;

            // check existence and approximate length of the 3 derived volume files by comparing  to size of source row file
            file = getRowCubeFile();
            long size = file.length();

            filename = getVolumeFilename(inline, derivedFileStemname);
            file = new File(stsDirectory, filename);
            if (!file.exists()) return false;
            if (Math.abs(file.length() - size) > 0.05 * size) return false;

            filename = getVolumeFilename(xline, derivedFileStemname);
            file = new File(stsDirectory, filename);
            if (!file.exists()) return false;
            if (Math.abs(file.length() - size) > 0.05 * size) return false;

            filename = getVolumeFilename(trace, derivedFileStemname);
            file = new File(stsDirectory, filename);
            if (!file.exists()) return false;
            if (Math.abs(file.length() - size) > 0.05 * size) return false;

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicBoundingBox.derivedFileOK() failed.", e, StsException.WARNING);
            return false;
        }
    }

    /** returns a floatBuffer if available, otherwise a byteBuffer */
    public StsMappedBuffer createMappedRowBuffer(String mode)
    {
        if (rowFloatFilename != null)
        {
            StsMappedBuffer buffer = StsMappedFloatBuffer.constructor(stsDirectory, rowFloatFilename, mode);
            if (buffer != null) return buffer;
        }
        return StsMappedByteBuffer.constructor(stsDirectory, rowCubeFilename, mode, dataMin, dataMax);
    }

    public StsMappedFloatBuffer createMappedFloatRowBuffer()
    {
        return StsMappedFloatBuffer.openReadWrite(stsDirectory, rowFloatFilename);
    }

    public StsMappedByteBuffer createMappedByteRowBuffer(String mode)
    {
        return StsMappedByteBuffer.constructor(stsDirectory, rowCubeFilename, mode, dataMin, dataMax);
    }

    public boolean getIsXLineCCW() { return isXLineCCW; }

    public void setIsXLineCCW(boolean isXLineCCW) { this.isXLineCCW = isXLineCCW; }

    public boolean getIsArchived() { return isArchived; }

    public void setIsArchived(boolean value) { isArchived = value; }

    public float getVerticalScalar() { return verticalScalar; }

    public float getHorizontalScalar() { return horizontalScalar; }

    public void setVerticalScalar(float vertScale) { verticalScalar = vertScale; }

    public void setHorizontalScalar(float horzScale) { horizontalScalar = horzScale; }

    public int getAttributeIndex(String name)
    {
        for (int n = 0; n < this.nAttributes; n++)
            if (name.equals(this.attributeNames[n]))
                return n;
        return -1;
    }

    public int getNAttributes()
    {
        return attributeNames.length;
    }

    public String[] getAttributes()
    {
        return attributeNames;
    }

    public String[] getTimeAttributes()
    {
        return StsSEGYFormat.getTimeAttributes(attributeNames);
    }

    public String[] getDistanceAttributes()
    {
        return StsSEGYFormat.getDistanceAttributes(attributeNames);
    }

    public String[] getOtherAttributes()
    {
        return StsSEGYFormat.getOtherAttributes(attributeNames);
    }

    public String[] getDatumAttributes()
    {
        return StsSEGYFormat.getDatumAttributes(attributeNames);
    }

    public String[] getVelocityAttributes()
    {
        return StsSEGYFormat.getVelocityAttributes(attributeNames);
    }

    public String[] getCoordinateAttributes()
    {
        return StsSEGYFormat.getCoordinateAttributes(attributeNames);
    }

    // TODO These methods should really be defined in an interface used by the export routines.

    /**
     * Used for export of volume. Returns a byte plane for given row plane number.
     * Must be overriden in subclasses being exported.
     */
    public byte[] readRowPlaneByteData(int nPlane)
    {
        StsException.systemError(this, "readRowPlaneByteData(nPlane)", "has not been implemented in subclass.");
        return null;
    }

    /**
     * Used for export of volume. Returns a float plane for given row plane number.
     * Must be overriden in subclasses being exported.
     */
    public float[] readRowPlaneFloatData(int nPlane)
    {
        StsException.systemError(this, "System error: readRowPlaneFloatData(nPlane)", "has not been implemented in subclass.");
        return null;
    }

    public boolean setupReadRowFloatBlocks()
    {
        StsException.systemError(this, "System error: setupReadRowFloatBlocks()", "has not been implemented in subclass.");
        return false;
    }

    public boolean setupRowByteBlocks()
    {
        StsException.systemError(this, "System error: setupRowByteBlocks()", "has not been implemented in subclass.");
        return false;
    }

    public void setVelocityTypeString(String typeString)
    {
        StsException.systemError(this, "System error: setVelocityTypeString(String)", "has not been implemented in subclass.");
    }

    public String getVelocityTypeString()
    {
        StsException.systemError(this, "System error: getVelocityTypeString()", "has not been implemented in subclass.");
        return null;
    }

    public void setVelocityType(byte type)
    {
        StsException.systemError(this, "System error: setVelocityType(byte)", "has not been implemented in subclass.");
    }

    public byte getVelocityType()
    {
        StsException.systemError(this, "System error: getVelocityType()", "has not been implemented in subclass.");
        return 0;
    }

    public float getMinDepthAtTime(float time)
    {
        StsException.systemError(this, "System error: getMinDepthAtTime(float)", "has not been implemented in subclass.");
        return 0;
    }

    public float getMaxDepthAtTime(float time)
    {
        StsException.systemError(this, "System error: getMaxDepthAtTime(float)", "has not been implemented in subclass.");
        return 0;
    }

    /**
     * Returns an [n][3] array of CDP x, y, and number by CDP.
     *
     * @return double[][]
     */
    public double[][] getLineXYCDPs(int row)
    {
        StsException.systemError(this, "System error: getLineXYCDPs()", "has not been implemented in subclass.");
        return null;
    }

    public int getColumnsInRow(int row)
    {
        return this.getColMax();
    }

    public boolean hasExportableData(int row)
    {
        return false;
    }

    /** override this method for 2d where nCroppedSlices depends on row */
    public int getNSlices(int row)
    {
        return nSlices;
    }

    public int getNSamples(int row)
    {
        return getNSlices(row);
    }
    public StsSeismicBoundingBox getLineBoundingBox(int row)
    {
        StsException.notImplemented(this, "getLineBoundingBox");
        return null;
    }

    /** some subclasses have multiple colorscales for different displays.  Implement addColorscale to add colorscales to a list. */
    public void addColorscale(StsColorscale colorscale)
    {
        StsException.notImplemented(this, "addColorscale");
    }

    public void setLinesIndexRanges(StsRotatedGridBoundingBox[] lines)
    {
        if (lines == null) return;
        nRows = lines.length;
        nCols = lines[0].nCols;
        for (int n = 1; n < lines.length; n++)
            nCols = Math.max(lines[n].nCols, nCols);
    }

    public boolean analyzeBinaryHdr(StsProgressPanel progressPanel) //throws IOException
    {
        return segyData.analyzeBinaryHdr(progressPanel);
}

    public boolean isAValidFileSize()
   {
        return segyData.isAValidFileSize();
   }

    public void initializeSegyIO()
    {
        segyIO = new StsSegyIO(segyData);
    }

    class KnownPoint
    {
        float inline = nullValue;
        float crossline = nullValue;
        double x = doubleNullValue;
        double y = doubleNullValue;

        KnownPoint()
        {
        }
    }

    public void setUseKnownPoints(boolean val)
    {
        useKnownPoints = val;
    }

    public void initializeKnownPoints()
    {
        knownPoints = new KnownPoint[3];
        for (int n = 0; n < 3; n++)
            knownPoints[n] = new KnownPoint();
        return;
    }

    public boolean checkKnownPointsOK()
    {
        useKnownPoints = knownPointsOK();
        return useKnownPoints;
    }

    private boolean knownPointsOK()
    {
        for (int i = 0; i < 3; i++)
        {
            if (knownPoints[i].x == doubleNullValue)
                return false;
            if (knownPoints[i].y == doubleNullValue)
                return false;
            if (knownPoints[i].inline == nullValue)
                return false;
            if (knownPoints[i].crossline == nullValue)
                return false;
        }
        return true;
    }

    public double getPoint1X()
    {
        return knownPoints[0].x;
    }

    public double getPoint1Y()
    {
        return knownPoints[0].y;
    }

    public float getPoint1Inline()
    {
        return knownPoints[0].inline;
    }

    public float getPoint1Xline()
    {
        return knownPoints[0].crossline;
    }

    public double getPoint2X()
    {
        return knownPoints[1].x;
    }

    public double getPoint2Y()
    {
        return knownPoints[1].y;
    }

    public float getPoint2Inline()
    {
        return knownPoints[1].inline;
    }

    public float getPoint2Xline()
    {
        return knownPoints[1].crossline;
    }

    public double getPoint3X()
    {
        return knownPoints[2].x;
    }

    public double getPoint3Y()
    {
        return knownPoints[2].y;
    }

    public float getPoint3Inline()
    {
        return knownPoints[2].inline;
    }

    public float getPoint3Xline()
    {
        return knownPoints[2].crossline;
    }

    public void setPoint1X(double value)
    {
        knownPoints[0].x = value;
    }

    public void setPoint1Y(double value)
    {
        knownPoints[0].y = value;
    }

    public void setPoint1Inline(float value)
    {
        knownPoints[0].inline = value;
    }

    public void setPoint1Xline(float value)
    {
        knownPoints[0].crossline = value;
    }

    public void setPoint2X(double value)
    {
        knownPoints[1].x = value;
    }

    public void setPoint2Y(double value)
    {
        knownPoints[1].y = value;
    }

    public void setPoint2Inline(float value)
    {
        knownPoints[1].inline = value;
    }

    public void setPoint2Xline(float value)
    {
        knownPoints[1].crossline = value;
    }

    public void setPoint3X(double value)
    {
        knownPoints[2].x = value;
    }

    public void setPoint3Y(double value)
    {
        knownPoints[2].y = value;
    }

    public void setPoint3Inline(float value)
    {
        knownPoints[2].inline = value;
    }

    public void setPoint3Xline(float value)
    {
        knownPoints[2].crossline = value;
    }
/*
    public boolean analyzeGeometryWithKnownPoints(StsSEGYFormat segyFormat)
    {
        StsSEGYFormat.TraceHeader[] surveyTraces = new StsSEGYFormat.TraceHeader[3];
        for (int n = 0; n < 3; n++)
        {
            surveyTraces[n] = segyFormat.constructTraceHeader();
            surveyTraces[n].x = (float) knownPoints[n].x * horizontalScalar;
            surveyTraces[n].y = (float) knownPoints[n].y * horizontalScalar;
            surveyTraces[n].xLine = knownPoints[n].crossline;
            surveyTraces[n].iLine = knownPoints[n].inline;
        }

        // compute rotation angle and bin spacings in line and xline directions
        if (!analyzeAngle(surveyTraces[0], surveyTraces[1], surveyTraces[2]))
        {
            return false;
        }
        initializeRange();

        xOrigin = knownPoints[0].x * horizontalScalar;
        yOrigin = knownPoints[0].y * horizontalScalar;

        return true;
    }
*/
    /**
     * Given an origin trace (arbitrary) and two other traces, compute the line and xline spacings
     * and rotation angle of grid.  Any 3 traces can by supplied as long as they are not in a straight line.
     * The line and xline numbering increments already need to have been computed.
     * We assume lines and xlines are orthogonal. We have two vectors from origin to each of the two points.
     * Each vector has components in the line and xline directions in which we know the index difference in
     * each component but not the index spacing.  We have 3 unknowns: yInc, xInc, and lineAngle
     * (the angle from the +X axis to the +line direction.  Using Pythagorean theorem, we solve two equations
     * for the two spacing unknowns.  The angle is then computed using the dot and cross products between
     * the same vector in the rotated and unrotated coordinate systems.
     */
    public boolean analyzeAngle(StsSEGYFormat.TraceHeader originTrace, StsSEGYFormat.TraceHeader traceA, StsSEGYFormat.TraceHeader traceB)
    {
        double dLineA = 0.0, dXLineA = 0.0;
        double dLineB = 0.0, dXLineB = 0.0;
        double dLineSqA = 0.0, dXLineSqA = 0.0;
        double dLineSqB = 0.0, dXLineSqB = 0.0;
        double lengthSqA = 0.0, lengthSqB = 0.0;
        double dxA = 0.0, dyA = 0.0, dxB = 0.0, dyB = 0.0;
        boolean traceAok = false, traceBok = false;

        if (traceA != null)
        {
            if (rowNumInc != 0.0f)
            {
                dLineA = (traceA.iLine - originTrace.iLine) / rowNumInc;
                dLineSqA = dLineA * dLineA;
            }
            if (colNumInc != 0.0f)
            {
                dXLineA = (traceA.xLine - originTrace.xLine) / colNumInc;
                dXLineSqA = dXLineA * dXLineA;
            }
            dxA = traceA.x - originTrace.x;
            dyA = traceA.y - originTrace.y;
            lengthSqA = dxA * dxA + dyA * dyA;
            if (lengthSqA > 1.0f) traceAok = true;
        }
        if (traceB != null)
        {
            if (rowNumInc != 0.0f)
            {
                dLineB = (traceB.iLine - originTrace.iLine) / rowNumInc;
                dLineSqB = dLineB * dLineB;
            }
            if (colNumInc != 0.0f)
            {
                dXLineB = (traceB.xLine - originTrace.xLine) / colNumInc;
                dXLineSqB = dXLineB * dXLineB;
            }
            dxB = traceB.x - originTrace.x;
            dyB = traceB.y - originTrace.y;
            lengthSqB = dxB * dxB + dyB * dyB;
            if (lengthSqB > 1.0f) traceBok = true;
        }

        if (!traceAok && !traceBok) return false;
        if (traceAok && traceBok)
        {
            // check if angle between two vectors is too close ( less than 1 degree)
            double crossIndexes = dXLineA * dLineB - dLineA * dXLineB;
            double dotIndexes = dXLineA * dXLineB + dLineA * dLineB;
            double indexAngle = StsMath.atan2(dotIndexes, crossIndexes);
            if (Math.abs(indexAngle) < 1) return false;

            double crossSq = dXLineSqA * dLineSqB - dXLineSqB * dLineSqA;

            //        if(Math.abs(crossSq) < 0.1*lengthSqA) return false;

            double yIncSq = (lengthSqB * dXLineSqA - lengthSqA * dXLineSqB) / crossSq;
            double yInc = Math.sqrt(yIncSq) * horizontalScalar;

            double xIncSq;
            if (dXLineSqA > dXLineSqB)
                xIncSq = (lengthSqA - yIncSq * dLineSqA) / dXLineSqA;
            else
                xIncSq = (lengthSqB - yIncSq * dLineSqB) / dXLineSqB;
            double xInc = Math.sqrt(xIncSq) * horizontalScalar;

            // check if xLines are 90 deg CCW from lines; if not change sign of yInc
            double crossLines = dxA * dyB - dyA * dxB;

            // compute angle from same vector in rotated coordinate system to vector in unrotated coordinate system
            // this is the rotation angle from +X in unrotated coordinates to +Line direction
            double cosL = (xInc * dXLineA * dxA + yInc * dLineA * dyA);
            double sinL = (xInc * dXLineA * dyA - yInc * dLineA * dxA);
            double lineAngle = StsMath.atan2(cosL, sinL);
            this.yInc = (float) yInc;
            this.xInc = (float) xInc;
            this.angle = (float) lineAngle;

            isXLineCCW = (crossIndexes * crossLines >= 0);
        }
        else if (traceAok)
        {
            double lengthA = Math.sqrt(lengthSqA);
            if (dLineA != 0.0)
            {
                this.yInc = (float) (lengthA / dLineA);
                this.xInc = 0.0f;
                this.angle = (float) StsMath.atan2(dxA, dyA);
            }
            else if (dXLineA != 0.0)
            {
                this.yInc = 0.0f;
                this.xInc = (float) (lengthA / dXLineA);
                this.angle = (float) StsMath.atan2(dxA, dyA);
            }
            else
                return false;
        }
        else if (traceBok)
        {
            double lengthB = Math.sqrt(lengthSqB);
            if (dLineB != 0.0)
            {
                this.yInc = (float) (lengthB / dLineB);
                this.xInc = 0.0f;
                this.angle = (float) StsMath.atan2(dxB, dyB);
            }
            else if (dXLineB != 0.0)
            {
                this.yInc = 0.0f;
                this.xInc = (float) (lengthB / dXLineB);
                this.angle = (float) StsMath.atan2(dxB, dyB);
            }
            else
                return false;
        }
        if(yInc <= 0.0f || xInc <= 0.0f || Float.isNaN(xInc) || Float.isNaN(yInc)) return false;
        initializeRange(originTrace);
        if (!isXLineCCW)
        {
            moveOriginToLL();
        }
        return true;
    }

       /**
     * Flip row numbering if xLine is CW from line as
     * origin has been moved to lower-left from upper-left.
     */
    public void checkFlipRowNumOrder()
    {
        if(isXLineCCW) return;
        float temp = rowNumMin;
        rowNumMin = rowNumMax;
        rowNumMax = temp;
        rowNumInc = -rowNumInc;
    }

    public void initializeRange(StsSEGYFormat.TraceHeader originTrace)
    {
        xMin = xInc*Math.round((colNumMin - originTrace.xLine)/colNumInc);
        xMax = xMin + (nCols - 1) * xInc;
        yMin = yInc*Math.round((rowNumMin - originTrace.iLine)/rowNumInc);
        yMax = (nRows - 1) * yInc;
    }

   // public void setCdpInterval(float cdpInterval) {this.cdpInterval = cdpInterval;}

    //public float getCdpInterval() {return cdpInterval;}

    public void setOverrideSampleFormat(int format)
    {
        segyData.setOverrideSampleFormat(format);
    }

    public int getSampleFormat()
    {
        return segyData.getSampleFormat();
    }

    public int getHeaderSampleFormat()
    {
        return segyData.getHeaderSampleFormat();
    }

    public void setTextHeader(byte[] textHeader)
    {
        segyData.setTextHeader(textHeader);
    }

    public void setBinaryHeader(byte[] binaryHeader)
    {
        segyData.setBinaryHeader(binaryHeader);
    }

    public byte[] getTraceHeader()
    {
        return segyData.getTraceHeader();
    }

    public void setTraceHeader(byte[] traceHeader)
    {
        segyData.setTraceHeader(traceHeader);
    }

    public int getTextHeaderSize()
    {
        return segyData.getTextHeaderSize();
    }

    public void setTextHeaderSize(int textHeaderSize)
    {
        segyData.setTextHeaderSize(textHeaderSize);
    }

    public int getBinaryHeaderSize()
    {
        return segyData.getBinaryHeaderSize();
    }

    public void setBinaryHeaderSize(int binaryHeaderSize)
    {
        segyData.setBinaryHeaderSize(binaryHeaderSize);
    }

    public int getTraceHeaderSize()
    {
        return segyData.getTraceHeaderSize();
    }

    public void setTraceHeaderSize(int traceHeaderSize)
    {
        segyData.setTraceHeaderSize(traceHeaderSize);
    }

    public int getBytesPerSample()
    {
        return segyData.getBytesPerSample();
    }

    public void setBytesPerSample(int bytesPerSample)
    {
        segyData.setBytesPerSample(bytesPerSample);
    }

    public int getFileHeaderSize()
    {
        return segyData.getFileHeaderSize();
    }

    public void setFileHeaderSize(int fileHeaderSize)
    {
        segyData.setFileHeaderSize(fileHeaderSize);
    }

    public int getBytesPerTrace()
    {
        return segyData.getBytesPerTrace();
    }

    public void setBytesPerTrace(int bytesPerTrace)
    {
        segyData.setBytesPerTrace(bytesPerTrace);
    }

    public long getSegyFileSize()
    {
        return segyData.getSegyFileSize();
    }

    public void setSegyFileSize(long segyFileSize)
    {
        segyData.setSegyFileSize(segyFileSize);
    }

    public int getNTotalTraces()
    {
        if (segyData == null) return 0;
        return segyData.getNTotalTraces();
    }

    public long getNTotalSamples()
    {
        if (segyData == null) return 0;
        return segyData.getNTotalSamples();
    }

    public void setNTotalTraces(int nTotalTraces)
    {
        segyData.setNTotalTraces(nTotalTraces);
    }

    public float getSampleSpacing()
    {
        return segyData.getSampleSpacing();
    }

    public void setSampleSpacing(float sampleSpacing)
    {
        segyData.setSampleSpacing(sampleSpacing);
    }

    public void setClipString(String value)
    {
        for(int n = 0; n < 3; n++)
            if(value == CLIP_STRINGS[n])
            {
                clipType = (byte)n;
                return;
            }
        clipType = 0;
    }

    public String getClipString() { return CLIP_STRINGS[clipType]; }

    public void setHasUserNull(boolean value)
    {
        hasUserNull = value;
        if(!hasUserNull) userNull = nullValue;
    }

    public boolean getHasUserNull() { return hasUserNull; }

    public void setUserNull(float value) { userNull = value; }
    public float getUserNull() { return userNull; }

    public byte[] getTextHeader()
    {
        return segyData.getTextHeader();
    }

    public byte[] getBinaryHeader()
    {
        return segyData.getBinaryHeader();
    }

    public byte[] getTraceHeaderBinary(int nTrace)
    {
        return segyData.getTraceHeaderBinary(nTrace);
    }

    public StsSEGYFormat getSegyFormat()
    {
        return segyData.getSegyFormat();
    }

    public void setSegyFormat(StsSEGYFormat segyFormat)
    {
        segyData.setSegyFormat(segyFormat);
    }

    public void readFileHeader() throws IOException
    {
        segyData.readFileHeader();
    }


    public void setIsLittleEndian(boolean isLittleEndian)
    {
        segyData.setIsLittleEndian(isLittleEndian);
    }

    public boolean getIsLittleEndian()
    {
        return segyData.getIsLittleEndian();
    }

    public int getTotalTraces()
    {
        return segyData.nTotalTraces;
    }

    public int getNTraces()
    {
        return segyData.nTotalTraces;
    }

    public double getBinaryHeaderValue(StsSEGYFormatRec rec)
    {
        return segyData.getBinaryHeaderValue(rec);
    }
    public String[] getAttributeNames() { return attributeNames; }

    static public float calcDiskRequired(StsSeismicBoundingBox[] volumes)
    {
        float size = 0;
        if (volumes == null) return size;
        for (int i = 0; i < volumes.length; i++)
            size += volumes[i].getSegyFileSize() - StsSEGYFormat.defaultHeaderSize;
        return size / 1000000.0f;
    }

	public int analyzeTraces(double scanPercentage, StsProgressPanel progressPanel) // throws IOException
	{
		try
		{
			initializeSegyIO();

            int result = segyIO.analyzeTraces(segyData.randomAccessSegyFile, (float)scanPercentage, progressPanel);
			dataMin = segyIO.dataMin;
			dataMax = segyIO.dataMax;
            if( result == StsSegyIO.ANALYSIS_OK)
                calculateHistogram(segyIO.histogramSamples, segyIO.nHistogramSamples);
			return result;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "analyzeTraces", e);
			return StsAnalyzeTraces.ANALYSIS_UNKNOWN_EXCEPTION;
		}
	}

    //TODO use the StsToolkit calculateHistogram instead (need to add scale, scaleOffset to argument options)
    public void calculateHistogram(float[] histogramSamples, int nHistogramSamples)
    {
        StsException.systemDebug(this, "calculateHistogram(histogramSamples, nHistogramSamples");
        if(nHistogramSamples == 0) return;
        int[] dataCnt = new int[255];
        initializeScaling();
        dataHist = new float[255];
        float scale = 100.0f;
        if(hasUserNull)
        {
            int nNonNullSamples = 0;
            for (int n = 0; n < histogramSamples.length; n++)
            {
                float value = histogramSamples[n];
                if(value == userNull) continue;
                int scaledInt = scaleFloatToInt(value);
                dataCnt[scaledInt]++;
                nNonNullSamples++;
            }
            if(nNonNullSamples == 0) return;
            scale /= nNonNullSamples;
        }
        else
        {
            for (int n = 0; n < histogramSamples.length; n++)
            {
                int scaledInt = scaleFloatToInt(histogramSamples[n]);
                dataCnt[scaledInt]++;
            }
        }
        for (int i = 0; i < 255; i++)
            dataHist[i] = scale*dataCnt[i];
    }

    /**
     * ScaleType is used in conversion of floats to bytes.
     * SCALE_VOLUME: dataMin and dataMax define the float range.
     * SCALE_NORMALIZED:  float data is already normalized, 0 to 1.
     * SCALE_TRACE: Auto-scale each trace using min and max of trace
     * The default is SCALE_VOLUME. Override this in subclasses for NORMALIZED and TRACE.
     */

    public byte getScaleType()
    {
        return SCALE_VOLUME;
    }

    public void checkScaling()
    {
        byte scaleType = getScaleType();
        if (scaleType != SCALE_VOLUME) return;
        initializeScaling();
    }

    public void initializeScaling()
    {
        scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
        scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, dataMin);
    }

    public final int scaleFloatToInt(float value)
    {
        int scaledInt = Math.round(scale * value + scaleOffset);
        return StsMath.minMax(scaledInt, 0, 254);
    }

    public boolean scaleFloatsToBytes(float[] floats, byte[] bytes)
    {
        byte scaleType = getScaleType();
        if (scaleType == SCALE_VOLUME)
            return scaleToBytes(floats, bytes);
        else if (scaleType == SCALE_NORMALIZED)
            return scaleNormalizedFloatsToBytes(floats, bytes);
        else if(scaleType == SCALE_TRACE_NORMALIZED)
            return scaleNormalizedTraceFloatsToBytes(floats, bytes);
        else
            return normalizeAndScaleFloatsToBytes(floats, bytes);
    }

    private boolean scaleToBytes(float[] floats, byte[] bytes)
    {
        if( clipType == CLIP_NONE && !hasUserNull)
        {
            for(int n = 0; n < floats.length; n++)
            {
                float value = floats[n];
                bytes[n] = StsMath.floatToUnsignedByte254WithScale(value, scale, dataMin, dataMax);
            }
            return true;
        }
        for(int n = 0; n < floats.length; n++)
        {
            float value = floats[n];
            if(hasUserNull && value == userNull)
            {
                bytes[n] = StsParameters.nullByte;
                continue;
            }
            if( clipType == CLIP_NULL )
            {
                if(value < dataMin || value > dataMax)
                {
                    bytes[n] = StsParameters.nullByte;
                    continue;
                }
            }
            else if( clipType == CLIP_MINMAX )
            {
                if(value < dataMin)
                {
                    bytes[n] = 0;
                    continue;
                }
                else if(value > dataMax)
                {
                    bytes[n] = (byte)254;
                    continue;
                }
            }
            bytes[n] = StsMath.floatToUnsignedByte254WithScale(value, scale, scaleOffset);
        }
        return true;
    }

    static public boolean scaleNormalizedFloatsToBytes(float[] floatData, byte[] byteData)
    {
        if (floatData == null) return false;
        for (int n = 0; n < floatData.length; n++)
        {
            if (floatData[n] == StsParameters.nullValue)
                byteData[n] = StsParameters.nullByte;
            else
                byteData[n] = StsMath.unsignedIntToUnsignedByte254((int) (254 * floatData[n]));
        }
        return true;
    }

    static public boolean scaleNormalizedTraceFloatsToBytes(float[] floatData, byte[] byteData)
    {
        if (floatData == null) return false;
        for (int n = 0; n < floatData.length; n++)
        {
            if (floatData[n] == StsParameters.nullValue)
                byteData[n] = StsParameters.nullByte;
            else
                byteData[n] = StsMath.unsignedIntToUnsignedByte254((int) (127 + 127 * floatData[n]));
        }
        return true;
    }

    static public boolean normalizeAndScaleFloatsToBytes(float[] floats, byte[] bytes)
    {
        if (floats == null) return false;
        int nValues = floats.length;
        float maxAmplitude = 0.0f;
        for (int n = 0; n < nValues; n++)
        {
            if (floats[n] != StsParameters.nullValue)
                maxAmplitude = Math.max(maxAmplitude, Math.abs(floats[n]));
        }
        float scale;
        if (maxAmplitude == 0.0f)
            scale = 0.0f;
        else
            scale = 127 / maxAmplitude;

        for (int n = 0; n < nValues; n++)
        {
            if (floats[n] == StsParameters.nullValue)
                bytes[n] = StsParameters.nullByte;
            else
                bytes[n] = (byte) (127 + scale * floats[n]);
        }
        return true;
    }

    public boolean scaleFloatsToBytes(FloatBuffer floats, ByteBuffer bytes, int length)
    {
        byte scaleType = getScaleType();
        if (scaleType == SCALE_VOLUME)
            return scaleToBytes(floats, bytes, length);
        else if (scaleType == SCALE_NORMALIZED)
            return scaleNormalizedFloatsToBytes(floats, bytes, length);
        else if(scaleType == SCALE_TRACE_NORMALIZED)
            return scaleNormalizedTraceFloatsToBytes(floats, bytes, length);
        else
            return normalizeAndScaleFloatsToBytes(floats, bytes, length);
    }

    private boolean scaleToBytes(FloatBuffer floats, ByteBuffer bytes, int length)
    {
        try
        {
            if(clipType == CLIP_NONE && !hasUserNull)
            {

                for(int n = 0; n < length; n++)
                {
                    float value = floats.get();
                    bytes.put(StsMath.floatToUnsignedByte254WithScale(value, scale, dataMin, dataMax));
                }
                return true;
            }
            else
            {
                for(int n = 0; n < length; n++)
                {
                    float value = floats.get();
                    if(value == StsParameters.nullValue || hasUserNull && value == userNull)
                        bytes.put(StsParameters.nullByte);
                    else if(clipType == CLIP_NONE)
                    {
                        bytes.put(StsMath.floatToUnsignedByte254WithScale(value, scale, dataMin, dataMax));
                    }
                    else if( clipType == CLIP_NULL)
                    {

                        if(value < dataMin || value > dataMax)
                            bytes.put(StsParameters.nullByte);
                        else
                            bytes.put(StsMath.floatToUnsignedByte254WithScale(value, scale, dataMin, dataMax));
                    }
                    else if( clipType == CLIP_MINMAX)
                    {
                        bytes.put(StsMath.floatToUnsignedByte254WithScale(value, scale, dataMin, dataMax));
                    }
                }
                return true;
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsMath.class, "floatToUnsignedByte254WithScale", e);
            return false;
        }
    }

    static public boolean scaleNormalizedFloatsToBytes(FloatBuffer floats, ByteBuffer bytes, int length)
    {
        if (floats == null) return false;
        for (int n = 0; n < length; n++)
        {
            float value = floats.get();
            if (value == StsParameters.nullValue)
                bytes.put(StsParameters.nullByte);
            else
                bytes.put(StsMath.unsignedIntToUnsignedByte254((int) (254 * value)));
        }
        return true;
    }
    static public boolean scaleNormalizedTraceFloatsToBytes(FloatBuffer floats, ByteBuffer bytes, int length)
     {
         if (floats == null) return false;
         for (int n = 0; n < length; n++)
         {
             float value = floats.get();
             if (value == StsParameters.nullValue)
                 bytes.put(StsParameters.nullByte);
             else
                 bytes.put(StsMath.unsignedIntToUnsignedByte254((int) (127 + 127 * value)));
         }
         return true;
     }

    static public boolean normalizeAndScaleFloatsToBytes(FloatBuffer floats, ByteBuffer bytes, int length)
    {
        if (floats == null) return false;
        float maxAmplitude = 0.0f;
        for (int n = 0; n < length; n++)
        {
            float value = floats.get();
            if (value != StsParameters.nullValue)
                maxAmplitude = Math.max(maxAmplitude, Math.abs(value));
        }
        float scale;
        if (maxAmplitude == 0.0f)
            scale = 0.0f;
        else
            scale = 127 / maxAmplitude;

        floats.rewind();
        for (int n = 0; n < length; n++)
        {
            float value = floats.get();
            if (value == StsParameters.nullValue)
                bytes.put(StsParameters.nullByte);
            else
                bytes.put((byte) (127 + scale * value));
        }
        return true;
    }

    public boolean isVolumeScaling()
    {

        return getScaleType() == SCALE_VOLUME;
    }

    /** status handling procedure.
     *  Volume starts out with status STATUS_OK (5).
     *
     *  On fileFormat wizardStep, analyzeBinaryHdr is first run.  If sampleSpacing is <= 0 or nSamples <= 0,
     *  status is set to HEADER_BAD (4).  If these are ok, status is set to HEADER_OK (6).
     *  Next on fileFormat wizardStep, checkValidFileSize is run.  If fileSize is not consistent with fileHeaderSize,
     *  nSamples, sampleFormat, and traceHeaderSize, status is set to FILE_BAD (3).  If ok, status is set to FILE_OK (7).
     *  If the user changes any parameters on fileFormat wizardStep, status of selected files is reset to STATUS_OK.
     *  On the traceDefinition wizardStep, analyzeTraces is run on each file whose status is not bad for any reason (bad is status < 5).
     *  If successful, status is set to TRACES_OK (8); if not, status is set to TRACES_BAD (2).
     *  Next on traceDefinition wizardStep, analyzeGrid is run for each file whose status is not bad.
     *  If successful, status is set to GEOMETRY_OK; if not, status is set to GEOMETRY_BAD.
     *  For 3d, an additional step is run to analye the xy geometry; if ok, status is set to GRID_OK; otherwise GRID_BAD.
     *  If any parameters are changed on traceDefinition wizardStep, status is reset to FILE_OK if status is <= FILE_BAD.
     *
     *
     * @return
     */
    public boolean isStatusOK() { return statusString == STATUS_OK_STR; }
    public void resetStatus() { status = STATUS_OK; }

    /** trace records have changed, so reset status to successful trace analysis status (STATUS_TRACES_OK)
     *  if status indicates volume failed grid or geometry analysis header&trace analysis was successful (status <= STATUS_GRID_BAD)
     */
    public void resetTraceStatus() { if(status <= STATUS_GRID_BAD || status > STATUS_TRACES_OK) status = STATUS_TRACES_OK; }

    public void setDataRangeType(String dataRangeType) { this.dataRangeType = dataRangeType; }
    public void adjustDataRange(float[] floats)
    {
        if (floats == null || floats.length == 0) return;
        for (int n = 0; n < floats.length; n++)
        {
            if (floats[n] < dataMin)
                dataMin = floats[n];
            else if (floats[n] > dataMax)
                dataMax = floats[n];
        }
    }

    public boolean hasOutputFiles()
	{
        File file;
        file = new File(stsDirectory + rowCubeFilename);
		if(!file.exists()) return false;
		file = new File(stsDirectory + colCubeFilename);
		if(!file.exists()) return false;
		file = new File(stsDirectory + sliceCubeFilename);
		if(!file.exists()) return false;
        if(rowFloatFilename != null)
        {
            file = new File(stsDirectory + rowFloatFilename);
            if(!file.exists()) return false;
        }
        return true;
	}

    public void constructHistogram(float[] data, int nValues)
    {
        clearHistogram();
        float scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
        float scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, dataMin);
        for(int n = 0; n < nValues; n++)
        {
            int index = Math.round(scale*data[n] + scaleOffset);
            accumulateHistogram(index);
        }
        calculateHistogram();
    }

    public void accumulateHistogram(float data, float scale, float scaleOffset)
    {
        int index = Math.round(scale*data + scaleOffset);
        accumulateHistogram(index);
    }

    public void accumulateHistogram(float[] data, float scale, float scaleOffset)
    {
        for(int n = 0; n < data.length; n++)
            accumulateHistogram(data[n], scale, scaleOffset);
    }

    public void accumulateHistogram(byte byteIndex)
    {
        int histIndex = StsMath.signedByteToUnsignedInt(byteIndex);
        if (histIndex > 254) histIndex = 254;
        dataCnt[histIndex]++;
        ttlHistogramSamples++;
    }

    public void accumulateHistogram(int histIndex)
    {
        histIndex = StsMath.minMax(histIndex, 0, 254);
        dataCnt[histIndex]++;
        ttlHistogramSamples++;
    }

    public void calculateHistogram()
    {
        System.out.println("calculateHistogram");
        for (int i = 0; i < 255; i++)
        {
            dataHist[i] = ((float) dataCnt[i] / (float) ttlHistogramSamples) * 100.0f;
        }
    }

    public float[] getHistogram() { return dataHist; }

    public void clearHistogram()
    {
        System.out.println("ClearHistogram");
        if (dataHist == null) return;
        for (int i = 0; i < 255; i++)
        {
            dataCnt[i] = 0;
            dataHist[i] = 0.0f;
        }
        ttlHistogramSamples = 0;
    }

    public String getRowColNumLabel()
    {
        return new String("inline range: " + rowNumMin + " " + rowNumMax + " crossline range: " + colNumMin + " " + colNumMax);
    }

    public String getVolumeTypeString() { return StsParameters.getSampleTypeString(volumeType); }
    public void setVolumeTypeString(byte volumeType) { this.volumeType = volumeType; }

    // for StsXYSurfaceGridable interface compatability
/*
	public float interpolateBilinearZ(float[] xyz, boolean computeIfNull, boolean setPoint)
	{
		return 0.0f;
	}

	public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint)
	{
		return 0.0f;
	}

	public StsPoint getPoint(int row, int col)
	{
		return new StsPoint(0.0, 0.0, 0.0);
	}

	public float[] getXYZ(int row, int col)
	{
		return new float[] {0.0f, 0.0f, 0.0f};
	}

	public float[] getNormal(float rowF, float colF)
	{
		return new float[] {0.0f, 0.0f, 1.0f};
	}

	public float[] getNormal(int row, int col)
	{
		return new float[] {0.0f, 0.0f, 1.0f};
	}

	public float computePointZ(int row, int col)
	{
		return 0.0f;
	}

	public float getComputePointZ(int row, int col)
	{
		return 0.0f;
	}

    public void checkConstructGridNormals() { }
*/
}

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.lang.reflect.*;
import java.nio.*;
import java.util.*;

/**
 * StsSurface represents a basic surface which defines only a grid with no associated topological
 * elements such as faults, zones, etc.  @see com.Sts.DBTypes.StsModelSurface for a surface which
 * contains these topological elements.
 * <p/>
 * StsSurface includes a rotational angle unlike most gridding formats, so for gridding formats which
 * do not define an angle, it must be specified if other than zero.
 * <p/>
 * Regardless of how the imported surface is defined, when StsSurface is constructed it has a
 * consistent geometry:  xOrigin and yOrigin (doubles) define the origin for the grid.  Rows are
 * aligned in increasing local x coordinate direction and columns are aligned in increasing local
 * y direction.  Z values for the grid are ordered sequentially along rows; i.e., the first point
 * is row 0 col 0, the second point is row 0 col 1, etc.  The angle for the grid is counter-clockwise
 * from the globable X-axis to the local X-axis and is between 0 and 360.
 * The xMin and yMin are offsets in the local coordinate system from the project
 * origin to this grid origin.
 * See <a href="surfaceCoorSystem.html">Surface Coor System</a> for details.
 */
public class StsSurface extends StsRotatedGridBoundingBox implements StsSelectable, StsXYSurfaceGridable, StsTreeObjectI,
    ItemListener, StsTextureSurfaceFace, StsDistanceTransformFace, StsCultureDisplayable
{
    /** indicates the surfaces has nulls requiring special display routines */
    protected boolean hasNulls = false;
    /**
     * The standard null value; input surfaces which have a different null value
     * have those values converted to this value
     */
    public float nullZValue = StsParameters.nullValue;
    /** The minimum row number which has contains a non-null value */
    protected int iMinNonNull;
    /** The maximum row number which has contains a non-null value */
    protected int iMaxNonNull;
    /** The minimum column number which has contains a non-null value */
    protected int jMinNonNull;
    /** The maximum column number which has contains a non-null value */
    protected int jMaxNonNull;
    protected String prefixString = "";
    /** The color of this surface described in red-green-blue-alpha values from 0 to 1. */
    protected StsColor stsColor;
    /** z values in array [nRows][nCols]; time, depth, or seismic depth */
    protected float[][] pointsZ = null;
    /** adjusted z values: depth if pointsZ are time or seismic depth; not used if pointsZ are depth */
    protected float[][] adjPointsZ = null;
    /** null types in array [nRows][nCols] see listing of types */
    public byte[][] pointsNull = null;
    /** polygons which define regions to be nulled as points on a fault */
    protected StsFaultPolygon[] faultPolygons;
    /** grid lines are displayed: black if fill is true and surface color if not. */
    protected boolean displayGrid = false;
    /** filled and lighted surface is displayed in the surface color. */
    protected boolean displayFill = true;
    /** display property color on surface. */
    protected boolean displayProperty = false;
    /** pixel mode */
    protected boolean isPixelMode = false;
    /** shader */
    transient protected boolean contourColor = false; // doesn't do much with subsampled textures
    /** surface is currently lighted */
    protected boolean lighting = true;
    /** surface is currently smoothly lighted */
    protected boolean smoothNormals = true;
    /** marker associated with this surface */
    protected StsMarker marker = null;
    /** zDomain specifies whether surface is originally in time or depth */
    protected byte zDomainOriginal = StsParameters.TD_NONE;
    /** zDomain(s) supported; set to two domains (time & depth or approxDepth & depth) if a velocity model is available and has been applied. */
    protected byte zDomainSupported = StsParameters.TD_NONE;
    /** If surface has both depth and time, zMin is time min */
    protected float depthMin;
    /** If surface has both depth and time, zMax is time max */
    protected float depthMax;
    /** Turn on/off screen door transparency */
    protected boolean transparency = false;
    /** List of surface attributes associated with this surface */
    StsObjectRefList surfaceAttributes = null;
    /* name of current surfaceTexture object*/
    String surfaceTextureName = null;
    /** min property value for current texture */
    transient float propertyMin = 0.0f;
    /** max property value for current texture */
    transient float propertyMax = 1.0f;
    /**
     * zDomain currently being displayed. Changing domains requires building new display lists and textures;
     * in which case zDomainDisplayed is set to none, display() method deletes display lists, rebuilds displays
     * for current project zDomain and sets ZDomainDisplayed to this zDomain.
     */
    transient protected byte zDomainDisplayed = StsParameters.TD_NONE;

    transient float relativeRotationAngle = 0.0f;
    transient double rowXInc, rowYInc, colXInc, colYInc;
    transient protected float[][][] normals = null;
    /** offset of surface vertically from its original position */
    transient float offset = 0.0f;

    transient public float[][] weights = null; // scratch array used in surface interpolation
    transient public float[][] dZdX, dZdY; // scratch arrays used in gradient-weighted interpolation
    transient boolean initialized = false;

    /** array holding texture data: fetched from textureDisplayable object and deleted after transfer to OpenGL */
    //    transient byte[] textureData = null;
    /** color display list number: fetched from textureDisplayable object */
    transient int textureColorListNum = 0;
    transient boolean colorListChanged = true;

    /** shader used for displays */
    transient int shader = StsJOGLShader.NONE;

    /** Texture selected for display on this surface. */
    transient StsSurfaceTexture newSurfaceTexture = null;
    /** Current texture being isVisible on surface. */
    transient StsSurfaceTexture currentSurfaceTexture = null;
    /** list of available textureDisplayable objects for surface */
    transient StsSurfaceTexture[] surfaceTextureList = null;

    /** Tiles on which texture is generated */
    transient StsTextureTiles textureTiles = null;
    /** Texture has been changed: replace texture subImage2D with new texture */
    transient boolean textureChanged = true;
    /** Surface geometry has been changed: reconstruct it */
    transient boolean geometryChanged = true;
    boolean useShader = false;
    /** Texture color display list */
    transient int colorDisplayListNum = 0;
    /** Display lists currently being used for surface geometry */
    transient boolean usingDisplayLists = false;
    /**
     * Construct a regular grid. If the number of rows and/or cols exceed
     * maxRowCols, construct a fine grid decimated to meet this criteria.
     * rowStride is the increment between rows in coarse grid.
     */
    transient int rowStride = 1;
    /**
     * Construct a regular grid. If the number of rows and/or cols exceed
     * maxRowCols, construct a fine grid decimated to meet this criteria.
     * colStride is the increment between cols in coarse grid.
     */
    transient int colStride = 1;
    /** indicates we've built a coarse grid in addition to the fine grid */
    transient boolean hasCoarseGrid = false;

    /** flag (set any time) which deletes displayLists on next draw attempt */
    //    transient protected boolean deleteDisplayLists = false;
    /** display list number for surface fill */
    transient private int surfDisplayListNum = 0;
    /** display list number for surface grid lines */
    transient private int gridDisplayListNum = 0; // display list number (> 0)

    // transient StsList tStrips;
    transient byte[] units = null;

    /**
     * Used by runnable threads, so we have associated it with this class so it is accessible
     * from methods which run a thread and update this progress bar
     */
    transient StsProgressBarDialog progressBarDialog = null;
    /** default texture on surface */
    transient StsSurfaceTexture surfaceDefaultTexture = new SurfaceDefaultTexture(this);
    /** depth texture on surface */
    transient SurfaceDepthTexture surfaceDepthTexture = null;
    /** gradient texture on surface */
    transient SurfaceGradTexture surfaceGradTexture = null;
    /** textures for any seismic volumes displayable on this surface.  Texture data is filled only when needed. */
    transient StsSurfaceTexture[] seismicTextures;
    /** current curvature attribute */
    transient public StsSurfaceCurvatureAttribute curvatureAttribute = null;
    /** transient curvature attributes on surface built in this session (must be saved to be retained in project) */
    transient StsSurfaceCurvatureAttribute[] transientCurvatureAttributes = new StsSurfaceCurvatureAttribute[0];
    /** textures for any surface attributes displayable on this surface.  Texture data is filled only when needed. */
    transient StsSurfaceTexture[] attributeTextures;
    /** Maximum gradient range for gradTexture display */
    transient float gradRangeMax;
    /** Maximum curvature range for gradCurvature display */
    //transient float curvRangeMin;
    /** Minimum curvature range for gradCurvature display */
    //transient float curvRangeMax;
    /** indicates that minCorrelFilter needs to be used in distanceTransform interpolation operation */
    transient boolean applyMinCorrelFilter = false;

    transient static StsTimer timer;

    // Maintain the native units
    transient protected byte nativeHorizontalUnits = StsParameters.DIST_NONE;
    transient protected byte nativeVerticalUnits = StsParameters.DIST_NONE;

    // Convenience copies of useful flags
    static public final byte SURF_PNT = StsParameters.SURF_PNT;
    static public final byte SURF_BOUNDARY = StsParameters.SURF_BOUNDARY;

    static public final byte SURF_GAP = StsParameters.SURF_GAP;
    static public final byte SURF_GAP_SET = StsParameters.SURF_GAP_SET;
    static public final byte SURF_GAP_FILLED = StsParameters.SURF_GAP_FILLED;
    //	static public final byte NULL_GAP_OR_BOUNDARY = StsParameters.NULL_GAP_OR_BOUNDARY;
    static public final byte SURF_GAP_NOT_FILLED = StsParameters.SURF_GAP_NOT_FILLED;
    static private final int maxNullIndex = 7;

    static public final byte PLUS = StsParameters.PLUS;
    static public final byte MINUS = StsParameters.MINUS;
    static public final float roundOff = StsParameters.roundOff;
    static public final int ROW = StsParameters.ROW;
    static public final int COL = StsParameters.COL;

    static public final float halfSqrt2f = StsParameters.halfSqrt2f;

    static final boolean debug = false;
    static boolean fillGaps = false; // mainDebug only: fill in gaps in surface

    /** Convenience copies of TStrip flags */
    static public final byte STRIP_INVALID = StsParameters.STRIP_INVALID;
    static public final byte STRIP_BOTH = StsParameters.STRIP_BOTH;
    static public final byte STRIP_BOT = StsParameters.STRIP_BOT;
    static public final byte STRIP_TOP = StsParameters.STRIP_TOP;

    static public final int nullInteger = StsParameters.nullInteger;
    static public final float largeFloat = StsParameters.largeFloat;

    // Class version ID: rerun serial version when the class is changed
    static final long serialVersionUID = -5294107068658877877L;
    static protected StsObjectPanel objectPanel = null;
    static boolean useTrueCellCenteredDrawing = true;

   /** Max number of rows or cols in coarse level grid */
    static private int maxDisplayRowOrCols = 300;

    // State/type flags
    static public final byte NONE = 0;
    static public final byte IMPORTED = 1;
    static public final byte MODEL = 2;
    static public final byte CHECKED = 3;
    static public final byte INITIALIZED = 4;
    static public final byte GAPPED = 5;
    static public final byte BUILT = 6;
    static public final byte INTERPOLATED = 7;
    // display fields
    static StsComboBoxFieldBean surfaceTextureListBean;
    static StsFloatFieldBean propertyMaxBean, propertyMinBean;
    static StsFloatFieldBean offsetBean;
    static StsDateFieldBean bornField = new StsDateFieldBean(StsSurface.class, "bornDate", "Born Date:");
    static StsDateFieldBean deathField = new StsDateFieldBean(StsSurface.class, "deathDate", "Death Date:");
    static public StsFieldBean[] displayFields = null;

    static StsEditableColorscaleFieldBean colorscaleBean;
    static public StsFieldBean[] propertyFields = null;

    static public final String seismicGrp = "grid-seismic";
    static public final String zmapGrp = "grid";
    static public final String[] fileGroups = new String[] { seismicGrp, zmapGrp };

    static boolean debugDistTrans = false;
    static int debugDistTransRow = 41;
    static int debugDistTransCol = 591;
    /**
     * Development flag:  we are currently always using textures, even for plain surface display,
     * So the non-texture code is currently being skipped.  We left it in in case we need to reuse at some point.
     */
    static final boolean useTextures = true;
    static final byte nullByte = StsParameters.nullByte;

    /** constructor for DB */
    public StsSurface()
    {
    }

    public StsSurface(boolean persistent)
    {
        super(persistent);
    }

    static public StsSurface constructSurface(String name, StsColor stsColor, byte type, StsRotatedGridBoundingBox boundingBox, byte zDomain, float[][] pointsZ)
    {
        try
        {
            return new StsSurface(name, stsColor, type, boundingBox, zDomain, pointsZ);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurface.construct() failed.", e,
                StsException.WARNING);
            return null;
        }
    }

    /**
     * constructor for surface that may have nulls called by static public constructor
     *
     * @param name     name of this surface
     * @param stsColor color in rgba 0 to 1 float values
     * @param type     kind of surface: IMPORTED or MODEL are current types
     */
    private StsSurface(String name, StsColor stsColor, byte type, StsRotatedGridBoundingBox boundingBox, byte zDomain, float[][] pointsZ) throws StsException
    {
        this(name, stsColor, type, boundingBox.nCols, boundingBox.nRows, boundingBox.xOrigin, boundingBox.yOrigin, boundingBox.xInc, boundingBox.yInc,
            boundingBox.xMin, boundingBox.yMin, boundingBox.angle, pointsZ, false, StsSurface.nullValue, zDomain, null);
    }

    static public StsSurface constructSurface(String name, StsColor stsColor,
                                              byte type,
                                              int nCols, int nRows,
                                              double xOrigin, double yOrigin,
                                              float xInc, float yInc,
                                              float xMin, float yMin,
                                              float angle, float[][] pointsZ,
                                              boolean hasNulls,
                                              float nullZValue, byte zDomain, byte vUnits, byte hUnits,
                                              StsProgressPanel progressPanel)
    {
        try
        {
            //nativeVerticalUnits = vUnits;
            //nativeHorizontalUnits = hUnits;

            return new StsSurface(name, stsColor, type, nCols, nRows, xOrigin,
                yOrigin, xInc, yInc, xMin, yMin, angle,
                pointsZ, hasNulls, nullZValue, zDomain, progressPanel); //, hUnits, vUnits);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.construc() failed.", e,
                StsException.WARNING);
            return null;
        }
    }

    static public StsSurface constructSurface(String name, StsColor stsColor, byte type, int nCols, int nRows,
                                              double xOrigin, double yOrigin, float xInc, float yInc,
                                              float xMin, float yMin, float angle, float[][] pointsZ,
                                              boolean hasNulls, float nullZValue, byte zDomain,
                                              StsProgressPanel progressPanel)
    {
        try
        {
            return new StsSurface(name, stsColor, type, nCols, nRows, xOrigin,
                yOrigin, xInc, yInc, xMin, yMin, angle,
                pointsZ, hasNulls, nullZValue, zDomain, progressPanel);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.construc() failed.", e,
                StsException.WARNING);
            return null;
        }
    }

    /** constructor for surface that may have nulls */
    private StsSurface(String name, StsColor stsColor, byte type, int nCols,
                       int nRows,
                       double xOrigin, double yOrigin, float xInc, float yInc,
                       float xMin, float yMin, float angle,
                       float[][] pointsZ, boolean hasNulls, float nullZValue,
                       byte zDomain, StsProgressPanel progressPanel) throws StsException
    {
        initialize(name, stsColor, type, nCols, nRows, xOrigin, yOrigin, xInc,
            yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue,
            zDomain, zDomain, progressPanel);
        initSurfaceTextureList();

    }

    private StsSurface(String name, StsColor stsColor, byte type, int nCols,
                       int nRows,
                       double xOrigin, double yOrigin, float xInc, float yInc,
                       float xMin, float yMin, float angle,
                       float[][] pointsZ, boolean hasNulls, float nullZValue,
                       byte zDomain, StsProgressPanel progressPanel, byte hUnits, byte vUnits) throws StsException
    {
        initialize(name, stsColor, type, nCols, nRows, xOrigin, yOrigin, xInc,
            yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue,
            zDomain, zDomain, progressPanel);
        initSurfaceTextureList();

        nativeVerticalUnits = vUnits;
        nativeHorizontalUnits = hUnits;
    }

    protected StsSurface(boolean persistent, String name)
    {
        super(false);
        setName(name);
    }

    static public StsSurface constructTempSurface(String name)
    {
        return new StsSurface(false, name);
    }

    public boolean initialize(String name, StsColor stsColor, byte type,
                              int nCols, int nRows,
                              double xOrigin, double yOrigin, float xInc,
                              float yInc, float xMin, float yMin, float angle,
                              float[][] pointsZ, boolean hasNulls,
                              float nullZValue, byte zDomain, byte zDomainSupported, StsProgressPanel progressPanel)
    {
        // save inputs
        setName(name);
        this.stsColor = stsColor;
        this.type = type;

        displayFill = getSurfaceClass().getDefaultDisplayFill();
        displayGrid = getSurfaceClass().getDefaultDisplayGrid();

        return initialize(nCols, nRows, xOrigin, yOrigin, xInc, yInc, xMin,
            yMin, angle, pointsZ, hasNulls, nullZValue, zDomain, zDomainSupported, progressPanel);
    }

    public boolean initialize(StsGridDefinition gridDef, byte zDomain)
    {
        return initialize(gridDef.getNCols(), gridDef.getNRows(), gridDef.getXOrigin(), gridDef.getYOrigin(),
            gridDef.getXInc(), gridDef.getYInc(), gridDef.getXMin(), gridDef.getYMin(),
            gridDef.getAngle(), null, true, StsParameters.nullValue, zDomain, zDomainSupported, null);
    }

    public boolean initialize(int nCols, int nRows, double xOrigin, double yOrigin, float xInc, float yInc, float xMin, float yMin,
                              float angle, float[][] pointsZ, boolean hasNulls, float nullZValue,
                              byte zDomain, byte zDomainSupported, StsProgressPanel progressPanel)
    {
        this.surfaceAttributes = StsObjectRefList.constructor(10, 1, "surfaceAttributes", this);
        this.surfaceTextureName = getName();
        this.nCols = nCols;
        this.nRows = nRows;
        computeRowColStride();
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
        this.xInc = xInc;
        this.yInc = yInc;
        this.xMin = xMin;
        this.yMin = yMin;
        this.angle = angle;
        this.pointsZ = pointsZ;
        this.hasNulls = hasNulls;
        this.nullZValue = nullZValue;
        this.zDomainOriginal = zDomain;
        this.zDomainSupported = zDomainSupported;

        if (!checkGrid(currentModel.win3d))
        {
            return false;
        }

        setAngle();
        checkResetOrigin();
        if (pointsZ != null)
        {
            initializePointsZ();
            setOffsetRange();
        }
        //        setLighting();

        StsProject project = currentModel.getProject();

        if (!project.isOriginSet())
        {
            project.setOriginAndAngle(this.xOrigin, this.yOrigin, this.angle);
            this.xMin = 0.0f;
            this.yMin = 0.0f;
            this.xMax = (nCols - 1) * xInc;
            this.yMax = (nRows - 1) * yInc;
        }
        else
        {
            float[] relativeOrigin = project.computeRelativeOrigin(this.xOrigin, this.yOrigin);
            if (relativeOrigin != null)
            {
                this.xMin += relativeOrigin[0];
                this.yMin += relativeOrigin[1];
                this.xMax = xMin + (nCols - 1) * xInc;
                this.yMax = yMin + (nRows - 1) * yInc;
            }
        }

        if (pointsZ != null)
        {
            //            StsStatusArea status = StsStatusArea.getStatusArea();
            StsMessageFiles.logMessage("Processing grid...");

            interpolateDistanceTransform(progressPanel);
            //            StsInterpolationRadiusWeighted.constructor(this).interpolate();
            //          StsSurfaceInterpolationWeightedPlane.getInstance(this).interpolateSurface();
            constructGrid();
        }
        //initSurfaceTextureList();

        //		project.adjustRange(getRelativeRange());

        //        System.out.println("StsGrid: nulls ? " + hasNulls + " nx " + nX + " ny " + nY + "\n" +
        //        					"Xmin " + xMin + " Xmax " + xMax + " Xinc " + xInc + "\n" +
        //        					"Ymin " + yMin + " Ymax " + yMax + " Yinc " + yInc);
        return true;
    }

    public void computeRowColStride()
    {
        if (nRows > maxDisplayRowOrCols)
        {
            rowStride = 2;
            hasCoarseGrid = true;
        }
        if (nCols > maxDisplayRowOrCols)
        {
            colStride = 2;
            hasCoarseGrid = true;
        }
        // mainDebug
        /*
           rowStride = 1;
           colStride = 1;
           hasCoarseGrid = true;
        */
    }

    static public StsSurface constructThicknessGrid(StsSurface topSurface, StsSurface botSurface, StsGridDefinition gridDef, byte zDomain)
    {
        String name = "thickness_" + topSurface.getName() + "-" + botSurface.getName();
        StsSurface thicknessSurface = new StsSurface(false, name);
        thicknessSurface.initialize(gridDef, zDomain);
        thicknessSurface.initSurfaceTextureList();
        thicknessSurface.initializePoints();
        thicknessSurface.constructThicknessGrid(topSurface, botSurface);

        // for debugging: remove when finished
        thicknessSurface.setStsColor(topSurface.getStsColor());
        thicknessSurface.addToModel();

        return thicknessSurface;
    }

    private void constructThicknessGrid(StsSurface topSurface, StsSurface botSurface)
    {
        int nRows = topSurface.nRows;
        int nCols = topSurface.nCols;
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                if (topSurface.getPointNull(row, col) == SURF_PNT && botSurface.getPointNull(row, col) == SURF_PNT)
                {
                    setPointNull(row, col, SURF_PNT);
                    float topZ = topSurface.getPointZ(row, col);
                    float botZ = botSurface.getPointZ(row, col);
                    // jbw why ????
                    if (botZ < topZ )
                    {
                        botZ = topZ;
                        botSurface.setPointFilled(row, col, botZ, SURF_PNT);
                    }
                    float dZ = botZ - topZ;
                    setPointFilled(row, col, dZ, SURF_PNT);
                }
                else
                    setPointNull(row, col, SURF_BOUNDARY);
            }
        }
        interpolateDistanceTransform();
        //        StsSurfaceInterpolationRadiusWeighted.constructor(this, NULL_BOUNDARY).interpolate();
    }


    public boolean isDisplayable()
    {
        return isPersistent();
    }

    /** Add a group and surface attributes to this surface */
    public void addSurfaceAttributes(StsObject[] surfaceAttributes)
    {
        this.surfaceAttributes.add(surfaceAttributes);
        //    	initSurfaceTextureDisplayableList();
    }

    public void initSurfaceTextureList()
    {
        surfaceTextureList = getSurfaceTextureList();
        //        if(textureDisplayableListBean == null) return;
        //        textureDisplayableListBean.setListItems(textureDisplayableList);
        if (surfaceTextureName != null)
        {
            for (int i = 0; i < surfaceTextureList.length; i++)
            {
                if (surfaceTextureName.equals(surfaceTextureList[i].getName()))
                {
                    setSurfaceTexture(surfaceTextureList[i], null);
                    if (surfaceTextureListBean != null)
                        surfaceTextureListBean.setSelectedItem(newSurfaceTexture);
                    break;
                }
            }
        }
    }

    public StsSurfaceTexture[] getSurfaceTextureList()
    {
        constructTextureDisplayableList();
        return surfaceTextureList;
    }
    
    public StsSurfaceCurvatureAttribute createCurvatureAttribute(byte curveType, int filterSize, StsProgressPanel progressPanel)
    {
    	curvatureAttribute = new StsSurfaceCurvatureAttribute(this, curveType, filterSize);
    	//doesn't exist so calculate
    	boolean runOk = curvatureAttribute.createAttribute(progressPanel);
    	if (runOk)
    		setNewCurvatureAttribute(curvatureAttribute);
    	transientCurvatureAttributes = (StsSurfaceCurvatureAttribute[])StsMath.arrayAddElement(transientCurvatureAttributes, curvatureAttribute);
        return curvatureAttribute;
    }
    /** set the current surfaceTextureDisplayable, but not if its the already this one */
    public void setNewCurvatureAttribute(StsSurfaceCurvatureAttribute curvatureAttribute)
    {
        setNewSurfaceTexture(curvatureAttribute.surfaceTexture);
    }
    
    public void saveCurvatureAttribute()
    {
    	if (curvatureAttribute == null) return;
        if(curvatureAttribute.isPersistent()) return;
        curvatureAttribute.setFiltered(true);
        curvatureAttribute.addToModel();
        String filename = new String("grid" + "." + "bin"  + "." + this.name);
        String attFilename = new String(filename + "." + curvatureAttribute.getName() + ".0");
        curvatureAttribute.writeBinaryFile(attFilename);

        surfaceAttributes.add(curvatureAttribute);
        transientCurvatureAttributes = (StsSurfaceCurvatureAttribute[])StsMath.arrayDeleteElement(transientCurvatureAttributes, curvatureAttribute);
    }

    /** called by StsSurface and StsModelSurface to build/rebuild the textureDisplayableList */
    public void constructTextureDisplayableList()
    {
        if (debug) StsException.systemDebug(this, "initSurfaceTextureList", " building list for: " + getName());
        StsSeismicVolume[] seismicVolumes = (StsSeismicVolume[]) currentModel.getCastObjectList(StsSeismicVolume.class);
        int nSurfaceAttributes = getNumberSurfaceAttributes();
        int nListItems = seismicVolumes.length + nSurfaceAttributes + 5 + transientCurvatureAttributes.length;
        int n = 0;
        surfaceTextureList = new StsSurfaceTexture[nListItems];
        // defaultSurfaceTexture = new StsSurfaceTexture(this);
        surfaceTextureList[n++] = surfaceDefaultTexture;
        if (surfaceDepthTexture == null)
            surfaceDepthTexture = new SurfaceDepthTexture(this);
        surfaceTextureList[n++] = surfaceDepthTexture;
        if (surfaceGradTexture == null)
            surfaceGradTexture = new SurfaceGradTexture(this);
        surfaceTextureList[n++] = surfaceGradTexture;
//        if (surfaceCurvTexture == null)
//        	surfaceCurvTexture = new SurfaceCurvTexture(this);
//        if (surfaceCurvTexture != null)
//        	surfaceTextureList[n++] = surfaceCurvTexture;
        StsSurfaceTexture correlCoefSurfaceTexture = getCorrelCoefTexture();
        if (correlCoefSurfaceTexture != null)
            surfaceTextureList[n++] = correlCoefSurfaceTexture;

        //saved surface attributes
        for (int nn = 0; nn < nSurfaceAttributes; nn++)
        {
            StsSurfaceAttribute surfaceAttribute = (StsSurfaceAttribute) surfaceAttributes.getElement(nn);
            surfaceTextureList[n++] = surfaceAttribute.getSurfaceTexture();
        }
    	//curvature textures created in this session
    	for (StsSurfaceCurvatureAttribute surfaceCurvature : transientCurvatureAttributes)
        {
    		surfaceTextureList[n++] = surfaceCurvature.surfaceTexture;
        }

        for (int nn = 0; nn < seismicVolumes.length; nn++)
            surfaceTextureList[n++] = seismicVolumes[nn].getSurfaceTexture(this);

        surfaceTextureList = (StsSurfaceTexture[]) StsMath.trimArray(surfaceTextureList, n);

//        for (StsSurfaceTexture surfaceTextureDisplayable : surfaceTextureList)
//            surfaceTextureDisplayable.initializeColorscaleActionListener();
    }

    public StsSurfaceCurvatureAttribute[] getCurvatureAttributes()
    {
        Object[] attributes = surfaceAttributes.getTrimmedList();
        StsSurfaceCurvatureAttribute[] atts = new  StsSurfaceCurvatureAttribute[attributes.length + transientCurvatureAttributes.length];
        int cnt = 0;
        for(int i=0; i<attributes.length; i++)
        {
            if(attributes[i] instanceof StsSurfaceCurvatureAttribute)
                atts[cnt++] = (StsSurfaceCurvatureAttribute)attributes[i];
        }
        for(int i=0; i<transientCurvatureAttributes.length; i++)
            atts[cnt++] = transientCurvatureAttributes[i];
        return (StsSurfaceCurvatureAttribute[])StsMath.trimArray(atts, cnt);
        //return (StsSurfaceCurvatureAttribute[])StsMath.arrayAddArray(attributes, transientCurvatureAttributes, StsSurfaceCurvatureAttribute.class);
    }

    private int getNumberSurfaceAttributes()
    {
        if(surfaceAttributes == null) return 0;
        return surfaceAttributes.getSize();
    }

    protected StsSurfaceTexture getCorrelCoefTexture()
    {
        StsHorpickClass horpickClass = (StsHorpickClass) currentModel.getStsClass(StsHorpick.class);
        if (horpickClass == null) return null;
        StsHorpick horpick = horpickClass.getHorpickWithSurface(this);
        if (horpick == null) return null;
        return horpick.getCorrelCoefsSurfaceTexture(this);
    }

    protected StsSurfaceTexture getSurfaceAttributeTextureDisplayable(int index)
    {
        if (surfaceAttributes == null) return null;
        StsSurfaceAttribute attribute = (StsSurfaceAttribute) surfaceAttributes.getElement(index);
        StsSurfaceTexture attributeTexture = attribute.getSurfaceTexture();
        attributeTextures = checkAddTexture(attributeTexture, attributeTextures);
        return attributeTexture;
    }

    private StsSurfaceTexture[] checkAddTexture(StsSurfaceTexture surfaceTexture, StsSurfaceTexture[] currentSurfaceTextures)
    {
        for (StsSurfaceTexture currentTexture : currentSurfaceTextures)
        {
            if (currentTexture == surfaceTexture) return currentSurfaceTextures;
        }
        currentSurfaceTextures = (StsSurfaceTexture[]) StsMath.arrayAddElement(currentSurfaceTextures, surfaceTexture);
        return currentSurfaceTextures;
    }

    public float getCultureZ(float x, float y) { return getZ(new StsPoint(x, y)); }

    public boolean isPlanar() { return false; }

    /*
    public float getCorrelCoefMinFilter()
    {
     StsHorpickClass horpickClass = (StsHorpickClass) currentModel.getStsClass(StsHorpick.class);
     StsHorpick horpick = horpickClass.getHorpickWithSurface(this);
     if(horpick == null) return 0.0f;
     return horpick.getMinCorrelFilter();
    }

    public void setCorrelCoefMinFilter(float minCorrelFilter)
    {
     StsHorpickClass horpickClass = (StsHorpickClass) currentModel.getStsClass(StsHorpick.class);
     StsHorpick horpick = horpickClass.getHorpickWithSurface(this);
     if(horpick == null) return;
     horpick.setMinCorrelFilter(minCorrelFilter);
    }
    */
    public void applyMinCorrelFilter(boolean applyMinCorrelFilter, StsSurfaceTexture correlCoefTexture)
    {
        if (correlCoefTexture == null)
        {
            return;
        }
        this.applyMinCorrelFilter = applyMinCorrelFilter;
        // setIsVisible(true);
        // setDisplayFill(true);
        if (newSurfaceTexture != correlCoefTexture)
        {
            setNewSurfaceTexture(newSurfaceTexture);
        }
        interpolateDistanceTransform();
    }

    public void checkInitializeGrid(byte nullType, boolean hasNulls)
    {
        if (pointsZ != null && pointsNull != null)
        {
            return;
        }
        initializeGrid(nullType, hasNulls);
    }

    public void initializeGrid(byte nullType, boolean hasNulls)
    {
        pointsZ = new float[nRows][nCols];
        setPointsZ();
        pointsNull = new byte[nRows][nCols];
        setPointsNull(nullType);
        this.hasNulls = hasNulls;
    }

    public void deleteGrid()
    {
        pointsZ = null;
        pointsNull = null;
        initializeSurface();
    }

    public void setLighting(boolean lighting)
    {
        geometryChanged = true;
        this.lighting = lighting;
        dbFieldChanged("lighting", lighting);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getLighting()
    {
        return lighting;
    }

    public void nullsFilledChanged()
    {
        textureChanged = true;
        colorListChanged = true;
    }

    public void setAngle()
    {
        super.setAngle();
        if (angle == 0.0)
        {
            colXInc = xInc;
            colYInc = 0.0;
            rowXInc = 0.0;
            rowYInc = yInc;
        }
        else
        {
            colXInc = cosXY * xInc;
            colYInc = sinXY * xInc;
            rowXInc = -sinXY * yInc;
            rowYInc = cosXY * yInc;
        }
    }

    public void setIsPixelMode(boolean b)
    {
        if (isPixelMode == b)
        {
            return;
        }
        isPixelMode = b;
        geometryChanged = true;
        if (newSurfaceTexture != null)
        {
            textureChanged = true;
        }
        dbFieldChanged("isPixelMode", isPixelMode);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getIsPixelMode() { return isPixelMode; }

    public void setSmoothNormals(boolean b)
    {
        if (smoothNormals == b)
        {
            return;
        }
        smoothNormals = b;
        geometryChanged = true;
        normals = null;
        if (newSurfaceTexture != null)
        {
            textureChanged = true;
        }
        dbFieldChanged("smoothNormals", smoothNormals);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getSmoothNormals() { return smoothNormals; }

    public void setContourColor(boolean b)
    {
        if (contourColor == b)
        {
            return;
        }
        contourColor = b;
        geometryChanged = true;
        if (newSurfaceTexture != null)
        {
            textureChanged = true;
            textureTiles = null;
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getContourColor() { return contourColor; }
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_WITH_SPECULAR_LIGHTS; }

    public boolean getUseShader()
    {
        return useShader && currentSurfaceTexture != surfaceDefaultTexture;
    }

    public int getShader() { return textureTiles.shader; }

    public void setUseShader(boolean b)
    {
        if (useShader == b) return;
        useShader = b;
        if (newSurfaceTexture != null)
        {
            textureChanged = true;
            textureTiles = null;
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public void setDisplayFill(boolean b)
    {
        if (displayFill == b)
        {
            return;
        }
        displayFill = b;
        dbFieldChanged("displayFill", displayFill);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getDisplayFill() { return displayFill; }

    /*
       public void setDisplayProperty(boolean b)
       {
           if(displayProperty == b) return;
           displayProperty = b;
           currentModel.win3dDisplayAll();
       }
    */
    public void setDisplayGrid(boolean b)
    {
        if (displayGrid == b)
        {
            return;
        }
        displayGrid = b;
        dbFieldChanged("displayGrid", displayGrid);
        currentModel.viewObjectRepaint(this, this);

    }


    //    public boolean getDisplayProperty() { return displayProperty; }
    public boolean getDisplayGrid()
    {
        return displayGrid;
    }


    public boolean toggleSurfacePickingOn()
    {
        boolean displayFill = getDisplayFill();
        boolean displayGrid = getDisplayGrid();
        if (!displayFill && !displayGrid) setDisplayGrid(true);
        return checkIsLoaded();
    }

    public void toggleSurfacePickingOff()
    {
        setIsVisible(false);
    }

    public void setStsColor(StsColor stsColor)
    {
        if (this.stsColor != null && this.stsColor.equals(stsColor)) return;
        this.stsColor = stsColor;
        dbFieldChanged("stsColor", this.stsColor);
        textureChanged = true;
        currentModel.viewObjectRepaint(this, this);
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setPointsZ(float[][] pointsZ)
    {
        this.pointsZ = pointsZ;
    }


    public void setPointsZ()
    {
        if (pointsZ == null) return;
        for (int row = 0; row < nRows; row++)
            Arrays.fill(pointsZ[row], nullValue);
    }

    public void setAdjPointsZ(float[][] adjPointsZ)
    {
        this.adjPointsZ = adjPointsZ;
    }

    public float[][] getPointsZ()
    {
        return pointsZ;
    }

    public float[][] getAdjPointsZ()
    {
        return adjPointsZ;
    }

    public float[][] getCreateAdjPointsZ()
    {
        if (adjPointsZ == null)
            initializeDepthFromTime();
        return adjPointsZ;
    }

    public float[][] getPoints()
    {
        if (isDepth && adjPointsZ != null)
            return adjPointsZ;
        else
            return pointsZ;
    }

    public void setPointsNull(byte[][] pointsNull)
    {
        this.pointsNull = pointsNull;
    }

    public byte[][] getPointsNull()
    {
        return pointsNull;
    }

    public void setNullZValue(float nullZValue)
    {
        this.nullZValue = nullZValue;
    }

    public float getNullZValue()
    {
        return nullZValue;
    }

    public void setZDomainOriginal(byte zDomain)
    {
        zDomainOriginal = zDomain;
        zDomainSupported = zDomain;
    }

    public byte getZDomainOriginal()
    {
        return zDomainOriginal;
    }

    public void setBornDate(String born)
    {
        if (!StsDateFieldBean.validateDateInput(born))
        {
            bornField.setValue(StsDateFieldBean.convertToString(bornDate));
            return;
        }
        super.setBornDate(born);
        return;
    }

    public void setDeathDate(String death)
    {
        if (!StsDateFieldBean.validateDateInput(death))
        {
            deathField.setValue(StsDateFieldBean.convertToString(deathDate));
            return;
        }
        super.setDeathDate(death);
        return;
    }

    public float getXSize()
    {
        return xMax - xMin;
    }

    public float getYSize()
    {
        return yMax - yMin;
    }

    public float getOffset()
    {
        return offset;
    }

    public float[] getRelativeOrigin()
    {
        return new float[]
            {
                xMin, yMin};
    }

    /*
    public float[] getRelativeRange()
    {
     float xMin = 0.0f;
     float xMax = 0.0f;
     if(rowXInc > 0.0f)
     {
      xMin = (float)rowXInc*iMinNonNull;
      xMax = (float)rowXInc*iMaxNonNull;
     }
     else if(rowXInc < 0.0f)
     {
      xMin = (float)rowXInc*iMaxNonNull;
      xMax = (float)rowXInc*iMinNonNull;
     }
     if(colXInc > 0.0f)
     {
      xMin += (float)colXInc*jMinNonNull;
      xMax += (float)colXInc*jMaxNonNull;
     }
     else if(colXInc < 0.0f)
     {
      xMin += (float)colXInc*jMaxNonNull;
      xMax += (float)colXInc*jMinNonNull;
     }

     float yMin = 0.0f;
     float yMax = 0.0f;
     if(rowYInc > 0.0f)
     {
      yMin = (float)rowYInc*iMinNonNull;
      yMax = (float)rowYInc*iMaxNonNull;
     }
     else if(rowYInc < 0.0f)
     {
      yMin = (float)rowYInc*iMaxNonNull;
      yMax = (float)rowYInc*iMinNonNull;
     }
     if(colYInc > 0.0f)
     {
      yMin += (float)colYInc*jMinNonNull;
      yMax += (float)colYInc*jMaxNonNull;
     }
     else if(colYInc < 0.0f)
     {
      yMin += (float)colYInc*jMaxNonNull;
      yMax += (float)colYInc*jMinNonNull;
     }
     return new float[] { relativeXOrigin+xMin, relativeXOrigin+xMax,
           relativeYOrigin+yMin, relativeYOrigin+yMax,
           zMin, zMax };
    }
    */
    public void setHasNulls(boolean hasNulls)
    {
        this.hasNulls = hasNulls;
    }

    public void setDepthMin(float depthMin)
    {
        this.depthMin = depthMin;
        dbFieldChanged("depthMin", depthMin);
    }

    public void setDepthMax(float depthMax)
    {
        this.depthMax = depthMax;
        dbFieldChanged("depthMax", depthMax);
    }

    //TODO Suggest changing name to getMin since getZMin may imply depth
    public float getZMin()
    {
        if (isDepth) return depthMin;
        else return zMin;
    }

    //TODO Suggest changing name to getMax since getZMax may imply depth
    public float getZMax()
    {
        if (isDepth) return depthMax;
        else return zMax;
    }

    public float getTimeMin()
    {
        return zMin;
    }

    public float getTimeMax()
    {
        return zMax;
    }

    public void setOffset(float offset)
    {
        this.offset = offset;
        geometryChanged = true;
        if (newSurfaceTexture != null)
        {
            textureChanged = true;
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getHasNulls()
    {
        return hasNulls;
    }

    // y is assumed to be in rotated coordinate system
    public float getRowF(float y)
    {
        return StsMath.minMax((y - yMin) / yInc, 0.0f, (float) (nRows - 1));
    }

    // x is assumed to be in rotated coordinate system
    public float getColF(float x)
    {
        return StsMath.minMax((x - xMin) / xInc, 0.0f, (float) (nCols - 1));
    }

    public StsGridPoint getGridPoint(int row, int col)
    {
        StsGridPoint gridPoint = new StsGridPoint(this);
        gridPoint.setRowCol(row, col);
        float z = pointsZ[row][col];
        gridPoint.setZ(z);
        gridPoint.setT(z);
        return gridPoint;
    }

    /** Only used in horizonPicker to set seed for a new surface. */
    public void setGridPoint(StsGridPoint gridPoint)
    {
        int row = gridPoint.row;
        int col = gridPoint.col;
        float z = gridPoint.getZorT();
        byte nullType = gridPoint.nullType;
        pointsZ[row][col] = z;
        setPointNull(row, col, nullType);
    }

    public byte getPointNull(int row, int col)
    {
        if (!insideGrid(row, col))
        {
            return SURF_GAP_NOT_FILLED;
        }
        if (pointsNull == null)
        {
            return SURF_PNT;
        }
        return pointsNull[row][col];
    }

    public boolean isPointNotNull(int row, int col)
    {
        if (pointsNull == null)
        {
            return false;
        }
        return pointsNull[row][col] == SURF_PNT;
    }

    public boolean isPointNull(int row, int col)
    {
        if (pointsNull == null) return false;
        return pointsNull[row][col] != SURF_PNT;
    }

    public boolean isPointNullGap(int row, int col)
    {
        byte pointNull = getPointNull(row, col);
        return pointNull == SURF_GAP || pointNull == SURF_GAP_SET || pointNull == SURF_GAP_FILLED;
    }

    public boolean isPointNullThisType(int row, int col, byte nullType)
    {
        if (pointsNull == null)
        {
            return false;
        }
        return pointsNull[row][col] == nullType;
    }

    public void setPrefixString(String s)
    {
        prefixString = s;
    }

    protected boolean checkGrid(StsWin3d win3d)
    {
        initNonNullRange();

        boolean gridOK = true;

        if (nCols <= 0)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,
                "Number columns incorrect: " + nCols +
                    " for grid " + getName());
            gridOK = false;
        }
        if (nRows <= 0)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,
                "Number rows incorrect: " + nRows +
                    " for grid " + getName());
            gridOK = false;
        }
        if (xInc == 0.0f)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,
                "X increment is zero" +
                    " for grid " + getName());
            gridOK = false;
        }
        if (yInc == 0.0f)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,
                "Y increment is zero" +
                    " for grid " + getName());
            gridOK = false;
        }
        return gridOK;
    }

    public StsSurface cloneSurface()
    {
        StsSurface newSurface = new StsSurface();
        newSurface.copySurface(this);
        return newSurface;
    }

    public StsSurface copySurface(StsSurface surface)
    {
        // Not copied: zMin, zMax
        this.nRows = surface.nRows;
        this.nCols = surface.nCols;
        computeRowColStride();
        this.xOrigin = surface.xOrigin;
        this.yOrigin = surface.yOrigin;
        this.xInc = surface.xInc;
        this.yInc = surface.yInc;
        this.xMin = surface.xMin;
        this.yMin = surface.yMin;
        this.xMax = surface.xMax;
        this.yMax = surface.yMax;
        this.angle = surface.angle;
        this.hasNulls = surface.hasNulls;
        this.nullZValue = surface.nullZValue;
        this.pointsZ = new float[nRows][nCols];
        this.pointsNull = new byte[nRows][nCols];
        if (!checkGrid(currentModel.win3d))
        {
            return null;
        }
        return this;
    }

    protected void checkResetOrigin()
    {
        if (xInc < 0.0f && yInc < 0.0f)
        {
            xOrigin += colXInc * (nCols - 1) + rowXInc * (nRows - 1);
            yOrigin += colYInc * (nCols - 1) + rowYInc * (nRows - 1);
            xInc = -xInc;
            yInc = -yInc;
        }
        else if (xInc < 0.0f && yInc > 0.0f)
        {
            xOrigin += colXInc * (nCols - 1);
            yOrigin += colYInc * (nCols - 1);
            xInc = -xInc;
        }
        else if (xInc > 0.0f && yInc < 0.0f)
        {
            xOrigin += rowXInc * (nRows - 1);
            yOrigin += rowYInc * (nRows - 1);
            yInc = -yInc;
        }
    }

    protected boolean initializePointsZ()
    {
        int row, col;

        if (pointsZ == null)
        {
            return true;
        }

        /** Set nullFlag and zMin, zMax */

        pointsNull = new byte[nRows][nCols];
        zMin = largeFloat;
        zMax = -largeFloat;

        for (row = 0; row < nRows; row++)
        {
            for (col = 0; col < nCols; col++)
            {
                float z = pointsZ[row][col];

                if (z == nullZValue)
                {
                    setPointNull(row, col, SURF_BOUNDARY);
                    hasNulls = true;
                }
                else
                {
                    setPointNull(row, col, SURF_PNT);
                    if (z < zMin)
                    {
                        zMin = z;
                    }
                    if (z > zMax)
                    {
                        zMax = z;
                    }
                }
            }
        }

        if (zMin <= 0.0f && zMax <= 0.0f)
        {
            flipZValues();
        }
        if (zMin <= 0.0f || zMax <= 0.0f)
        {
            boolean isSubSea = StsYesNoDialog.questionValue(currentModel.win3d, "Z range is from " + zMin + " to " + zMax + ".\nAre these subsea values, i.e., negative below sea level?");
            if (isSubSea) flipZValues();
        }
        // if time values, range is stored in zMin zMax
        // if depth, stored in depthMin depthMax
        if (zDomainOriginal == StsParameters.TD_DEPTH)
        {
            depthMin = zMin;
            depthMax = zMax;
        }
        return true;
    }

    protected void initializePoints()
    {
        pointsNull = new byte[nRows][nCols];
        pointsZ = new float[nRows][nCols];

        byte currentNullType = 0;
        if (Main.debugPoint)
            currentNullType = pointsNull[Main.debugPointRow][Main.debugPointCol];
        for (int row = 0; row < nRows; row++)
        {
            Arrays.fill(pointsZ[row], nullZValue);
            Arrays.fill(pointsNull[row], SURF_BOUNDARY);
        }
        if (Main.debugPoint)
        {
            StsException.systemDebug(this, "setPointsNull",
                getName() + "pointsNull[" + Main.debugPointRow + "][" + Main.debugPointCol + "] changed to " + StsParameters.getSurfacePointTypeName(StsParameters.SURF_BOUNDARY) +
                    " from " + StsParameters.getSurfacePointTypeName(currentNullType));
        }
    }

    public void computeZRange()
    {
        zMin = largeFloat;
        zMax = -largeFloat;

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                if (pointsNull[row][col] != SURF_PNT)
                {
                    continue;
                }
                float z = pointsZ[row][col];
                if (z == nullZValue)
                {
                    continue;
                }
                if (z < zMin)
                {
                    zMin = z;
                }
                if (z > zMax)
                {
                    zMax = z;
                }
            }
        }
    }

    private void setOffsetRange()
    {
        //        StsFloatFieldBean offsetBean = (StsFloatFieldBean) StsFieldBean.
        //            getBeanWithFieldName(displayFields, "offset");
        StsProject project = currentModel.getProject();
        // float projectZMin = project.getZorTMin();
        // float projectZMax = project.getZorTMax();
        // float dZMin = Math.min(zMin, projectZMin);
        // float dZMax = Math.max(zMax, projectZMax);
        float dz = zMax - zMin;
        if (offsetBean != null)
        {
            offsetBean.setValueAndRangeFixStep(0.0, project.getZorTMin(), project.getZorTMax(), project.getZorTInc());
            offsetBean.setContinuousMode(false);
        }
    }

    /** reorderX and reorderY called if inc is negative;
     *  min and max also need to be flipped, but check in
     *  case its already been done.
     */
    /*
    private void reorderX()
    {
     xInc = -xInc;

//        if(xMin <= xMax) return;

     float temp = xMin;
     xMin = xMax;
     xMax = temp;
    }

    private void reorderY()
    {
     yInc = -yInc;

     if(yMin <= yMax) return;

     float temp = yMin;
     yMin = yMax;
     yMax = temp;
    }
    */
    /*
    private void reverseRows(float[][] v, int nRows, int nCols)
    {
     if(v == null) return;
     for(int row = 0; row < nRows; row++)
     {
      float[] rowV = v[row];
      int first = 0;
      int last = nCols-1;
      while(first < last)
      {
       float temp = rowV[first];
       rowV[first] = rowV[last];
       rowV[last] = temp;
       first++; last--;
      }
     }
    }

    private void reverseCols(float[][] v, int nRows, int nCols)
    {
     if(v == null) return;
     int first = 0;
     int last = nRows-1;
     while(first < last)
     {
      float[] temp = v[first];
      v[first] = v[last];
      v[last] = temp;
      first++; last--;
     }
    }
    */
    /** reverse the z values */
    private void flipZValues()
    {
        float temp = zMin;
        zMin = -zMax;
        zMax = -temp;

        if (hasNulls)
        {
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (pointsNull[row][col] == SURF_PNT)
                    {
                        pointsZ[row][col] *= -1;
                    }
                }
            }
        }
        else
        {
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    pointsZ[row][col] *= -1;
                }
            }
        }
    }

    public float[] getXYZnotNull(int row, int col)
    {
        if (pointsZ == null || !insideGrid(row, col))
        {
            return null;
        }

        int n = getIndex(row, col);
        if (pointsNull != null && pointsNull[row][col] != SURF_PNT)
        {
            return null;
        }
        return getXYZorT(row, col);
    }

    public float[] getXY(int row, int col)
    {
        return new float[]
            {
                (float) xMin + col * xInc, (float) yMin + row * yInc};
    }

    public float[] getXYZorT(int row, int col)
    {
        return new float[]{ xMin + col * xInc, yMin + row * yInc, getZorT(row, col) };
    }

    public final float getZorT(int row, int col)
    {
        if (isDepth && adjPointsZ != null)
            return adjPointsZ[row][col];
        else
            return pointsZ[row][col];
    }

    public final float getZorT(byte zDomain, int row, int col)
    {
        if (zDomain == StsProject.TD_DEPTH && adjPointsZ != null)
            return adjPointsZ[row][col];
        else
            return pointsZ[row][col];
    }


    public float getZNotNull(int row, int col)
    {
        if (pointsZ == null || !insideGrid(row, col))
        {
            return nullZValue;
        }
        if (pointsNull != null && pointsNull[row][col] != SURF_PNT)
        {
            return nullZValue;
        }
        return getZorT(row, col);
    }

    /** Only used in surface construction from markers. */
    public void setGridPoint(int row, int col, float z)
    {
        pointsZ[row][col] = z;
        setPointNull(row, col, SURF_PNT);
    }

    public StsPoint[] getRectanglePoints(int rowMin, int rowMax, int colMin,
                                         int colMax)
    {
        int row, col;

        int nRows = rowMax - rowMin + 1;
        int nCols = colMax - colMin + 1;
        int nPoints = 2 * (nRows + nCols);

        StsPoint[] points = new StsPoint[nPoints];
        int n = 0;
        for (row = rowMin; row <= rowMax; row++)
        {
            points[n++] = getPointZorT(row, colMin);
        }
        for (col = colMin; col <= colMax; col++)
        {
            points[n++] = getPointZorT(rowMax, col);
        }
        for (row = rowMax; row >= rowMin; row--)
        {
            points[n++] = getPointZorT(row, colMax);
        }
        for (col = colMax; col >= colMin; col--)
        {
            points[n++] = getPointZorT(rowMin, col);
        }
        return points;
    }

    /** null out grid cells that are within fault polygons */
    public void applyFaultPolygons(StsFaultPolygon[] faultPolygons)
    {
        if (faultPolygons == null)
        {
            return;
        }
        ArrayList[] rowColTypeList = new ArrayList[faultPolygons.length];
        for (int i = 0; i < faultPolygons.length; i++)
        {
            rowColTypeList[i] = new ArrayList(5);
        }
        for (int col = 0; col < nCols; col++)
        {
            rowLoop:
            for (int row = 0; row < nRows; row++)
            {
                for (int i = 0; i < faultPolygons.length; i++)
                {
                    float[] xy = getXY(row, col);
                    double absoluteX = xy[0] + xOrigin;
                    double absoluteY = xy[1] + yOrigin;
                    if (faultPolygons[i].contains(absoluteX, absoluteY))
                    {
                        rowColTypeList[i].add(new int[]
                            {row, col,
                                (int) pointsNull[row][col]});
                        setPointNull(row, col, SURF_GAP);
                        continue rowLoop;
                    }
                }
            }
        }

        boolean changed = false;
        for (int i = 0; i < faultPolygons.length; i++)
        {
            int nNulled = rowColTypeList[i].size();
            if (nNulled < 1)
            {
                continue;
            }

            int[] rows = new int[nNulled];
            int[] cols = new int[nNulled];
            byte[] types = new byte[nNulled];
            for (int j = 0; j < nNulled; j++)
            {
                int[] rowColType = (int[]) rowColTypeList[i].get(j);
                rows[j] = rowColType[0];
                cols[j] = rowColType[1];
                types[j] = (byte) rowColType[2];
                changed = true;
            }
            rowColTypeList[i].clear();
            faultPolygons[i].setRows(rows);
            faultPolygons[i].setCols(cols);
            faultPolygons[i].setNullTypes(types);
        }

        if (changed)
        {
            constructGrid();
        }
    }

    /** restore nulled out grid cells that are within fault polygons */
    public void removeFaultPolygons(StsFaultPolygon[] faultPolygons)
    {
        if (faultPolygons == null)
        {
            return;
        }
        boolean changed = false;
        for (int i = 0; i < faultPolygons.length; i++)
        {
            int[] rows = faultPolygons[i].getRows();
            if (rows == null)
            {
                continue;
            }
            int[] cols = faultPolygons[i].getCols();
            byte[] types = faultPolygons[i].getNullTypes();
            for (int j = 0; j < rows.length; j++)
            {
                setPointNull(rows[j], cols[j], types[j]);
                changed = true;
            }
            faultPolygons[i] = null;
        }

        if (changed)
        {
            constructGrid();
        }
    }

    protected void setGridPoints(float[] xmkr, float[] ymkr, float[] zmkr, float maxDist, float power)
    {
        if (xmkr == null || ymkr == null || zmkr == null)
        {
            return;
        }
        float incDist = (float) Math.sqrt(xInc * xInc + yInc * yInc);
        for (int col = 0; col < nCols; col++)
        {
            rowLoop:
            for (int row = 0; row < nRows; row++)
            {
                float[] xy = getXY(row, col);
                float x = xy[0];
                float y = xy[1];
                float wtSum = 0.0f;
                float zSum = 0.0f;
                for (int i = 0; i < xmkr.length; i++)
                {
                    float dx = x - xmkr[i];
                    float dy = y - ymkr[i];
                    float dist = (float) Math.sqrt(dx * dx + dy * dy);
                    if (dist < incDist) //== 0.0f)
                    {
                        setGridPoint(row, col, zmkr[i]);
                        continue rowLoop;
                    }
                    if (dist > maxDist)
                    {
                        continue;
                    }
                    float pow = (float) Math.pow(dist, power);
                    if (pow == 0.0f)
                    {
                        continue;
                    }
                    float wt = 1.0f / pow;
                    wtSum += wt;
                    zSum += zmkr[i] * wt;
                }
                if (wtSum != 0.0f)
                {
                    setGridPoint(row, col, zSum / wtSum);
                }
                else
                {
                    setGridPoint(row, col, nullZValue);
                    setPointNull(row, col, SURF_GAP);
                }
                float z = (wtSum == 0.0f) ? nullZValue : zSum / wtSum;
                setGridPoint(row, col, z);
            }
        }
    }

    /* check for z null */
    public boolean notNull(int row, int col)
    {
        if (!insideGrid(row, col))
        {
            return false;
        }
        int n = getIndex(row, col);
        return pointsNull[row][col] == SURF_PNT &&
            pointsZ[row][col] != nullZValue;
    }

    public boolean insideGrid(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

    public boolean insideGrid(float rowF, float colF)
    {
        if (rowF < 0.0f || rowF > nRows - 1)
        {
            return false;
        }
        if (colF < 0.0f || colF > nCols - 1)
        {
            return false;
        }
        return true;
    }

    public boolean setup()
    {
        try
        {
            computeRowColStride();
            textureChanged = true;
            geometryChanged = true;
            if (pointsZ != null)
            {
                interpolateDistanceTransform();
                constructGrid();
            }
            setAngle(); // sets transients rowXInc, etc
            setRelativeRotationAngle();
        }
        catch (Exception e)
        {}

        return true;
    }

    public boolean pointsLoaded(byte surfaceType)
    {
        if (pointsZ != null)
        {
            return true;
        }

        //        if(surfaceType == StsSurface.IMPORTED)
        //            return retrieveValuesFromDataStore();
        //        else
        return false;
    }

    /** we can unload imported surfaces only: too dangerous for model surfaces */
    public void unload()
    {
        setIsVisible(false);
        if (getType() == StsSurface.IMPORTED)
        {
            // these are the only arrays currently used by imported surfaces
            // don't unload for now: needed for modelSurface
            //            pointsZ = null;
            //           pointsNull = null;
            //            normals = null;
        }
    }

    public void constructGrid()
    {
        textureChanged = true;
        //        addTextureToDeleteList(currentModel.win3d.glPanel3d);
        normals = null;
        geometryChanged = true;

        try
        {
            //			if(hasNulls && !fillGaps) computeTStrips();
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.constructGrid() failed.", e,
                StsException.WARNING);
        }
    }

    private void initNonNullRange()
    {
        iMinNonNull = 0;
        iMaxNonNull = nRows - 1;
        jMinNonNull = 0;
        jMaxNonNull = nCols - 1;
    }

    public int getMinNonNullRow()
    {
        return iMinNonNull;
    }

    public int getMaxNonNullRow()
    {
        return iMaxNonNull;
    }

    public int getMinNonNullCol()
    {
        return jMinNonNull;
    }

    public int getMaxNonNullCol()
    {
        return jMaxNonNull;
    }

    /*
    public float getNonNullXMin()
    {
     return xMin + jMinNonNull*xInc;
    }

    public float getNonNullXMax()
    {
     return xMin + jMaxNonNull*xInc;
    }

    public float getNonNullYMin()
    {
     return yMin + iMinNonNull*yInc;
    }

    public float getNonNullYMax()
    {
     return yMin + iMaxNonNull*yInc;
    }
    */
    /* interpolate gap nulls & set boundary nulls */

    public void interpolateDistanceTransform()
    {
        interpolateDistanceTransform(null);
    }

    public void interpolateDistanceTransform(StsProgressPanel progressPanel)
    {
        StsDistanceTransformInterpolation distanceTransform = new StsDistanceTransformInterpolation(this, debugDistTrans);
        distanceTransform.interpolateDistanceTransform(progressPanel);
    }

    public float[][] initializeDistances()
    {
        if (pointsNull == null)
        {
            return null;
        }

        float[][] distances = new float[nRows][nCols];
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                if (pointsNull[row][col] == SURF_PNT)
                {
                    distances[row][col] = 0.0f;
                }
                else
                {
                    distances[row][col] = largeFloat;
                }
            }
        }
        return distances;
    }

    public float[] getDistanceParameters()
    {
        return new float[]
            {
                yInc, xInc};
    }

    public float doDistanceTransform(int row, int col, StsDistanceTransformPoint[] points)
    {
        float weight;

        float nearestDistance = largeFloat;
        boolean debug = debugDistanceTransform(row, col);
        float z = 0;
        for (int n = 1; n < 5; n++)
        {

            float distance = points[n].getDistance();
            int nearestRow = points[n].row;
            int nearestCol = points[n].col;
            if (debug)
            {
                System.out.println("    row: " + nearestRow + " col: " +
                    nearestCol + " distance: " + distance +
                    " z: " + pointsZ[nearestRow][nearestCol]);
            }
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                z = pointsZ[nearestRow][nearestCol];
            }
        }

        float d0 = points[0].distance;

        if (nearestDistance == largeFloat)
        {
            if (debug)
            {
                System.out.println("    no available points to interpolate");
            }
            return d0;
        }

        if (d0 != largeFloat)
        {
            float z0 = pointsZ[row][col];
            float w0 = 1.0f / d0;
            float w1 = 1.0f / nearestDistance;
            float f = w0 / (w0 + w1);
            if (debug)
            {
                System.out.print("  interpolation between old z: " + z0 +
                    " and new z: " + z + " interpolation factor: " +
                    f + " nearest distance: " + nearestDistance);
            }
            z = f * z0 + (1 - f) * z;
            if (debug)
            {
                System.out.println(" resulting z: " + z);
            }
            if (d0 < nearestDistance)
            {
                nearestDistance = d0;
            }
        }
        else
        {
            if (debug)
            {
                System.out.println("    new z value at this point: " + z +
                    " nearest distance: " + nearestDistance);
            }
        }
        //        if(wtSum == 0.0f) return largeFloat;

        //        float z;
        //        if(wtSum == 1.0f)
        //            z = zSum;
        //        else
        //            z = zSum/wtSum;


        pointsZ[row][col] = z;
        setPointNull(row, col, SURF_GAP_FILLED);
        /*
        if (nearestDistance < maxInterpolationDistance)
            setPointNull(row, col, NULL_GAP_FILLED);
        else
            setPointNull(row, col, NULL_BOUNDARY);
        */
        return nearestDistance;
    }

    boolean debugDistanceTransform(int row, int col)
    {
        if (!debugDistTrans)
        {
            return false;
        }
        return row == debugDistTransRow && col == debugDistTransCol;
    }

    /* interpolate gap nulls & set boundary nulls */
    public void interpolateNullZPointsRadiusWeighted()
    {
        StsInterpolationRadiusWeighted interpolator = StsInterpolationRadiusWeighted.constructor(this);
        interpolator.interpolate();
    }

    /* interpolate gap nulls & set boundary nulls */
    /*
    public void interpolateNullZPoints()
    {
     int i, j, n0, n, nn, ii;
     float rdGap;
     float f;
     boolean inGap;
     float z;
     int gapMin, gapMax;
     float zGapMin, zGapMax;
     float newZ, newWt;
     byte pointNull;

     // More than maxGapWidth number of nulls is considered a gap and
     //  flagged as NULL_GAP, otherwise NULL_BOUNDARY.
     //
     int maxGapWidth = (int)(0.1f*Math.max(nRows, nCols));

     if(pointsNull == null) return;

     boolean mainDebug = false; // StsTrace.getTrace();

     weights = new float[nRows][nCols];

     iMinNonNull = nRows;
     iMaxNonNull = -1;
     jMinNonNull = nCols;
     jMaxNonNull = -1;

     int imin, imax, jmin, jmax;
     byte nullType;

     for (i=0, n0=0; i<nRows; i++, n0 = n0+nCols) // iterate over rows
     {
      // Find first non-null in row

      jmin = -1;
      jmax = nCols-1;
      n = n0;
      for(j = 0; j < nCols; j++, n++)
      {
       if (pointsNull[i][j] == NOT_NULL)
       {
        jmin = j;
        if (j<jMinNonNull) jMinNonNull = j;
        if (i<iMinNonNull) iMinNonNull = i;
        break;
       }
      }

      // If a non-null was found, set values to left of first non-null
      //  to same Z value and NULL_BOUNDARY type

      if(jmin >= 0)
      {
       n = n0;
       z = pointsZ[i][jmin];
       for (j=0; j<jmin; j++, n++)   // west
       {
        setPointNull(i, j, NULL_BOUNDARY);
        pointsZ[i][j] = z;
       }

       // If we found a non-null searching from the left, search from
       //  right back to jmin to find last non-null in row . Since we
       //  found a non-null searching from left, we are guaranteed to
       //  find one searching from the right.
       //

       for(j = nCols-1, n = n0+nCols-1; j >= jmin; j--, n--)
       {
        if (pointsNull[i][j] == NOT_NULL)
        {
         jmax = j;
         if (j>jMaxNonNull) jMaxNonNull = j;
         if (i>iMaxNonNull) iMaxNonNull = i;
         break;
        }
       }

       // Set values to right of last non-null to same
       //  Z value and NULL_BOUNDARY type.
       //

       n = n0 + jmax;
       z = pointsZ[i][jmax];
       for (j=jmax+1; j<nCols; j++, n++) // east
       {
        setPointNull(i, j, NULL_BOUNDARY);
        pointsZ[i][j] = z;
       }

       // For interior gaps, linearly interpolate z values from ends
       //  and assign weights to be used in interpolation using columns
       //  in next section
       //

       inGap = false;
       gapMin = -1;
       for(j = jmin, n = n0+jmin; j <= jmax; j++, n++)
       {
        nullType = pointsNull[i][j];
        if(!inGap)
        {
         if(nullType == NOT_NULL)
          gapMin = j;
         else
          inGap = true;
        }
        else // inGap == true
        {
         if(nullType == NOT_NULL && gapMin >= 0)
         {
          gapMax = j;
          rdGap = 1.0f/(gapMax - gapMin);
          zGapMin = pointsZ[i][gapMin];
          zGapMax = pointsZ[i][gapMax];

          if(gapMax - gapMin + 1 <= maxGapWidth)
           pointNull = NULL_GAP;
          else
           pointNull = NULL_GAP_OR_BOUNDARY;

          for(int jj = gapMin+1; jj < gapMax; jj++)
          {
           f = (jj - gapMin)*rdGap;
           pointsZ[i][jj] = zGapMin + f*(zGapMax - zGapMin);
           setPointNull(i, jj, pointNull);
           weights[i][jj] = 1.0f/(jj - gapMin) + 1.0f/(gapMax - jj);
          }

          inGap = false;
          gapMin = j;
         }
        }
       }
      }
      if (mainDebug) System.out.println("null boundary:  points i="
       + i + "  jWest=0-" + (jmin-1) + "  jEast=" + (jmax+1)
       + "-" + (nCols-1));
      }

     // Now search up each column until we find first NOT_NULL or NULL_BOUNDARY
     //  (set in loops above).  Set the values below in column to this z-value
     //  and set nullType to NULL_BOUNDARY
     //
     for (j=0; j<nCols; j++)
     {
      imin = -1;
      imax = nRows-1;

      for (i=0; i<nRows; i++)
      {
       nullType = pointsNull[i][j];

       if(nullType == NOT_NULL)
      //             if(nullType == NOT_NULL || nullType == NULL_BOUNDARY || nullType == NULL_GAP_OR_BOUNDARY)
       {
//                    if(nullType == NULL_GAP_OR_BOUNDARY) pointsNull[n] = NULL_BOUNDARY;
        imin = i;
        break;
       }
      }

      // If a non-null was found, set values below first non-null
      //  to same Z value and NULL_BOUNDARY type

      if(imin >= 0)
      {
       z = pointsZ[imin][j];
       for (i=0; i<imin; i++)
       {
        nullType = checkSetBoundaryOrGapNull(pointsNull[i][j], NULL_BOUNDARY);
        setPointNull(i, j, nullType);
        pointsZ[i][j] = z;
//					if(nullType == NULL_BOUNDARY) pointsZ[i][j] = z;
       }

       // If we found a non-null searching from the bottom, search from
       //  top down to imin to find last non-null in col. Since we
       //  found a non-null searching from bottom, we are guaranteed to
       //  find one searching from the top.
       //

       imax = nRows-1;
       for(i = nRows-1; i >= imin; i--)
       {
        nullType = pointsNull[i][j];

        if(nullType == NOT_NULL)
//                    if(nullType == NOT_NULL || nullType == NULL_BOUNDARY || nullType == NULL_GAP_OR_BOUNDARY)
        {
//                        if(nullType == NULL_GAP_OR_BOUNDARY) pointsNull[n] = NULL_BOUNDARY;
         imax = i;
         break;
        }
       }

       // Set values above last non-null in col to same
       //  Z value and NULL_BOUDARY type.
       //

       z = pointsZ[imax][j];
       for (i=imax+1; i<nRows; i++)
       {
        nullType = checkSetBoundaryOrGapNull(pointsNull[i][j], NULL_BOUNDARY);
        setPointNull(i, j, nullType);
        pointsZ[i][j] = z;
//					if(nullType == NULL_BOUNDARY) pointsZ[i][j] = z;
       }

       // For interior gaps, linearly interpolate z values from ends
       //  and interpolate along with col interpolated values from above
       //

       inGap = false;
       gapMin = -1;
       zGapMin = 0.0f;

       for(i = imin; i <= imax; i++)
       {
        nullType = pointsNull[i][j];

        if(!inGap)
        {
         if(nullType == NOT_NULL)
         {
          gapMin = i;
          zGapMin = pointsZ[i][j];
         }
         else
          inGap = true;
        }
        else // ingap == true
        {
         if(nullType == NOT_NULL && gapMin >= 0)
         {
          gapMax = i;
          zGapMax = pointsZ[i][j];
          rdGap = 1.0f/(gapMax - gapMin);

          if(gapMax - gapMin + 1 <= maxGapWidth)
           pointNull = NULL_GAP;
          else
           pointNull = NULL_GAP_OR_BOUNDARY;

          for(ii = gapMin+1; ii < gapMax; ii++)
          {
           f = (ii - gapMin)*rdGap;
           nn = index(ii, j);

           newZ = zGapMin + f*(zGapMax - zGapMin);
           newWt = 1.0f/(ii - gapMin) + 1.0f/(gapMax - ii);
           pointsZ[ii][j] = (newZ*newWt + pointsZ[ii][j]*weights[ii][j])/(newWt + weights[ii][j]);
           nullType = checkSetBoundaryOrGapNull(pointsNull[ii][j], pointNull);
           setPointNull(ii, j, nullType);
          }

          inGap = false;
          gapMin = i;
         }
        }
       }
      }

      if (mainDebug) System.out.println("null boundary:  points j="
        + j + "  iSouth=0-" + (imin-1) + "  iNorth=" + (imax+1)
        + "-" + (nRows-1));
     } // for j

     // now fill in the corners

           if(iMinNonNull < nRows)
           {
               if(jMinNonNull < nCols)
               {
             // lower left
                   z = pointsZ[iMinNonNull][jMinNonNull];
                   for(i = 0; i < iMinNonNull; i++)
                       for(j = 0; j < jMinNonNull; j++)
                           pointsZ[i][j] = z;
               }
               if(jMaxNonNull >= 0)
               {
                   // lower right
                   z = pointsZ[iMinNonNull][jMaxNonNull];
                   for(i = 0; i < iMinNonNull; i++)
                       for(j = jMaxNonNull+1; j < nCols; j++)
                           pointsZ[i][j] = z;
               }
           }
           if(iMaxNonNull >= 0)
           {
               if(jMinNonNull < nCols)
               {
                   // upper left
                   z = pointsZ[iMaxNonNull][jMinNonNull];
                   for(i = iMaxNonNull+1; i < nRows; i++)
                       for(j = 0; j < jMinNonNull; j++)
                           pointsZ[i][j] = z;
               }
               if(jMaxNonNull < nCols)
               {
                   // upper right
                   z = pointsZ[iMaxNonNull][jMaxNonNull];
                   for(i = iMaxNonNull+1; i < nRows; i++)
                       for(j = jMaxNonNull+1; j < nCols; j++)
                           pointsZ[i][j] = z;
               }
           }
     weights = null; // free scratch array

     // This method interpolates in 8 cardinal directions: good but expensive.
     // Method above interpolates and takes weighted average of 4 cardinal
     // directions: not as good but much faster.

     // interpolateGapsOctantSearch();
    }
    */
    // If one nullType is BOUNDARY and the other is BOUNDARY or GAP_OR_BOUNDARY,
    //  the resulting type is BOUNDARY; otherwise GAP.
    //
    /*
    private byte checkSetBoundaryOrGapNull(byte currentPointNull, byte newPointNull)
    {
     if(newPointNull != NULL_BOUNDARY && newPointNull != NULL_GAP && newPointNull != NULL_GAP_OR_BOUNDARY)
     {
      StsException.outputException(new StsException(StsException.DEBUG,
    "StsGrid.CheckSetboundaryOrGapNULL() called with illegal newPointNull type: " +
      StsParameters.nullName(newPointNull)));
      return currentPointNull;
     }

     if(currentPointNull == NULL_BOUNDARY)
     {
    if(newPointNull == NULL_BOUNDARY || newPointNull == NULL_GAP_OR_BOUNDARY)
       return NULL_BOUNDARY;
      else if(newPointNull == NULL_GAP)
       return NULL_GAP;
      else
       return currentPointNull;
     }
     else if(currentPointNull == NULL_GAP_OR_BOUNDARY)
     {
    if(newPointNull == NULL_BOUNDARY || newPointNull == NULL_GAP_OR_BOUNDARY)
       return NULL_BOUNDARY;
      else if(newPointNull == NULL_GAP)
       return NULL_GAP;
      else
       return currentPointNull;
     }
     else if(currentPointNull == NULL_GAP)
      return NULL_GAP;
     else
      return currentPointNull;
    }

    private void interpolateGapsOctantSearch()
    {
     int i, j, n;
     float z;

     // gap interpolation
     //
     // note:  compare against NULL_GAP to interpolate only null gaps,
     // otherwise compare against NOT_NULL to interpolate all nulls in map bounds;
     // the later "fills in" boundaries, which may not be desirable, but it
     // also guarantees all fault gaps are interpolated
     for (i=iMinNonNull; i<=iMaxNonNull; i++)
     {
      for (j=jMinNonNull; j<=jMaxNonNull; j++)
      {
       if (pointsNull[i][j] != NOT_NULL)
       {
        z = getZFromNeighbors(i, j);
        if (z != nullZValue) pointsZ[i][j] = z;

        if (mainDebug) System.out.println("z(" + i + "," + j + ") = "
          + z + " (null gap interpolated)");
       }
      }
     }
    }

//  interpolate z value in points from a weighted average from
     the nearest neighbors in 8 directions
    private float getZFromNeighbors(int i0, int j0)
    {
     boolean foundOne = false;
     int i, j, n, nn;

     // If this point is not Null, we are done
     if (pointsNull[i0][j0] == NOT_NULL) return pointsZ[i0][j0];

     boolean[] isNeighbor = new boolean[8];  // N first, then clockwise
     int[] dI = new int[] { 1, -1, 0, 0, 1, 1, -1, -1 };
     int[] dJ = new int[] { 0, 0, 1, -1, 1, -1, 1, -1 };

     // Set null values
     for (i=0; i<8; i++)
      isNeighbor[i] = false;

     // north
     for (i=i0+1; i<nRows; i++, i++)
     {
      if(pointsNull[i][j0] == NOT_NULL)
      {
       isNeighbor[0] = true;
       foundOne = true;
       break;
      }
     }

     // south
     for (i=i0-1; i>=0; i--)
     {
      if(pointsNull[i][j0] == NOT_NULL)
      {
       isNeighbor[1] = true;
       foundOne = true;
       break;
      }
     }

     // east
     for (j=j0+1; j<nCols; j++)
     {
      if(pointsNull[i0][j] == NOT_NULL)
      {
       isNeighbor[2] = true;
       foundOne = true;
       break;
      }
     }

     // west
     for (j=j0-1; j>=0; j--)
     {
      if(pointsNull[i0][j] == NOT_NULL)
      {
       isNeighbor[3] = true;
       foundOne = true;
       break;
      }
     }

     // north-east
     for (i=i0+1, j=j0+1; i<nRows && j<nCols; i++, j++)
     {
      if(pointsNull[i][j] == NOT_NULL)
      {
       isNeighbor[4] = true;
       foundOne = true;
       break;
      }
     }
     // north-west
     for (i=i0+1, j=j0-1; i<nRows && j>=0; i++, j--)
     {
      if(pointsNull[i][j] == NOT_NULL)
      {
       isNeighbor[5] = true;
       foundOne = true;
       break;
      }
     }
     // south-east
     for (i=i0-1, j=j0+1; i>=0 && j<nCols; i--, j++)
     {
      if(pointsNull[i][j] == NOT_NULL)
      {
       isNeighbor[6] = true;
       foundOne = true;
       break;
      }
     }
     // south-west
     for (i=i0-1, j=j0-1; i>=0 && j>=0; i--, j--)
     {
      if(pointsNull[i][j] == NOT_NULL)
      {
       isNeighbor[7] = true;
       foundOne = true;
       break;
      }
     }

     if(!foundOne) return nullZValue;

     // compute weighted average based on distance away

     float zTotal = 0.0f;
     float wtTotal = 0.0f;
     for(n = 0; n < 8; n++)
     {
      boolean neighbor = isNeighbor[n];
      if(!neighbor) continue;

      int di = dI[n];
      int dj = dJ[n];
      float weight = 1.0f/(float)Math.sqrt(di*di + dj*dj);

      zTotal += pointsZ[i0][j0] * weight;
      wtTotal += weight;
     }
     return zTotal/wtTotal;
    }
    */
    //    public void setUnits(byte[] units)
    //    {
    //        this.units = units;
    //        return;
    //    }

    public void setNativeUnits(byte vertUnits, byte horzUnits)
    {
        nativeVerticalUnits = vertUnits;
        nativeHorizontalUnits = horzUnits;
    }

    public void toProjectUnits()
    {
        float vScalar = 1.0f;
        float hScalar = 1.0f;
        if (getZDomainOriginal() == StsParameters.TD_DEPTH)
        {
            vScalar = currentModel.getProject().getDepthScalar(nativeVerticalUnits);
        }
        else
        {
            vScalar = currentModel.getProject().getTimeScalar(nativeVerticalUnits);
        }

        hScalar = currentModel.getProject().getXyScalar(nativeHorizontalUnits);

        if ((vScalar != 1.0f) || (hScalar != 1.0f))
        {
            StsMessageFiles.infoMessage("Native surface units not equal to project units, translating surface");
            xOrigin = xOrigin * hScalar;
            yOrigin = yOrigin * hScalar;
            xInc = xInc * hScalar;
            yInc = yInc * hScalar;
            xMin = xMin * hScalar;
            yMin = yMin * hScalar;
            xMax = xMax * hScalar;
            yMax = yMax * hScalar;
            zMin = zMin * vScalar;
            zMax = zMax * vScalar;

            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (pointsZ[row][col] != nullZValue)
                    {
                        pointsZ[row][col] = pointsZ[row][col] * vScalar;
                    }
                }
            }
        }
        else
        {
            StsMessageFiles.infoMessage("Native surface units equal to project units.");
        }

    }

    /*
       public boolean writeBinaryFile(String dirname, String filename, byte[] units)
       {
           this.units = units;
           return writeBinaryFile(dirname, filename);
       }
    */
    public boolean writeBinaryFile(String dirname, String filename)
    {
        StsBinaryFile binaryFile = null;

        try
        {
            if (pointsZ == null)
            {
                return false;
            }

            binaryFile = getBinaryFile(dirname, filename);
            binaryFile.openWrite(false);
            writeBinaryHeader(binaryFile);

            if (!hasNulls)
            {
                for (int row = 0; row < nRows; row++)
                {
                    binaryFile.setFloatValues(pointsZ[row]);
                }
            }
            else
            {
                for (int row = 0; row < nRows; row++)
                {
                    binaryFile.setFloatValues(pointsZ[row], pointsNull[row],
                        SURF_PNT, nullZValue);
                }
            }
            binaryFile.close();

            /** Output the surface attributes associated with this surface */
            for (int i = 0; i < surfaceAttributes.getSize(); i++)
            {
                StsSurfaceAttribute surfaceAtt = (StsSurfaceAttribute) surfaceAttributes.getElement(i);
                String attFilename = new String(filename + "." + surfaceAtt.getName() + ".0");
                surfaceAtt.writeBinaryFile(attFilename);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurface.writeBinaryFile() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean writeBinaryHeader(StsBinaryFile binaryFile)
    {
        try
        {
            binaryFile.setByteValues(new byte[]{nativeVerticalUnits, nativeHorizontalUnits});
            binaryFile.setIntegerValues(new int[]{nRows, nCols});
            binaryFile.setDoubleValues(new double[]{getXOrigin(), getYOrigin()});
            binaryFile.setFloatValues(new float[]{xInc, yInc, angle, zMin, zMax, nullZValue,
                rowNumMin, rowNumMax, rowNumInc, colNumMin, colNumMax, colNumInc});
            binaryFile.setBooleanValues(new boolean[]{hasNulls});
            binaryFile.setByteValues(new byte[]{zDomainOriginal});
            return true;
        }
        catch (Exception ex)
        {
            StsException.outputException("StsSurface.writeBinaryFile() failed writing header for surface (" + getName() + ")", ex, StsException.WARNING);
            return false;
        }
    }

    private StsBinaryFile getBinaryFile(String dirname, String filename)
    {
        try
        {
            StsFile file = StsFile.constructor(dirname, filename);
            if (file == null)
            {
                return null;
            }
            //            String urlName = new String("file:" + binaryFileDir + File.separator + binaryFilename);
            //            URL url = new URL(urlName);
            //            return new StsBinaryFile(url);
            return new StsBinaryFile(file);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.getBinaryFile() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private boolean getGridNull(byte nullType, boolean mapFillOn)
    {
        if (mapFillOn)
        {
            return (nullType != SURF_BOUNDARY);
        }
        else
        {
            return (nullType == SURF_PNT);
        }
    }

    public boolean changeNullTypes(byte currentType, byte newType)
    {
        boolean changed = false;
        if (pointsNull != null)
        {
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (pointsNull[row][col] == currentType)
                    {
                        setPointNull(row, col, newType);
                        changed = true;
                    }
                }
            }
        }
        return changed;
    }

    public boolean changeNullTypeIfNot(int row, int col, byte excludeType,
                                       byte newType)
    {
        if (pointsNull[row][col] == excludeType)
        {
            return false;
        }
        setPointNull(row, col, newType);
        return true;
    }
/*
    public void fillRemainingGapPoints()
    {
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                if (pointsNull[row][col] == SURF_GAP)
                {
                    setPointNull(row, col, SURF_GAP_FILLED);
                }
            }
        }
    }
*/
    /** Draw any map edges on section */
    private void drawOnCursor3d(StsGLPanel3d glPanel3d, StsPoint point0,
                                StsPoint point1, boolean drawGaps,
                                boolean drawDotted)
    {
        float ftemp;
        int itemp;
        float dColF, dRowF;
        float slope;
        int col, row;
        float colF, rowF;
        int i, j;
        final int lineStart = 0;
        final int lineDraw = 1;
        final int lineNull = 2;
        int lineStatus = lineStart;
        StsColor drawColor;
        boolean fillOn, gridOn, edgeOnGrid;

        if (!getIsVisible())
        {
            return;
        }
        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
        {
            return;
        }

        StsGridPoint g0 = new StsGridPoint(point0, this);
        StsGridPoint g1 = new StsGridPoint(point1, this);
        if (!StsGridDefinition.clipLineToGrid(g0, g1, this, false))
        {
            return;
        }

        StsSurfaceClass surfaceClass = (StsSurfaceClass) currentModel.
            getCreateStsClass(StsSurface.class);

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        int nColsm1 = nCols - 1;
        int nRowsm1 = nRows - 1;

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPushMatrix();
            gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
        }

        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(StsGraphicParameters.edgeLineWidthHighlighted);
        //        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);

        if (drawDotted)
        {
            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
        }

        StsColor nullColor = StsColor.WHITE;

        //	set_line_smoothing(current_window, TRUE);

        float[] point;
        float[] lastPoint = null;

        if (g0.col == g1.col) /* map edge is along a col	*/
        {
            if (g0.rowF > g1.rowF)
            {
                ftemp = g0.rowF;
                g0.rowF = g1.rowF;
                g1.rowF = ftemp;

                itemp = g0.row;
                g0.row = g1.row;
                g1.row = itemp;
            }

            g0.row = (int) g0.rowF;
            g1.row = (int) g1.rowF;

            if (g0.col + 1 <= nColsm1)
            {
                dColF = g0.colF - g0.col;
            }
            else
            {
                g0.col = nColsm1 - 1;
                dColF = 1.0f;
            }

            boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
            gridOn = displaySurfaces && displayGrid;
            fillOn = displaySurfaces && displayFill;
            edgeOnGrid = (dColF <= 0.0f || dColF >= 1.0f);

            /** if cursor off of grid and grid is on or fill is on, draw in black */
            if (drawDotted)
            {
                if ((gridOn && edgeOnGrid) ^ fillOn)
                {
                    drawColor = StsColor.BLACK;
                }
                else
                {
                    drawColor = getStsColor();
                }
            }
            else
            {
                drawColor = StsColor.BLACK;
            }

            drawColor.setGLColor(gl);

            /** If cursor edge on grid, view shift more */
            double viewShift;
            if (gridOn && edgeOnGrid)
            {
                viewShift = StsGraphicParameters.cursorGridShift;
            }
            else
            {
                viewShift = StsGraphicParameters.gridShift;
            }

            if (drawDotted)
            {
                viewShift += 1.0f;
            }
            glPanel3d.setViewShift(gl, viewShift);

            for (i = g0.row; i < g1.row; i++)
            {
                point = interpolateGridPoint(i, g0.col, COL, dColF);
                if (point != null)
                {
                    if (lineStatus == lineNull) // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if (drawGaps)
                        {
                            if (drawDotted)
                            {
                                nullColor.setGLColor(gl);
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint, 0);
                                gl.glVertex3fv(point, 0);
                                gl.glEnd();
                                drawColor.setGLColor(gl);
                            }
                            else
                            {
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint, 0);
                                gl.glVertex3fv(point, 0);
                                gl.glEnd();
                            }
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex3fv(point, 0); // draw regular line start or segment
                    lineStatus = lineDraw;
                    lastPoint = point; // save the point
                }
                else if (lineStatus == lineDraw) // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw)
            {
                gl.glEnd(); // wrap up the last regular segment
            }
        }
        else if (g0.row == g1.row) /* map edge is along a row	*/
        {
            if (g0.col > g1.col)
            {
                ftemp = g0.colF;
                g0.colF = g1.colF;
                g1.colF = ftemp;

                itemp = g0.col;
                g0.col = g1.col;
                g1.col = itemp;
            }

            if (g0.row + 1 <= nRowsm1)
            {
                dRowF = g0.rowF - g0.row;
            }
            else
            {
                g0.row = nRowsm1 - 1;
                dRowF = 1.0f;
            }

            boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
            gridOn = displaySurfaces && displayGrid;
            fillOn = displaySurfaces && displayFill;
            edgeOnGrid = (dRowF <= 0.0f || dRowF >= 1.0f);

            /** if cursor off of grid and grid is on or fill is on, draw in black */
            if (drawDotted)
            {
                if ((gridOn && edgeOnGrid) ^ fillOn)
                {
                    drawColor = StsColor.BLACK;
                }
                else
                {
                    drawColor = getStsColor();
                }
            }
            else
            {
                drawColor = StsColor.BLACK;
            }

            drawColor.setGLColor(gl);

            /** If cursor edge on grid, view shift more */
            double viewShift;
            if (gridOn && edgeOnGrid)
            {
                viewShift = StsGraphicParameters.cursorGridShift;
            }
            else
            {
                viewShift = StsGraphicParameters.gridShift;
            }

            if (drawDotted)
            {
                viewShift += 1.0f;
            }
            glPanel3d.setViewShift(gl, viewShift);
            //			gl.glBegin(GL.GL_LINE_STRIP);

            for (j = g0.col; j < g1.col; j++)
            {
                point = interpolateGridPoint(g0.row, j, ROW, dRowF);
                if (point != null)
                {
                    if (lineStatus == lineNull) // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if (drawGaps)
                        {
                            if (drawDotted)
                            {
                                nullColor.setGLColor(gl);
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint, 0);
                                gl.glVertex3fv(point, 0);
                                gl.glEnd();
                                drawColor.setGLColor(gl);
                            }
                            else
                            {
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint, 0);
                                gl.glVertex3fv(point, 0);
                                gl.glEnd();
                            }
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex3fv(point, 0); // draw regular line start or segment
                    lineStatus = lineDraw;
                    lastPoint = point; // save the point
                }
                else if (lineStatus == lineDraw) // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw)
            {
                gl.glEnd(); // wrap up the last regular segment
            }
        }
        else
        {
            /** if cursor off of grid and grid is on or fill is on, draw in black */
            if (drawDotted)
            {
                fillOn = surfaceClass.getDisplaySurfaces() && displayFill;
                if (fillOn)
                {
                    drawColor = StsColor.BLACK;
                }
                else
                {
                    drawColor = getStsColor();
                }
            }
            else
            {
                drawColor = StsColor.BLACK;
            }

            drawColor.setGLColor(gl);

            int nc = Math.abs(g1.col - g0.col); // number of columns spanned
            int nr = Math.abs(g1.row - g0.row); // number of rows spanned

            if (nc > nr)
            {
                if (g0.col > g1.col)
                {
                    ftemp = g0.colF;
                    g0.colF = g1.colF;
                    g1.colF = ftemp;

                    itemp = g0.col;
                    g0.col = g1.col;
                    g1.col = itemp;

                    ftemp = g0.rowF;
                    g0.rowF = g1.rowF;
                    g1.rowF = ftemp;

                    itemp = g0.row;
                    g0.row = g1.row;
                    g1.row = itemp;
                }

                dColF = g0.colF - g0.col;
                dRowF = g0.rowF - g0.row;

                slope = (g1.rowF - g0.rowF) / (g1.colF - g0.colF);

                col = g0.col + 1;

                rowF = g0.rowF + slope * (g0.col + 1 - g0.colF);
                row = (int) rowF;

                dRowF = rowF - row;

                gl.glBegin(GL.GL_LINE_STRIP);
                for (i = col; i < g1.col; i++)
                {
                    point = interpolateGridPoint(row, i, ROW, dRowF);
                    if (point != null)
                    {
                        if (lineStatus == lineNull) // draw the fault gap
                        {
                            lineStatus = lineStart;
                            if (drawGaps)
                            {
                                if (drawDotted)
                                {
                                    nullColor.setGLColor(gl);
                                    gl.glBegin(GL.GL_LINE_STRIP);
                                    gl.glVertex3fv(lastPoint, 0);
                                    gl.glVertex3fv(point, 0);
                                    gl.glEnd();
                                    drawColor.setGLColor(gl);
                                }
                                else
                                {
                                    gl.glBegin(GL.GL_LINE_STRIP);
                                    gl.glVertex3fv(lastPoint, 0);
                                    gl.glVertex3fv(point, 0);
                                    gl.glEnd();
                                }
                            }
                        }
                        if (lineStatus == lineStart) // start of regular line
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                        }

                        gl.glVertex3fv(point, 0); // draw regular line start or segment
                        lineStatus = lineDraw;
                        lastPoint = point; // save the point
                    }
                    else if (lineStatus == lineDraw) // start of fault gap
                    {
                        gl.glEnd();
                        lineStatus = lineNull;
                    }
                    rowF += slope;
                    row = (int) rowF;
                    dRowF = rowF - row;
                }
                if (lineStatus == lineDraw)
                {
                    gl.glEnd(); // wrap up the last regular segment
                }
            }
            else
            {
                if (g0.row > g1.row)
                {
                    ftemp = g0.colF;
                    g0.colF = g1.colF;
                    g1.colF = ftemp;

                    itemp = g0.col;
                    g0.col = g1.col;
                    g1.col = itemp;

                    ftemp = g0.rowF;
                    g0.rowF = g1.rowF;
                    g1.rowF = ftemp;

                    itemp = g0.row;
                    g0.row = g1.row;
                    g1.row = itemp;
                }

                dColF = g0.colF - g0.col;
                dRowF = g0.rowF - g0.row;

                slope = (g1.colF - g0.colF) / (g1.rowF - g0.rowF);

                row = g0.row + 1;

                colF = g0.colF + slope * (g0.row + 1 - g0.rowF);
                col = (int) colF;

                dColF = colF - col;

                gl.glBegin(GL.GL_LINE_STRIP);
                for (i = row; i < g1.row; i++)
                {
                    point = interpolateGridPoint(i, col, COL, dColF);
                    if (point != null)
                    {
                        if (lineStatus == lineNull) // draw the fault gap
                        {
                            lineStatus = lineStart;
                            if (drawGaps)
                            {
                                if (drawDotted)
                                {
                                    nullColor.setGLColor(gl);
                                    gl.glBegin(GL.GL_LINE_STRIP);
                                    gl.glVertex3fv(lastPoint, 0);
                                    gl.glVertex3fv(point, 0);
                                    gl.glEnd();
                                    drawColor.setGLColor(gl);
                                }
                                else
                                {
                                    gl.glBegin(GL.GL_LINE_STRIP);
                                    gl.glVertex3fv(lastPoint, 0);
                                    gl.glVertex3fv(point, 0);
                                    gl.glEnd();
                                }
                            }
                        }
                        if (lineStatus == lineStart) // start of regular line
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                        }

                        gl.glVertex3fv(point, 0); // draw regular line start or segment
                        lineStatus = lineDraw;
                        lastPoint = point; // save the point
                    }
                    else if (lineStatus == lineDraw) // start of fault gap
                    {
                        gl.glEnd();
                        lineStatus = lineNull;
                    }
                    colF += slope;
                    col = (int) colF;
                    dColF = colF - col;
                }
                if (lineStatus == lineDraw)
                {
                    gl.glEnd(); // wrap up the last regular segment
                }
            }
        }

        /** Reset view shift and turn dotted-lines off */
        glPanel3d.resetViewShift(gl);
        if (drawDotted)
        {
            gl.glDisable(GL.GL_LINE_STIPPLE);
        }
        gl.glEnable(GL.GL_LIGHTING);

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPopMatrix();
        }

        //	set_line_smoothing(current_window, FALSE);
    }

    /** Draw any map edges on section */
    private void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean drawGaps)
    {
        final int lineStart = 0;
        final int lineDraw = 1;
        final int lineNull = 2;
        int lineStatus = lineStart;
        StsColor drawColor;
        boolean fillOn, gridOn;

        if (!getIsVisible() || pointsZ == null) return;
        StsSurfaceClass surfaceClass = (StsSurfaceClass) currentModel.getCreateStsClass(StsSurface.class);

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        int nColsm1 = nCols - 1;
        int nRowsm1 = nRows - 1;

        gl.glDisable(GL.GL_LIGHTING);
        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);

        StsColor nullColor = StsColor.WHITE;

        //	set_line_smoothing(current_window, TRUE);

        float[] point;
        float[] lastPoint = null;

        if (dirNo == StsCursor3d.XDIR) /* map edge is along a col	*/
        {
            float colF = getColCoor(dirCoordinate);
            int col = Math.round(colF);

            boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
            gridOn = displaySurfaces && displayGrid;
            fillOn = displaySurfaces && displayFill;

            //            glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);

            /** if cursor off of grid and grid is on or fill is on, draw in black */
            if (gridOn ^ fillOn)
            {
                drawColor = StsColor.BLACK;
            }
            else
            {
                drawColor = getStsColor();
            }

            drawColor.setGLColor(gl);

            for (int i = 0; i < nRows; i++)
            {
                point = getXYZnotNull(i, col);
                if (point != null)
                {
                    if (lineStatus == lineNull) // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if (drawGaps)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex2f(lastPoint[1], lastPoint[2]);
                            gl.glVertex2f(point[1], point[2]);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex2f(point[1], point[2]); // draw regular line start or segment
                    lineStatus = lineDraw;
                    lastPoint = point; // save the point
                }
                else if (lineStatus == lineDraw) // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw)
            {
                gl.glEnd(); // wrap up the last regular segment
            }
        }
        else if (dirNo == StsCursor3d.YDIR) /* map edge is along a row	*/
        {
            float rowF = getRowCoor(dirCoordinate);
            int row = Math.round(rowF);

            boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
            gridOn = displaySurfaces && displayGrid;
            fillOn = displaySurfaces && displayFill;

            //            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);

            /** if cursor off of grid and grid is on or fill is on, draw in black */
            if (gridOn ^ fillOn)
            {
                drawColor = StsColor.BLACK;
            }
            else
            {
                drawColor = getStsColor();
            }

            drawColor.setGLColor(gl);

            gl.glBegin(GL.GL_LINE_STRIP);
            for (int j = 0; j < nCols; j++)
            {
                point = getXYZnotNull(row, j);
                if (point != null)
                {
                    if (lineStatus == lineNull) // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if (drawGaps)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex2f(lastPoint[0], lastPoint[2]);
                            gl.glVertex2f(point[0], point[2]);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex2f(point[0], point[2]); // draw regular line start or segment
                    lineStatus = lineDraw;
                    lastPoint = point; // save the point
                }
                else if (lineStatus == lineDraw) // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw)
            {
                gl.glEnd(); // wrap up the last regular segment
            }
        }

        /** Reset view shift and turn dotted-lines off */
        //        glPanel3d.resetViewShift(gl);
        gl.glDisable(GL.GL_LINE_STIPPLE);
        gl.glEnable(GL.GL_LIGHTING);
    }

    /** Draw any map edges on section */
    private void drawOn3dCurtain(StsGLPanel3d glPanel3d,
                                 StsSeismicCurtain seismicCurtain,
                                 boolean drawGaps, boolean drawDotted)
    {
        final int lineStart = 0;
        final int lineDraw = 1;
        final int lineNull = 2;
        int lineStatus = lineStart;
        StsColor drawColor;
        boolean fillOn, gridOn;

        if (!getIsVisible())
        {
            return;
        }

        StsSurfaceClass surfaceClass = (StsSurfaceClass) currentModel.
            getCreateStsClass(StsSurface.class);

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        int nColsm1 = nCols - 1;
        int nRowsm1 = nRows - 1;

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPushMatrix();
            gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
        }

        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);

        if (drawDotted)
        {
            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
        }

        StsColor nullColor = StsColor.WHITE;

        boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
        gridOn = displaySurfaces && displayGrid;
        fillOn = displaySurfaces && displayFill;

        /** if cursor off of grid and grid is on or fill is on, draw in black */
        if (drawDotted)
        {
            if (gridOn && fillOn)
            {
                drawColor = StsColor.BLACK;
            }
            else
            {
                drawColor = getStsColor();
            }
        }
        else
        {
            drawColor = StsColor.BLACK;
        }

        drawColor.setGLColor(gl);

        /** If cursor edge on grid, view shift more */
        double viewShift = StsGraphicParameters.gridShift;

        if (drawDotted)
        {
            viewShift += 1.0f;
        }
        glPanel3d.setViewShift(gl, viewShift);

        StsGridPoint[] gridPoints = seismicCurtain.getCellGridPoints();
        float[] lastPoint = null;
        float[] point = gridPoints[0].getPoint().v;
        int nPoints = gridPoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            interpolateBilinear(gridPoints[n]);
            if (gridPoints[n].nullType == SURF_PNT)
            {
                lastPoint = point;
                point = gridPoints[n].getXYZorT(isDepth);
                if (lineStatus == lineNull) // draw the fault gap
                {
                    lineStatus = lineStart;
                    if (drawGaps)
                    {
                        if (drawDotted)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint, 0);
                            gl.glVertex3fv(point, 0);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                        else
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint, 0);
                            gl.glVertex3fv(point, 0);
                            gl.glEnd();
                        }
                    }
                }
                if (lineStatus == lineStart) // start of regular line
                {
                    gl.glBegin(GL.GL_LINE_STRIP);
                }

                gl.glVertex3fv(point, 0); // draw regular line start or segment
                lineStatus = lineDraw;
                lastPoint = point; // save the point
            }
            else if (lineStatus == lineDraw) // start of fault gap
            {
                gl.glEnd();
                lineStatus = lineNull;
            }
        }
        if (lineStatus == lineDraw)
        {
            gl.glEnd(); // wrap up the last regular segment
        }

        /** Reset view shift and turn dotted-lines off */
        glPanel3d.resetViewShift(gl);
        if (drawDotted)
        {
            gl.glDisable(GL.GL_LINE_STIPPLE);
        }
        gl.glEnable(GL.GL_LIGHTING);

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPopMatrix();
        }

        //	set_line_smoothing(current_window, FALSE);
    }
    private void drawOn2dCurtain(StsGLPanel3d glPanel3d,
                                 StsSeismicCurtain seismicCurtain,
                                 boolean drawGaps, boolean drawDotted)
    {
        final int lineStart = 0;
        final int lineDraw = 1;
        final int lineNull = 2;
        int lineStatus = lineStart;
        StsColor drawColor;
        boolean fillOn, gridOn;

        if (!getIsVisible())
        {
            return;
        }

        StsSurfaceClass surfaceClass = (StsSurfaceClass) currentModel.
            getCreateStsClass(StsSurface.class);

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        int nColsm1 = nCols - 1;
        int nRowsm1 = nRows - 1;

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPushMatrix();
            gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
        }

        gl.glDisable(GL.GL_LIGHTING);

        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);

        if (drawDotted)
        {
            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
        }

        StsColor nullColor = StsColor.WHITE;

        boolean displaySurfaces = surfaceClass.getDisplaySurfaces();
        gridOn = displaySurfaces && displayGrid;
        fillOn = displaySurfaces && displayFill;

        /** if cursor off of grid and grid is on or fill is on, draw in black */
        if (drawDotted)
        {
            if (gridOn && fillOn)
            {
                drawColor = StsColor.BLACK;
            }
            else
            {
                drawColor = getStsColor();
            }
        }
        else
        {
            drawColor = StsColor.BLACK;
        }

        drawColor.setGLColor(gl);

        /** If cursor edge on grid, view shift more */
        double viewShift = StsGraphicParameters.gridShift;

        if (drawDotted)
        {
            viewShift += 1.0f;
        }
        glPanel3d.setViewShift(gl, viewShift);

        StsGridPoint[] gridPoints = seismicCurtain.getCellGridPoints();
        float[] lastPoint = null;
        float[] point = gridPoints[0].getPoint().v;
        int nPoints = gridPoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            interpolateBilinear(gridPoints[n]);
            if (gridPoints[n].nullType == SURF_PNT)
            {
                lastPoint = point;
                point = gridPoints[n].getXYZorT(isDepth);
                if (lineStatus == lineNull) // draw the fault gap
                {
                    lineStatus = lineStart;
                    if (drawGaps)
                    {
                        if (drawDotted)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint, 0);
                            gl.glVertex3fv(point, 0);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                        else
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint, 0);
                            gl.glVertex3fv(point, 0);
                            gl.glEnd();
                        }
                    }
                }
                if (lineStatus == lineStart) // start of regular line
                {
                    gl.glBegin(GL.GL_LINE_STRIP);
                }

                gl.glVertex3fv(point, 0); // draw regular line start or segment
                lineStatus = lineDraw;
                lastPoint = point; // save the point
            }
            else if (lineStatus == lineDraw) // start of fault gap
            {
                gl.glEnd();
                lineStatus = lineNull;
            }
        }
        if (lineStatus == lineDraw)
        {
            gl.glEnd(); // wrap up the last regular segment
        }

        /** Reset view shift and turn dotted-lines off */
        glPanel3d.resetViewShift(gl);
        if (drawDotted)
        {
            gl.glDisable(GL.GL_LINE_STIPPLE);
        }
        gl.glEnable(GL.GL_LIGHTING);

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPopMatrix();
        }

        //	set_line_smoothing(current_window, FALSE);
    }

    private float[] interpolateGridPoint(int row, int col, int rowOrCol,
                                         float f)
    {
        float[] point0, point1;

        point0 = getXYZnotNull(row, col);

        if (rowOrCol == ROW)
        {
            point1 = getXYZnotNull(row + 1, col);
        }
        else
        {
            point1 = getXYZnotNull(row, col + 1);
        }

        if (f <= 0.0f)
        {
            if (point0 == null)
            {
                return null;
            }
            return point0;
        }
        else if (f >= 1.0f)
        {
            if (point1 == null)
            {
                return null;
            }
            return point1;
        }
        else
        {
            if (point0 == null)
            {
                return null;
            }
            if (point1 == null)
            {
                return null;
            }

            float[] point = new float[3];

            if (rowOrCol == ROW)
            {
                point[0] = point0[0];
                point[1] = point0[1] + f * (point1[1] - point0[1]);
            }
            else
            {
                point[1] = point0[1];
                point[0] = point0[0] + f * (point1[0] - point0[0]);
            }

            point[2] = point0[2] + f * (point1[2] - point0[2]);

            return point;
        }
    }

    private void getBoundedGridFloatCoordinates(StsGridPoint gridPoint)
    {
        getUnboundedGridFloatCoordinates(gridPoint);
        gridPoint.rowF = StsMath.minMax(gridPoint.rowF, 0.0f, (float) nRows - 1);
        gridPoint.colF = StsMath.minMax(gridPoint.colF, 0.0f, (float) nCols - 1);
    }

    /**
     * Given x and y in gridPoint, set float and int values of row and col.
     * Because int row and col are the indices of the cell containing the point,
     * if point is off the grid, set the float values to the boundary and the
     * int values to the cell indices (i.e., nRows-2 for nRows-1 gridPoints).
     */

    public void defineRowCol(StsGridPoint gridPoint)
    {
        getUnboundedGridFloatCoordinates(gridPoint);

        if (gridPoint.rowF < 0.0f)
        {
            gridPoint.rowF = 0.0f;
            gridPoint.row = 0;
        }
        else if (gridPoint.rowF >= (float) (nRows - 1))
        {
            gridPoint.rowF = (float) (nRows - 1);
            gridPoint.row = nRows - 2;
        }
        else
        {
            gridPoint.row = (int) gridPoint.rowF;
        }

        if (gridPoint.colF < 0.0f)
        {
            gridPoint.colF = 0.0f;
            gridPoint.col = 0;
        }
        else if (gridPoint.colF >= (float) (nCols - 1))
        {
            gridPoint.colF = (float) (nCols - 1);
            gridPoint.col = nCols - 2;
        }
        else
        {
            gridPoint.col = (int) gridPoint.colF;
        }
    }

    /** Compute float values of row and col.  Good for rotated grids as well */
    private void getUnboundedGridFloatCoordinates(StsGridPoint gridPoint)
    {
        float x, y;
        float dx, dy;

        x = gridPoint.getPoint().getX();
        y = gridPoint.getPoint().getY();
        gridPoint.colF = (x - xMin) / xInc;
        gridPoint.rowF = (y - yMin) / yInc;
    }

    public float getRowCoor(float[] xy)
    {
        return (xy[1] - yMin) / yInc;
    }

    public float getRowCoor(float y)
    {
        return (y - yMin) / yInc;
    }

    public float getColCoor(float[] xy)
    {
        return (xy[0] - xMin) / xInc;
    }

    public float getColCoor(float x)
    {
        return (x - xMin) / xInc;
    }

    public float getYCoor(int row, int col)
    {
        return getYCoor((float) row, (float) col);

    }

    public float getYCoor(float rowF, float colF)
    {
        //       rowF = StsMath.minMax(rowF, 0.0f, (float) (nRows - 1));
        return (yMin + rowF * yInc);
    }

    public float getXCoor(int row, int col)
    {
        return getXCoor((float) row, (float) col);
    }

    public float getXCoor(float rowF, float colF)
    {
        //        colF = StsMath.minMax(colF, 0.0f, (float) (nCols - 1));
        return (xMin + colF * xInc);
    }

    public float[] getRowColCoors(StsPoint xy)
    {
        float[] rowCol = new float[2];
        rowCol[0] = getRowCoor(xy.v); // y to row coordinates
        rowCol[1] = getColCoor(xy.v); // x to col coordinates
        return rowCol;
    }

    public int[] getNearestRowColCoors(float x, float y)
    {
        int row = getNearestRowCoor(y);
        int col = getNearestColCoor(x);
        return new int[]{row, col};
    }

    /*
        public float getRowColZ(int row, int col)
        {
            if (pointsZ == null || !insideGrid(row, col))
            {
                return nullZValue;
            }
            return pointsZ[row][col];
        }
    */
    public float[] getRowGridlineValues(int nRow)
    {
        return pointsZ[nRow];
    }

    public float[] getColGridlineValues(int col)
    {
        float[] values = new float[nRows];
        for (int row = 0; row < nRows; row++)
        {
            values[row] = pointsZ[row][col];
        }
        return values;
    }

    public final float getAvgGridIncSq()
    {
        return xInc * xInc + yInc * yInc;
    }

    public boolean surfaceGridIsSame(StsSurface surface)
    {
        if (!StsMath.sameAs(surface.getXMin(), xMin))
        {
            return false;
        }
        if (!StsMath.sameAs(surface.getYMin(), yMin))
        {
            return false;
        }
        if (surface.getNRows() != nRows)
        {
            return false;
        }
        if (surface.getNCols() != nCols)
        {
            return false;
        }
        if (!StsMath.sameAs(surface.getXInc(), xInc))
        {
            return false;
        }
        if (!StsMath.sameAs(surface.getYInc(), yInc))
        {
            return false;
        }
        return true;
    }

    public void displaySurfacePosition(StsGridPoint gridPoint)
    {
        float xyz[] = gridPoint.getXYZorT();
        float x = xyz[0];
        float y = xyz[1];
        float z = xyz[2];

        double[] absXY = currentModel.getProject().getAbsoluteXYCoordinates(xyz);

        // Get ij location info
        //        StsGridPoint gridPoint = new StsGridPoint(p, this);
        //        getUnboundedGridFloatCoordinates(gridPoint);

        // display x-y-z location

        int i = gridPoint.row;
        int j = gridPoint.col;
        int n = getIndex(i, j);
        String message = null;
        String msg1 = getName() + "  - X: " + (int) (absXY[0]) + "  Y: " +
            (int) (absXY[1]) + "  Z: " + (int) z +
            "  - local X: " + (int) x + "  local Y: " + (int) y;
        String msg2 = "   - i: " + i + " j: " + j;

        String msg3 = "";
        StsRotatedGridBoundingBox rotatedBoundingBox = currentModel.getProject().getRotatedBoundingBox();
        if (rotatedBoundingBox.rowColLabelsSet())
        {
            float rowNum = rotatedBoundingBox.getNearestBoundedRowNumFromY(y);
            float colNum = rotatedBoundingBox.getNearestBoundedColNumFromX(x);
            msg3 = new String(" - Line: " + rowNum + " Crossline: " + colNum);
        }
        String msg4 = " Type: " + StsParameters.getSurfacePointTypeName(pointsNull[i][j]);
        StsMessageFiles.infoMessage(msg1 + msg2 + msg3 + msg4);
    }

    protected String getNullName(int row, int col)
    {
        if (pointsNull == null)
        {
            return StsGridPoint.nullNames[SURF_PNT];
        }
        byte nullType = pointsNull[row][col];
        return StsParameters.getSurfacePointTypeName(nullType);
    }

    public byte getNullType(int row, int col)
    {
        if (pointsNull == null)
            return StsGridPoint.SURF_PNT;
        else
            return pointsNull[row][col];
    }

    public float getZ(StsPoint point)
    {
        StsGridPoint gridPoint = new StsGridPoint(point, this);
        return getZ(gridPoint);
    }

    public float getZ(StsGridPoint gridPoint)
    {
        interpolateBilinear(gridPoint);
        return gridPoint.getZorT();
    }

    public float getInterpolatedZ(float x, float y)
    {
        StsPoint xyz = new StsPoint(x, y, 0.0f);
        StsGridPoint gridPoint = new StsGridPoint(xyz, this);
        if(!interpolateBilinear(gridPoint)) return nullValue;
        return xyz.getZorT();
    }

    public final float getInterpolatedZorT(byte zDomain, float x, float y)
    {
        StsPoint xyz = new StsPoint(x, y, 0.0f);
        StsGridPoint gridPoint = new StsGridPoint(xyz, this);
        if(!interpolateBilinear(gridPoint)) return nullValue;
        boolean isDepth = zDomain == StsProject.TD_DEPTH;
        return xyz.getZorT(isDepth);
    }

    public final StsPoint getInterpolatedPoint( float x, float y)
    {
        StsPoint xyz = new StsPoint(x, y, 0.0f);
        StsGridPoint gridPoint = new StsGridPoint(xyz, this);
        if(!interpolateBilinear(gridPoint)) return null;
        return gridPoint.point;
    }

    public final StsPoint getPoint(float rowF, float colF)
    {
        return getInterpolatedPoint(getXCoor(colF), getYCoor(rowF));
    }

    public final float[] XYZorT(float rowF, float colF)
    {
        StsPoint point = getPoint(rowF, colF);
        return point.getXYZorT();
    }

    public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint)
    {
        StsGridPoint gridPoint = new StsGridPoint(point, this);
        interpolateBilinear(gridPoint);

        float z = gridPoint.getZorT(isDepth);
        if (setPoint) point.setZ(z);
        float t = gridPoint.getZorT(isDepth);
        if (setPoint) point.setT(t);
        if (isDepth)
            return z;
        else
            return t;
    }

    public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint)
    {
        interpolateBilinear(gridPoint);
        return gridPoint.getPoint().getZorT(isDepth);
    }
    /*
        public StsGridPoint interpolateBilinear(float[] xyz)
        {
            StsGridPoint gridPoint = new StsGridPoint(xyz, this);
            if (!interpolateBilinear(gridPoint))
            {
                return null;
            }
            return gridPoint;
        }
    */
    /** Use z values at nulls, but set null flag to NULL_GAP if 3 or 4 points are null */
    /*
        private boolean interpolateBilinear(StsGridPoint gridPoint)
        {
            int i = 0, j = 0;
            int n;
            float dx, dy;
            boolean mainDebug = false;
            float z;
            //int nNulls = 0;

            try
            {
                if (pointsZ == null)
                {
                    return false;
                }

                float rowF = gridPoint.rowF;
                float colF = gridPoint.colF;
                if (!insideGrid(rowF, colF))
                {
                    gridPoint.pxyz[2] = nullZValue;
                    return false;
                }

                i = (int) rowF;
                if (i == nRows - 1)
                {
                    i--;
                }
                dy = rowF - i;

                j = (int) colF;
                if (j == nCols - 1)
                {
                    j--;
                }
                dx = colF - j;

                if (mainDebug)
                {
                    System.out.println("\txInt = " + rowF +
                                       ", dx = " + dx + ", yInt = " + colF +
                                       ", dy = " + dy);
                }

                float weight = 0.0f;
                float zWeighted = 0.0f;
                float w;
                float nullWeight = 0.0f;

                w = (1.0f - dy) * (1.0f - dx);
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i][j];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i][j] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i][j] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = dy * (1.0f - dx);
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i + 1][j];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i+1][j] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i + 1][j] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = (1.0f - dy) * dx;
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i][j + 1];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[ijE] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i][j + 1] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = dy * dx;
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i + 1][j + 1];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i+1][j+1] = " + z);
                        }
                        if (pointsNull != null &&
                            pointsNull[i + 1][j + 1] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                gridPoint.pxyz[2] = zWeighted / weight;

                if (nullWeight / weight < 0.5f) //(nNulls < 2)
                {
                    gridPoint.nullType = NOT_NULL;
                }
                else
                {
                    gridPoint.nullType = NULL_GAP;
                }

                return true;
            }
            catch (Exception e)
            {
                StsException.outputException(
                    "StsGrid.interpolateBilnearNulls() failed." +
                    " surface: " + getName() + " i: " + i + " j: " + j, e,
                    StsException.WARNING);
                gridPoint.pxyz[2] = nullZValue;
                return false;
            }
        }
    */

    /** Use z values at nulls, but set null flag to NULL_GAP if 3 or 4 points are null */
    /*
        private boolean interpolateBilinear(StsGridPoint gridPoint)
        {
            int i = 0, j = 0;
            int n;
            float dx, dy;
            boolean mainDebug = false;
            float z;
            //int nNulls = 0;

            float[][] pointsZ;
            if (!isDepth || this.adjPointsZ == null)
                pointsZ = this.pointsZ;
            else
                pointsZ = this.adjPointsZ;

            try
            {
                if (pointsZ == null)
                {
                    return false;
                }

                float rowF = gridPoint.rowF;
                float colF = gridPoint.colF;
                if (!insideGrid(rowF, colF))
                {
                    gridPoint.pxyz[2] = nullZValue;
                    return false;
                }

                i = (int) rowF;
                if (i == nRows - 1)
                {
                    i--;
                }
                dy = rowF - i;

                j = (int) colF;
                if (j == nCols - 1)
                {
                    j--;
                }
                dx = colF - j;

                if (mainDebug)
                {
                    System.out.println("\txInt = " + rowF +
                                       ", dx = " + dx + ", yInt = " + colF +
                                       ", dy = " + dy);
                }

                float weight = 0.0f;
                float zWeighted = 0.0f;
                float w;
                float nullWeight = 0.0f;

                w = (1.0f - dy) * (1.0f - dx);
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i][j];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i][j] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i][j] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = dy * (1.0f - dx);
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i + 1][j];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i+1][j] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i + 1][j] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = (1.0f - dy) * dx;
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i][j + 1];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[ijE] = " + z);
                        }
                        if (pointsNull != null && pointsNull[i][j + 1] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                w = dy * dx;
                if (w > StsParameters.roundOff)
                {
                    z = pointsZ[i + 1][j + 1];
                    if (z == nullZValue)
                    {
                        nullWeight += w; //nNulls++;
                    }
                    else
                    {
                        weight += w;
                        zWeighted += w * z;
                        if (mainDebug)
                        {
                            System.out.println("\tz[i+1][j+1] = " + z);
                        }
                        if (pointsNull != null &&
                            pointsNull[i + 1][j + 1] != NOT_NULL)
                        {
                            nullWeight += w; //nNulls++;
                        }
                    }
                }

                gridPoint.pxyz[2] = zWeighted / weight;

                if (nullWeight / weight < 0.5f) //(nNulls < 2)
                {
                    gridPoint.nullType = NOT_NULL;
                }
                else
                {
                    gridPoint.nullType = NULL_GAP;
                }

                return true;
            }
            catch (Exception e)
            {
                StsException.outputException(
                    "StsGrid.interpolateBilnearNulls() failed." +
                    " surface: " + getName() + " i: " + i + " j: " + j, e,
                    StsException.WARNING);
                gridPoint.pxyz[2] = nullZValue;
                return false;
            }
        }
    */
    private boolean interpolateBilinear(StsGridPoint gridPoint)
    {
        int i = 0, j = 0;
        int n;
        float dx, dy;
        boolean debug = false;
        float z;
        float zAdj;
        //int nNulls = 0;

        try
        {
            if (pointsZ == null)
            {
                return false;
            }

            boolean doAdjPoint = adjPointsZ != null;
            //            boolean doAdjPoint = isDepth && (adjPointsZ != null);

            float rowF = gridPoint.rowF;
            float colF = gridPoint.colF;
            if (!insideGrid(rowF, colF))
            {
                gridPoint.setZ(nullZValue);
                gridPoint.nullType = SURF_GAP;
                return false;
            }

            i = (int) rowF;
            if (i == nRows - 1)
            {
                i--;
            }
            dy = rowF - i;

            j = (int) colF;
            if (j == nCols - 1)
            {
                j--;
            }
            dx = colF - j;

            if (debug)
            {
                System.out.println("\txInt = " + rowF +
                    ", dx = " + dx + ", yInt = " + colF +
                    ", dy = " + dy);
            }

            float weight = 0.0f;
            float zWeighted = 0.0f;
            float w;
            float nullWeight = 0.0f;
            float zAdjWeighted = 0.0f;

            w = (1.0f - dy) * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i][j];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (doAdjPoint) zAdjWeighted += w * adjPointsZ[i][j];

                    if (debug)
                    {
                        System.out.println("\tz[i][j] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i][j] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = dy * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i + 1][j];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (doAdjPoint) zAdjWeighted += w * adjPointsZ[i + 1][j];
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i + 1][j] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = (1.0f - dy) * dx;
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i][j + 1];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (doAdjPoint) zAdjWeighted += w * adjPointsZ[i][j + 1];
                    if (debug)
                    {
                        System.out.println("\tz[ijE] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i][j + 1] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = dy * dx;
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i + 1][j + 1];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (doAdjPoint) zAdjWeighted += w * adjPointsZ[i + 1][j + 1];
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j+1] = " + z);
                    }
                    if (pointsNull != null &&
                        pointsNull[i + 1][j + 1] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            z = zWeighted / weight; // time or approx depth
            if (doAdjPoint)
            {
                zAdj = zAdjWeighted / weight; // adjusted seismic depth
                gridPoint.getPoint().setZ(zAdj);
                gridPoint.getPoint().setT(z);
            }
            else
                gridPoint.getPoint().setZorT(z);

            if (nullWeight / weight < 0.5f) //(nNulls < 2)
            {
                gridPoint.nullType = SURF_PNT;
            }
            else
            {
                gridPoint.nullType = SURF_GAP;
            }

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException(
                "StsGrid.interpolateBilnearNulls() failed." +
                    " surface: " + getName() + " i: " + i + " j: " + j, e,
                StsException.WARNING);
            gridPoint.setZ(nullZValue);
            return false;
        }
    }

    public float interpolateBilinearNoNulls(float rowF, float colF)
    {
        return interpolateBilinearNoNulls(pointsZ, rowF, colF);
    }

    public float interpolateDepthNoNulls(float rowF, float colF)
    {
        if (adjPointsZ == null) return nullValue;
        return interpolateBilinearNoNulls(adjPointsZ, rowF, colF);
    }

    public float interpolateTimeNoNulls(float rowF, float colF)
    {
        if (pointsZ == null) return nullValue;
        return interpolateBilinearNoNulls(pointsZ, rowF, colF);
    }

    public float interpolateBilinearNoNulls(float[][] pointsZ, float rowF, float colF)
    {
        float z;
        int i = 0, j = 0;
        int n;
        float dx, dy;
        boolean debug = false;
        //int nNulls = 0;

        try
        {
            if (pointsZ == null)
            {
                return nullZValue;
            }

            if (!insideGrid(rowF, colF))
            {
                return nullZValue;
            }

            i = (int) rowF;
            if (i == nRows - 1)
            {
                i--;
            }
            dy = rowF - i;

            j = (int) colF;
            if (j == nCols - 1)
            {
                j--;
            }
            dx = colF - j;

            if (debug)
            {
                System.out.println("\txInt = " + rowF +
                    ", dx = " + dx + ", yInt = " + colF +
                    ", dy = " + dy);
            }

            float weight = 0.0f;
            float zWeighted = 0.0f;
            float w;
            float nullWeight = 0.0f;

            w = (1.0f - dy) * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i][j];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i][j] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i][j] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = dy * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i + 1][j];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i + 1][j] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = (1.0f - dy) * dx;
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i][j + 1];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[ijE] = " + z);
                    }
                    if (pointsNull != null && pointsNull[i][j + 1] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }

            w = dy * dx;
            if (w > StsParameters.roundOff)
            {
                z = pointsZ[i + 1][j + 1];
                if (z == nullZValue)
                {
                    nullWeight += w; //nNulls++;
                }
                else
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j+1] = " + z);
                    }
                    if (pointsNull != null &&
                        pointsNull[i + 1][j + 1] != SURF_PNT)
                    {
                        nullWeight += w; //nNulls++;
                    }
                }
            }
            if (weight > 0.0f)
            {
                return zWeighted / weight;
            }
            else
            {
                return nullZValue;
            }
        }
        catch (Exception e)
        {
            StsException.outputException(
                "StsGrid.interpolateBilnearNulls() failed." +
                    " surface: " + getName() + " i: " + i + " j: " + j, e,
                StsException.WARNING);
            return nullZValue;
        }
    }

    private void interpolateBilinearNoNulls(StsGridPoint gridPoint)
    {
        int i = 0, j = 0;
        float dx, dy;
        boolean debug = false;

        try
        {
            i = gridPoint.row;
            dy = gridPoint.rowF - i;
            j = gridPoint.col;
            dx = gridPoint.colF - j;

            if (debug)
            {
                System.out.println("\txInt = " + gridPoint.rowF +
                    ", dx = " + dx + ", yInt = " +
                    gridPoint.colF +
                    ", dy = " + dy);
            }

            float weight = 0.0f;
            float z;
            float zWeighted = 0.0f;
            float w;

            w = (1.0f - dy) * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = getZNotNull(i, j);
                if (z != nullZValue)
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i][j] = " + z);
                    }
                }
            }
            w = dy * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = getZNotNull(i + 1, j);
                if (z != nullZValue)
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j] = " + z);
                    }
                }
            }
            w = (1.0f - dy) * dx;
            if (w > StsParameters.roundOff)
            {
                z = getZNotNull(i, j + 1);
                if (z != nullZValue)
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[ijE] = " + z);
                    }
                }
            }
            w = dy * dx;
            if (w > StsParameters.roundOff)
            {
                z = getZNotNull(i + 1, j + 1);
                if (z != nullZValue)
                {
                    weight += w;
                    zWeighted += w * z;
                    if (debug)
                    {
                        System.out.println("\tz[i+1][j+1] = " + z);
                    }
                }
            }
            if (weight > 0.0f)
            {
                gridPoint.setZ(zWeighted / weight);
            }
            else
            {
                gridPoint.setZ(nullZValue);
            }
        }
        catch (Exception e)
        {
            StsException.outputException(
                "StsGrid.interpolateBilnearNoNulls() failed." +
                    " surface: " + getName() + " i: " + i + " j: " + j, e,
                StsException.WARNING);
            gridPoint.setZ(nullZValue);
        }
    }

    protected int getIndex(int row, int col)
    {
        return row * nCols + col;
    }

    public StsPoint getPoint(int row, int col)
    {
        if (pointsZ == null || !insideGrid(row, col))
        {
            return null;
        }

        return new StsPoint(getXCoor(row, col), getYCoor(row, col), pointsZ[row][col]);
    }

    public final StsPoint getStsPoint(float rowF, float colF)
    {
        return getInterpolatedPoint(getXCoor(colF), getYCoor(rowF));
    }

    public final float[] getXYZorT(float rowF, float colF)
    {
        StsPoint point = getStsPoint(rowF, colF);
        return point.getXYZorT();
    }

    public StsPoint getPointZandT(int row, int col)
    {
        if (pointsZ == null || !insideGrid(row, col))
        {
            return null;
        }
        StsPoint point = new StsPoint(6);
        point.setX(getXCoor(row, col));
        point.setY(getYCoor(row, col));
        point.setT(pointsZ[row][col]);

        if (adjPointsZ != null)
            point.setZ(adjPointsZ[row][col]);
        else
            point.setZ(pointsZ[row][col]);

        return point;
    }

    public StsPoint getPointZorT(int row, int col)
    {
        if (pointsZ == null || !insideGrid(row, col))
            return null;
        StsPoint point = new StsPoint(3);
        point.setX(getXCoor(row, col));
        point.setY(getYCoor(row, col));
        point.setZ(getZorT(row, col));
        return point;
    }

    // for StsXYSurfaceGridable compatibility
    public float getComputePointZ(int row, int col)
    {
        return pointsZ[row][col];
    }

    public float getPointZ(int row, int col)
    {
        return pointsZ[row][col];
    }

    public void setPointFilled(int row, int col, float z)
    {
        pointsZ[row][col] = z;
        setPointNull(row, col, SURF_GAP_FILLED);
    }

    public void setPointFilled(int row, int col, float z, byte fillFlag)
    {
        pointsZ[row][col] = z;
        setPointNull(row, col, fillFlag);
    }

    public float[][][] getComputeNormals()
    {
        if (normals != null)
        {
            return normals;
        }
        computeNormals();
        return normals;
    }

    public float[] getNormal(int row, int col)
    {
        if (!insideGrid(row, col))
        {
            return null;
        }
        if (normals == null)
        {
            return null;
        }
        return normals[row][col];
    }

    public float[] getNormal(float rowF, float colF)
    {
        float[] norm0, norm1, norm2;

        try
        {
            int row = (int) rowF;
            float dR = rowF - row;

            int col = (int) colF;
            float dC = colF - col;

            if (dR == 0.0f)
            {
                if (dC == 0.0f)
                {
                    return normals[row][col];
                }
                else
                {
                    return StsMath.interpolate(normals[row][col],
                        normals[row][col + 1], dC);
                }
            }
            else
            {
                norm0 = StsMath.interpolate(normals[row][col],
                    normals[row + 1][col], dR);
                if (dC == 0.0f)
                {
                    return norm0;
                }
                else
                {
                    norm1 = StsMath.interpolate(normals[row][col + 1],
                        normals[row + 1][col + 1], dR);
                    return StsMath.interpolate(norm0, norm1, dC);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return null;
        }
    }

    public void setPointsNull(byte nullType)
    {
        byte currentNullType = 0;
        if (pointsNull == null) return;
        if (Main.debugPoint)
            currentNullType = pointsNull[Main.debugPointRow][Main.debugPointCol];
        for (int row = 0; row < nRows; row++)
            Arrays.fill(pointsNull[row], nullType);
        if (Main.debugPoint)
        {
            StsException.systemDebug(this, "setPointsNull",
                getName() + "pointsNull[" + Main.debugPointRow + "][" + Main.debugPointCol + "] changed to " + StsParameters.getSurfacePointTypeName(nullType) +
                    " from " + StsParameters.getSurfacePointTypeName(currentNullType));
        }
    }

    public void setPointNull(int row, int col, byte nullType)
    {
        if (pointsNull == null)
        {
            return;
        }

        if (Main.debugPoint && row == Main.debugPointRow && col == Main.debugPointCol)
        {
            byte currentNullType = pointsNull[row][col];
            StsException.systemDebug(this, "setPointNull",
                getName() + "pointsNull[" + row + "][" + col + "] changed to " + StsParameters.getSurfacePointTypeName(nullType) +
                    " from " + StsParameters.getSurfacePointTypeName(currentNullType));
        }
        pointsNull[row][col] = nullType;
    }

    public void setMarker(StsMarker m)
    {
        dbFieldChanged("marker", m);
        marker = m;
    }

    public StsMarker getMarker()
    {
        return marker;
    }

    public String getMarkerName()
    {
        return marker == null ? "none" : marker.getName();
    }

    public void setMarkerName(String name)
    {
        if (marker != null)
        {
            marker.setName(name);
        }
    }

    public void initializeDepthFromTime()
    {
        adjPointsZ = new float[nRows][nCols];
        this.dbFieldChanged("adjPointsZ", adjPointsZ);
        this.fieldChanged("zDomainSupported", StsProject.TD_TIME_DEPTH);
    }

    public float getAdjPointZ(int row, int col)
    {
        if (adjPointsZ == null) return nullValue;
        return adjPointsZ[row][col];
    }

    /*
        public StsPoint XcomputeWellIntersect(StsWell well)
        {
            return XcomputeWellIntersect(well, false);
        }

        public StsPoint XcomputeWellIntersect(StsWell well,
                                              boolean requireIntersect)
        {
            int n;
            float error, lastError;
            float topPointZ, botPointZ;
            float surfaceZ;
            float f;
            StsPoint wellPoint;
            float wellPointZ, lastWellPointZ;

            StsPoint[] wellPoints = well.getPoints();
            int noWellPoints = wellPoints.length;

            if (mainDebug)
            {
                System.out.println("Well: " + well.index() + " top-down search.");
            }
            topDownLoop:
                {
                lastError = nullZValue;
                lastWellPointZ = wellPoints[0].getZ();

                for (n = 0; n < noWellPoints; n++)
                {
                    wellPoint = wellPoints[n];
                    surfaceZ = interpolateBilinearZ(wellPoint, false, false);
                    if (mainDebug)
                    {
                        System.out.println(" z: " + surfaceZ);
                    }
                    if (surfaceZ == nullZValue)
                    {
                        continue;
                    }

                    wellPointZ = wellPoint.getZ();
                    error = surfaceZ - wellPointZ;

                    if (error == 0.0f)
                    {
                        topPointZ = wellPointZ;
                        break topDownLoop;
                    }

                    if (lastError != nullZValue)
                    {
                        if (lastError < 0.0f && error > 0.0f ||
                            lastError > 0.0f && error < 0.0f)
                        {
                            f = Math.abs(lastError / (lastError - error));
                            topPointZ = lastWellPointZ +
                                (wellPointZ - lastWellPointZ) * f;
                            break topDownLoop;
                        }
                    }

                    lastError = error;
                    lastWellPointZ = wellPointZ;
                }

                // Well didn't intersect map: use top or bottom point whichever is nearest surface
                if (requireIntersect)
                {
                    return null;
                }

                wellPoint = wellPoints[0];
                surfaceZ = interpolateBilinearZ(wellPoint, false, false);
                float topError = Math.abs(surfaceZ - wellPoint.getZ());

                wellPoint = wellPoints[noWellPoints - 1];
                surfaceZ = interpolateBilinearZ(wellPoint, false, false);
                if (mainDebug)
                {
                    System.out.println(" z: " + surfaceZ);
                }
                float botError = Math.abs(surfaceZ - wellPoint.getZ());

                if (topError < botError)
                {
                    return wellPoints[0];
                }
                else
                {
                    return wellPoints[noWellPoints - 1];
                }
            }

            bottomUpLoop:
                {
                lastError = nullZValue;
                lastWellPointZ = wellPoints[noWellPoints - 1].getZ();

                if (mainDebug)
                {
                    System.out.println("Well: " + well.index() +
                                       " bottom-up search.");
                }
                for (n = noWellPoints - 1; n >= 0; n--)
                {
                    wellPoint = wellPoints[n];
                    surfaceZ = interpolateBilinearZ(wellPoint, false, false);
                    if (surfaceZ == nullZValue)
                    {
                        continue;
                    }
                    wellPointZ = wellPoint.getZ();

                    error = surfaceZ - wellPointZ;

                    if (error == 0.0f)
                    {
                        botPointZ = wellPointZ;
                        break bottomUpLoop;
                    }

                    if (lastError != nullZValue)
                    {
                        if (lastError < 0.0f && error > 0.0f ||
                            lastError > 0.0f && error < 0.0f)
                        {
                            f = Math.abs(lastError / (lastError - error));
                            botPointZ = lastWellPointZ +
                                (wellPointZ - lastWellPointZ) * f;
                            break bottomUpLoop;
                        }
                    }

                    lastError = error;
                    lastWellPointZ = wellPointZ;
                }

                // Didn't find botPointZ, but topPointZ is OK
                return well.getPointAtZ(topPointZ, false);
            }

            float avgZ = 0.5f * (topPointZ + botPointZ);
            return well.getPointAtZ(avgZ, true);
        }
    */
    public void initializeGridBoundaryPoints(StsModel model, StsEdgeLoopRadialLinkGrid linkedGrid)
    {
        int row, col;

        if (linkedGrid == null)
        {
            return;
        }

        StsGridIterator gridIterator = new StsGridIterator(linkedGrid);

        for (row = 0; row < nRows; row++)
        {
            if (!gridIterator.hasRow(row))
            {
                for (col = 0; col < nCols; col++)
                {
                    setPointNull(row, col, SURF_BOUNDARY);
                }
            }
            else
            {
                int colMin = gridIterator.colMin;
                int colMax = gridIterator.colMax;

                if (colMax > colMin)
                {
                    for (col = 0; col < colMin; col++)
                    {
                        setPointNull(row, col, SURF_BOUNDARY);
                    }

                    for (col = colMax + 1; col < nCols; col++)
                    {
                        setPointNull(row, col, SURF_BOUNDARY);
                    }
                }
            }
        }
    }

    private boolean hasNeighborNSEW(int row, int col)
    {
        if (row > 0 && !isPointNull(row - 1, col))
        {
            return true;
        }
        if (row < nRows - 1 && !isPointNull(row + 1, col))
        {
            return true;
        }
        if (col > 0 && !isPointNull(row, col - 1))
        {
            return true;
        }
        if (col < nCols - 1 && !isPointNull(row, col + 1))
        {
            return true;
        }
        return false;
    }

    public void copyPointsZWithOffset(float[][] otherPointsZ, float offset)
    {
        pointsZ = new float[nRows][nCols];
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                pointsZ[row][col] = otherPointsZ[row][col] + offset;
    }

    public void copyPointsZWithOffset(float[][] otherPointsZ, float[][] otherAdjPointsZ, StsSeismicVelocityModel velocityModel, float offset)
    {
        pointsZ = new float[nRows][nCols];
        adjPointsZ = new float[nRows][nCols];
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
            {
                float t = otherPointsZ[row][col];
                float z = otherAdjPointsZ[row][col] + offset;
                adjPointsZ[row][col] = z;
                float x = getXCoor(row, col);
                float y = getYCoor(row, col);
                try
                {
                    t = (float)velocityModel.getT(x, y, z, t);
                }
                catch(Exception e) { }
                pointsZ[row][col] = t;
            }
    }

    public void copyPointsNull(byte[][] otherPointsNull)
    {
        pointsNull = new byte[nRows][nCols];
        for(int row = 0; row < nRows; row++)
            System.arraycopy(otherPointsNull[row], 0, pointsNull[row], 0, nCols);
    }
    /** returns true if congruent */
    public boolean copyPoints(StsSurface surface)
    {
        int i, j, n, ii, jj;
        float x, y, z;

        StsStatusArea statusArea = currentModel.win3d.statusArea;

        initializePoints();
        boolean gridIsSame = xyGridSameAs(surface);
        //TODO also check if they are congruent but not same (voxelGrids match, but range is different) and act accordingly
        if (gridIsSame)
        {
            StsMessageFiles.logMessage("Copying grid points...");

            // should check that offsets are integral
            int rowMin, rowMax, otherRowMin = 0;
            int colMin, colMax, otherColMin = 0;

            float rowMinF = getRowCoor(surface.yMin);
            rowMin = Math.round(rowMinF); // start row on modelSurface for first point of other surface
            int rowOffset = rowMin;
            if (Math.abs(rowMin - rowMinF) > 0.1f)
            {
                System.out.println("Not exactly concruent for rowMinF: " + rowMinF);
            }
            if (rowMin < 0)
            {
                otherRowMin = -rowMin;
                rowMin = 0;
            }
            else
            {
                otherRowMin = 0;
            }

            float colMinF = getColCoor(surface.xMin);
            colMin = Math.round(colMinF);
            if (Math.abs(colMin - colMinF) > 0.1f)
            {
                System.out.println("Not exactly concruent for colMinF: " +
                    colMinF);
            }
            if (colMin < 0)
            {
                otherColMin = -colMin;
                colMin = 0;
            }
            else
            {
                otherColMin = 0;
            }

            float rowMaxF = getRowCoor(surface.yMax); // end row on modelSurface for last point of other surface
            rowMax = Math.round(rowMaxF);
            if (Math.abs(rowMax - rowMaxF) > 0.1f)
            {
                System.out.println("Not exactly concruent for rowMaxF: " +
                    rowMaxF);
            }
            rowMax = Math.min(rowMax, nRows - 1);

            float colMaxF = getColCoor(surface.xMax);
            colMax = Math.round(colMaxF);
            if (Math.abs(colMax - colMaxF) > 0.1f)
            {
                System.out.println("Not exactly concruent for colMaxF: " +
                    colMaxF);
            }
            colMax = Math.min(colMax, nCols - 1);

            float[][] otherPointsZ = surface.getPointsZ();
            byte[][] otherPointsNull = surface.getPointsNull();
            zMin = otherPointsZ[rowMin][colMin];
            zMax = zMin;
            for (i = rowMin, n = 0, ii = otherRowMin; i <= rowMax; i++, ii++)
            {
                for (j = colMin, jj = otherColMin; j <= colMax; j++, jj++,
                    n++)
                {
                    z = otherPointsZ[ii][jj];
                    pointsZ[i][j] = z;
                    pointsNull[i][j] = otherPointsNull[ii][jj];

                    if (z < zMin)
                        zMin = z;
                    else if (z > zMax)
                        zMax = z;
                }
            }
        }
        else // not gridIsSame: interpolate to new grid
        {
            StsMessageFiles.logMessage("Interpolating grid points...");
            statusArea.addProgress();
            statusArea.setMaximum(nRows);

            StsGridPoint gridPoint = new StsGridPoint(this);

            y = yMin;
            for (i = 0, n = 0; i < nRows; i++)
            {
                gridPoint.setY(y, surface);
                x = xMin;
                for (j = 0; j < nCols; j++, n++)
                {
                    gridPoint.setX(x, surface);
                    surface.interpolateBilinear(gridPoint);
                    z = gridPoint.getZorT();
                    pointsZ[i][j] = z;
                    if (hasNulls)
                    {
                        setPointNull(i, j, gridPoint.nullType);
                        if (gridPoint.nullType == SURF_PNT)
                        {
                            if (z < zMin)
                            {
                                zMin = z;
                            }
                            if (z > zMax)
                            {
                                zMax = z;
                            }
                        }
                    }
                    else
                    {
                        setPointNull(i, j, SURF_PNT);

                        if (z < zMin)
                        {
                            zMin = z;
                        }
                        if (z > zMax)
                        {
                            zMax = z;
                        }
                    }
                    x += xInc;
                }
                statusArea.setProgress(i + 1);
                y += yInc;
            }
            statusArea.removeProgress();
        }
        // if time, values re stored in rotatedBoundingBox.zMin, zMax
        if (this.zDomainOriginal == StsParameters.TD_DEPTH)
        {
            depthMin = zMin;
            depthMax = zMax;
        }
        return gridIsSame;
    }

    /*
       public void persistPoints()
       {
           setPointsZ(pointsZ);
           setPointsNull(pointsNull);
       }
    */
    public void clipToBoundary(int rowMin, int rowMax, int colMin, int colMax)
    {
        int i, j;

        for (i = 0; i < rowMin; i++)
        {
            for (j = 0; j < nCols; j++)
            {
                setPointNull(i, j, SURF_BOUNDARY);
            }
        }

        for (i = rowMax + 1; i < nRows; i++)
        {
            for (j = 0; j < nCols; j++)
            {
                setPointNull(i, j, SURF_BOUNDARY);
            }
        }

        for (i = rowMin; i <= rowMax; i++)
        {
            for (j = 0; j < colMin; j++)
            {
                setPointNull(i, j, SURF_BOUNDARY);
            }

            for (j = colMax + 1; j < nCols; j++)
            {
                setPointNull(i, j, SURF_BOUNDARY);
            }
        }
    }

    public void computeGradients()
    {
        int i, j;

        dZdX = new float[nRows][nCols];
        dZdY = new float[nRows][nCols];

        for (i = 0; i < nRows; i++)
        {
            for (j = 0; j < nCols; j++)
            {
                if (getPointNull(i, j) == SURF_PNT)
                {
                    computeGradients(i, j);
                }
            }
        }
    }

    public void computeGradients(int i, int j)
    {
        float z, zp = 0, zm = 0;

        boolean debug = Main.debugPoint && i == Main.debugPointRow && j == Main.debugPointCol;

        z = getZorT(i, j);

        if (getPointNull(i - 1, j) == SURF_PNT)
        {
            zm = getZorT(i - 1, j);

            if (getPointNull(i + 1, j) == SURF_PNT)
            {
                zp = getZorT(i + 1, j);
                dZdY[i][j] = (zp - zm) / 2;
            }
            else
            {
                dZdY[i][j] = z - zm;
            }
        }
        else if (getPointNull(i + 1, j) == SURF_PNT)
        {
            zp = getZorT(i + 1, j);
            dZdY[i][j] = zp - z;
        }
        else
        {
            dZdY[i][j] = 0.0f;
        }

        if (getPointNull(i, j - 1) == SURF_PNT)
        {
            zm = getZorT(i, j - 1);

            if (getPointNull(i, j + 1) == SURF_PNT)
            {
                zp = getZorT(i, j + 1);
                dZdX[i][j] = (zp - zm) / 2;
            }
            else
            {
                dZdX[i][j] = z - zm;
            }
        }
        else if (getPointNull(i, j + 1) == SURF_PNT)
        {
            zp = getZorT(i, j + 1);
            dZdX[i][j] = zp - z;
        }
        else
        {
            dZdX[i][j] = 0.0f;
        }

        if(debug) StsException.systemDebug(this, "computeGradients", getName() + " i: " + i + " j: " + j + " z: " + z + " zm: " + zm + " zp: " + zp);
    }

    public float dZdX(int row, int col)
    {
        return dZdX[row][col];
    }

    public float dZdY(int row, int col)
    {
        return dZdY[row][col];
    }

    public boolean delete()
    {
        if (!super.delete())
        {
            return false;
        }
        deleteTransientArrays();
        return true;
    }

    public void deleteTransientArrays()
    {
        weights = null;
        dZdX = null;
        dZdY = null;
    }

    public void truncateBySurfaceAbove(StsSurface surfaceAbove)
    {
        float[][] zArrayAbove = surfaceAbove.getPointsZ();
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                float pointAbove = zArrayAbove[row][col];
                if (pointsZ[row][col] < pointAbove)
                {
                    pointsZ[row][col] = pointAbove + 0.1f;
                }
            }
        }
    }

    public float getMinColorscaleAlpha()
    {
        if (newSurfaceTexture == null) return 1.0f;
        if (newSurfaceTexture.getColorscale() == null) return 1.0f;
        float[][] rgbs = newSurfaceTexture.getColorscale().computeRGBAArray(false);
        if (rgbs == null) return 1.0f;
        float minAlpha = 1.0f;
        for (int i = 0; i < rgbs[0].length; i++)
        {
            if (rgbs[3][i] > 0.f)
                minAlpha = Math.min(minAlpha, rgbs[3][i]);
        }
        if (minAlpha < 0.1f) minAlpha = 0.1f;
        return minAlpha;
    }

    /** display this surface */
    protected void displaySurface(StsGLPanel3d glPanel3d)
    {
        if (pointsZ == null)
        {
            return;
        }
        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        if (relativeRotationAngle != 0.0f)
        {
            gl.glPushMatrix();
            gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
        }
        if (offset != 0.0f)
        {
            gl.glPushMatrix();
            gl.glTranslatef(0.0f, 0.0f, offset);
        }

        if (isDisplayFill())
        {
         //   stsColor.setGLColor(gl);
            //            checkSetTextureDisplayable(gl);
            //			System.out.println("lighting "+lighting);
            if (lighting)
            {
                gl.glEnable(GL.GL_LIGHTING);
                //                glPanel3d.getCurrentView().setMaterialTransparency(false);
            }
            float minAlpha = getMinColorscaleAlpha();
            gl.glAlphaFunc(GL.GL_GEQUAL, minAlpha); // jbw to minimize fade to black on missings
            gl.glFrontFace(GL.GL_CW); // change to left-handed coordinate system (z is down)
            //			gl.glFrontFace(GL.GL_CW); // change to left-handed coordinate system (z is down)
            displaySurfaceFill(glPanel3d);
            gl.glAlphaFunc(GL.GL_GEQUAL, 0.1f);
            gl.glFrontFace(GL.GL_CCW); //
            gl.glDisable(GL.GL_LIGHTING);
        }
        if (displayGrid)
        {
            StsColor gridColor = displayFill ? StsColor.BLACK : stsColor;
            gridColor.setGLColor(gl);
            gl.glDisable(GL.GL_LIGHTING);
            displayGridLines(glPanel3d);
        }
        if (relativeRotationAngle != 0.0f)
        {
            gl.glPopMatrix();
        }
        if (offset != 0.0f)
        {
            gl.glPopMatrix();
        }
        boolean debug = currentModel.getBooleanProperty("Surface Gap Points");
        if(debug) displayPoints(currentModel, glPanel3d);
    }

    public void displayPoints(StsModel model, StsGLPanel3d glPanel3d)
    {
        try
        {
            float[] point = new float[3];
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);

            glPanel3d.setViewShift(gl, 4.0*StsGraphicParameters.gridShift);

            for(int row = 0; row < nRows; row++)
            {
                for(int col = 0; col < nCols; col++)
                {
                    point[0] = getXCoor(row, col);
                    point[1] = getYCoor(row, col);
                    byte pointType = pointsNull[row][col];
                    point[2] = this.getZorT(row, col);
                    displayPoint(point, pointType, glPanel3d);
                }
            }

            gl.glEnable(GL.GL_LIGHTING);
            glPanel3d.resetViewShift(gl);
        }
        catch(Exception e){} // mainDebug, so don't process exception
    }
    /*
        static public final byte NOT_NULL = 1; // good
        static public final byte NULL_GAP_SET = 2; // nulled by gapping routine
        static public final byte NULL_GAP = 3; // originally null
        static public final byte NULL_BOUNDARY = 4; // originally null and is outside boundary
        static public final byte NULL_GAP_FILLED = 5; // NULL_GAP type filled in.
        static public final byte NULL_GAP_NOT_FILLED = 7; // Needs to be filled by interpolation
    */
    // Debug display of point status: BLUE - good grid point (NOT_NULL);
    //                              PURPLE - interpolated point not yet filled (NULL_GAP_SET);
    //                                CYAN - interpolated point (NULL_GAP);
    //                                 RED - extrapolated point (NULL_BOUNDARY);
    //                              ORANGE - point is filled hole (NULL_GAP_FILLED).
    //                             MAGENTA - extrapolated point beyond edge: cut (NULL_GAP_NOT_FILLED);
    //                                GREY - default color.

    protected StsColor getPointTypeColor(byte pointType) { return StsParameters.getSurfacePointTypeColor(pointType); }

    protected void displayPoint(float[] point, byte pointType, StsGLPanel3d glPanel3d)
    {
        StsGLDraw.drawPoint(point, getPointTypeColor(pointType), glPanel3d, 4, 4, 0.0);
    }

    private boolean isDisplayFill()
    {
        return this.displayFill; //  || newSurfaceTexture != surfaceDefaultTexture;
    }

    /*
       private void checkSetTextureDisplayable(GL gl)
       {
           if(actionTextureDisplayable != null)
           {
               if(newTextureDisplayable == actionTextureDisplayable) return;
               newTextureDisplayable = actionTextureDisplayable;
               textureChanged = true;
               return;
           }
           else if(displayProperty)
           {
               if(propertyTextureDisplayable == null)
               {
//                propertyTextureDisplayable = currentModel.getCurrentPropertyTextureDisplayable();
    if(newTextureDisplayable == propertyTextureDisplayable) return;
                   newTextureDisplayable = propertyTextureDisplayable;
               }
           }
           else // !displayProperty
           {
               if(newTextureDisplayable != null)
               {
                   newTextureDisplayable = null;
                   textureChanged = true;
               }
           }
       }
    */
    /*
    private void displaySurfaceFill(StsGLPanel3d glPanel3d)
    {
     //timer.start();

     GL gl = glPanel3d.getGL();
     if(gl == null) return;

           if(getLighting() != isLighted)
           {
               deleteSurfaceDisplayList = true;
               isLighted = !isLighted;
           }

           if(deleteSurfaceDisplayList)
           {
               deleteSurfaceDisplayList(gl);
               deleteSurfaceDisplayList = false;
           }

     boolean useDisplayLists = currentModel.getBooleanProperty("Use Display Lists");
     if (surfDisplayListNum == 0 && useDisplayLists)  // build display list
     {
      surfDisplayListNum = gl.glGenLists(1);
      if(surfDisplayListNum == 0)
      {
       currentModel.logMessage("System Error in StsGrid.displaySurface: " +
               "Failed to allocate a display list");
       return;
      }

      gl.glNewList(surfDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
      drawSurfaceFill(gl);
      gl.glEndList();

      // We are using normals only once with displayList, so null it since we are finished with it.
      normals = null;

      //timer.stop("display list surface setup: ");

     }
     else if(useDisplayLists) // use existing display list
     {
      gl.glCallList( surfDisplayListNum );
       //timer.stop("display list surface draw: ");
     }
     else
     {
               deleteSurfaceDisplayList(gl);
      drawSurfaceFill(gl);
     }
    }
    */
    public void initializeSurface()
    {
        textureChanged = true;
        geometryChanged = true;
    }

    public void saveSurface()
    {
        textureChanged = true;
        geometryChanged = true;
        dbFieldChanged("pointsZ", pointsZ);
        dbFieldChanged("pointsNull", pointsNull);
        //        interpolateNullZPoints();
    }

    public String getNativeVerticalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeVerticalUnits];
    }

    public byte getNativeVerticalUnits()
    {
        return nativeVerticalUnits;
    }

    public String getNativeHorizontalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeHorizontalUnits];
    }

    public byte getNativeHorizontalUnits()
    {
        return nativeHorizontalUnits;
    }

    /**
     * This puts texture display on delete list.  Operation is performed
     * at beginning of next draw operation.
     */
    /*
         public void addTextureToDeleteList()
         {
             if (textureTiles != null)
             {
                 StsTextureList.addTextureToDeleteList(this);
             }
             textureChanged = true;
         }
    */
    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if (textureTiles == null) return;
        textureTiles.deleteTextures(gl);
        textureChanged = true;
    }

    private void displaySurfaceFill(StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        if (pointsZ == null) return;
        /*
           if(debug)
           {
               System.out.println("textureChanged: " + textureChanged + " colorListChanged: " + colorListChanged);
               debugPrintTextureDisplayable("newTextureDisplayable", newTextureDisplayable);
               debugPrintTextureDisplayable("currentTextureDisplayable", currentTextureDisplayable);
           }
        */
        // make sure this newTextureDisplayable is still displayable
        checkSurfaceTexture();
        // if texture has changed, delete the old one
        if (newSurfaceTexture != currentSurfaceTexture)
        {
            deleteTexturesAndDisplayLists(gl);
            currentSurfaceTexture = newSurfaceTexture;
            colorListChanged = true;
        }
        if(currentSurfaceTexture == null) return;

        if (!initializeTextureTiles(glPanel3d, gl)) return;
        // currently because we are using texture for even lighted surfaces, any colorList change requires a rebuild of textures
        if (colorListChanged) textureChanged = true;

        if (textureChanged)
        {
            deleteTexturesAndDisplayLists(gl);
            if (currentSurfaceTexture != null)
            {
                colorListChanged = true;
            }
//            textureChanged = false;
        }
        if (geometryChanged || currentModel.useDisplayLists && !usingDisplayLists)
        {
            if (geometryChanged)
            {
                doCropChanged();
                geometryChanged = false;
            }
            usingDisplayLists = currentModel.useDisplayLists;
            textureTiles.constructSurface(this, gl, currentModel.useDisplayLists, true);
        }

        if (useTextures)
        {
            gl.glEnable(GL.GL_BLEND);
            gl.glEnable(GL.GL_TEXTURE_2D);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);

            if (transparency)
            {
                gl.glEnable(GL.GL_POLYGON_STIPPLE);
                gl.glPolygonStipple(StsGraphicParameters.getNextStipple(), 0);
            }

            if (lighting)
            {
                gl.glEnable(GL.GL_LIGHTING);
                gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL, GL.GL_SEPARATE_SPECULAR_COLOR);
            }
            else
            {
                gl.glDisable(GL.GL_LIGHTING);
            }

            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_SMOOTH);

                if (colorListChanged)
                {
                    boolean nullsFilled = getSurfaceClass().getNullsFilled();
                    if (currentSurfaceTexture != null)
                    {
                        textureColorListNum = currentSurfaceTexture.getColorDisplayListNum(gl, nullsFilled);
                    }
                    //            else if(hasNulls)
                    //                textureColorListNum = getNullsColorList(gl);
                    colorListChanged = false;
                }
                gl.glCallList(textureColorListNum);

            if (textureChanged)
            {
                byte[] textureData = currentSurfaceTexture.getTextureData();
                textureTiles.displayTiles(this, gl, isPixelMode, textureData, nullByte);
                textureData = null;
                textureChanged = false;
            }
            else
            {
                 textureTiles.displayTiles(this, gl, isPixelMode, (byte[]) null, nullByte);
            }

            if (shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);

            if (lighting)
            {
                gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL, GL.GL_SINGLE_COLOR);
            }
            else
            {
                gl.glEnable(GL.GL_LIGHTING);
            }
            if (transparency) gl.glDisable(GL.GL_POLYGON_STIPPLE);
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glDisable(GL.GL_BLEND);

            normals = null;
        }
        else if (currentModel.useDisplayLists)
        {
            if (surfDisplayListNum == 0) // build display list
            {
                surfDisplayListNum = gl.glGenLists(1);
                if (surfDisplayListNum == 0)
                {
                    StsMessageFiles.logMessage(
                        "System Error in StsGrid.displaySurface: " +
                            "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(surfDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
                drawSurfaceFillWithoutNulls(gl);
                gl.glEndList();

                /** We are using normals only once with displayList, so null it since we are finished with it. */
                normals = null;

                //timer.stop("display list surface setup: ");
            }
            gl.glCallList(surfDisplayListNum);
        }
        else
        {
            drawSurfaceFillWithoutNulls(gl);
        }
        //        gl.glEnable(GL.GL_LIGHTING);
    }

    public void doCropChanged()
    {
        if (textureTiles == null) return;
        //    addTextureToDeleteList(glPanel3d);
        textureTiles.cropChanged();
    }

    /*
       public void setActionTextureDisplayable(StsSurfaceTexture newTextureDisplayable)
       {
           actionTextureDisplayable = newTextureDisplayable;
           this.newTextureDisplayable = newTextureDisplayable;
       }
    */

    private void debugPrintTexture(String message, StsSurfaceTexture texture)
    {
        if (texture == null)
            System.out.println(message + ": null");
        else
            System.out.println(message + ": " + texture.getName());
    }

    private void checkSurfaceTexture()
    {
        if (newSurfaceTexture == null) return;
        if (newSurfaceTexture.isDisplayable()) return;
        constructTextureDisplayableList();
        setNewSurfaceTexture(surfaceTextureList[0]);
    }

    public StsSurfaceTexture getNewSurfaceTexture()
    {
        return newSurfaceTexture;
    }

    /** set the current surfaceTextureDisplayable, but not if its the already this one */
    public void setNewSurfaceTexture(StsSurfaceTexture newSurfaceTexture)
    {
        if (this.newSurfaceTexture == newSurfaceTexture) return;
        this.newSurfaceTexture = newSurfaceTexture;
        newSurfaceTexture.selected();
        StsColorscale colorscale = newSurfaceTexture.getColorscale();
        if (colorscaleBean != null)
        {
            colorscaleBean.setValueObject(colorscale);
            colorscaleBean.setHistogram(newSurfaceTexture.getHistogram());
        }
        float min = newSurfaceTexture.getDataMin();
        float max = newSurfaceTexture.getDataMax();
        float range = (max - min) * 5.0f;
        if (propertyMinBean != null)
        {
            propertyMinBean.setValueAndRange(min, min - range, max + range);
            propertyMaxBean.setValueAndRange(max, min - range, max + range);
        }
        setSurfaceTextureName(newSurfaceTexture.getName());
        textureChanged = true;
        colorListChanged = true;
        currentModel.viewObjectRepaint(this, this);
    }

    public float getPropertyMin() { return propertyMin; }

    public float getPropertyMax() { return propertyMax; }

    public void setDefaultSurfaceTexture(StsGLPanel glPanel)
    {
        setSurfaceTexture(surfaceDefaultTexture, glPanel);
    }

    /** set the current surfaceTextureDisplayable, but not if its the already this one */
    public void setSurfaceTexture(StsSurfaceTexture surfaceTexture, StsGLPanel glPanel)
    {
        if (currentSurfaceTexture == surfaceTexture) return;
//         if (glPanel != null && currentSurfaceTexture != null) StsTextureList.addTextureToDeleteList(this);
        setNewSurfaceTexture(surfaceTexture);
    }

    public void setSurfaceTextureName(String surfaceTextureName)
    {
        if (debug)
            StsException.systemDebug(this, "surfaceTextureName", " from: " + this.surfaceTextureName + " to: " + surfaceTextureName);
        if(this.surfaceTextureName == surfaceTextureName) return;
        this.surfaceTextureName = surfaceTextureName;
        dbFieldChanged("surfaceTextureName", surfaceTextureName);
    }

    public String getSurfaceTextureName()
    {
        return surfaceTextureName;
    }

    /*
       private void changeTextureDisplayable(String name)
       {
           if (textureDisplayableList==null)
               initTextureDisplayableList();
           for (int i = 0; i < textureDisplayableList.length; i++)
           {
               if (name.equals(textureDisplayableList[i].getName()))
                   setNewTextureDisplayable(textureDisplayableList[i]);
           }
       }
    */
    /** unconditionally set the surfaceTextureDisplayable to this one */
    /*
        public void changeTextureDisplayable(StsSurfaceTexture textureDisplayable)
        {
            this.newTextureDisplayable = textureDisplayable;
            textureChanged = true;
            currentModel.win3dDisplay();
        }
    */
    public void setTextureChanged()
    {
        textureChanged = true;
    }

    public void setTextureChanged(StsSurfaceTexture changedTexture)
    {
     //   if (changedTexture == currentSurfaceTexture)
        {
            textureChanged = true;
        }
    }

    private boolean initializeTextureTiles(StsGLPanel3d glPanel3d, GL gl)
    {
        StsCropVolume subVolume = currentModel.getProject().getCropVolume();
        //        if (textureTiles != null) return true;
        //        int shader = StsJOGLShader.NONE;
        //        if (contourColor && newTextureDisplayable != this)
        //           shader = StsJOGLShader.ARB_TLUT_WITH_SPECULAR_LIGHTS;
        if (textureTiles != null && textureTiles.isSameSize(this) && !shaderChanged()) return true;
        if (textureTiles != null) deleteTexturesAndDisplayLists(gl);
        textureTiles = StsTextureTiles.constructor(currentModel, this, StsParameters.ZDIR, this, isPixelMode, subVolume);
        if (textureTiles == null) return false;
        return true;
    }


    private boolean shaderChanged()
    {
        boolean usingShader = textureTiles.shader != StsJOGLShader.NONE;
        boolean useShader = contourColor && StsJOGLShader.canUseShader;
        if (useShader == usingShader) return false;
        textureChanged = true;
        if (useShader)
            shader = StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
        else
            shader = StsJOGLShader.NONE;
        return true;
    }

    /*
       public void deleteColorDisplayList(GL gl)
       {
           if(colorDisplayListNum == 0) return;
    if(newTextureDisplayable != null) newTextureDisplayable.deleteColorDisplayList(gl);
           else
               gl.glDeleteLists(colorDisplayListNum, 1 );
           colorDisplayListNum = 0;
       }
    */

    public byte[] getTextureData()
    {
        return surfaceDefaultTexture.getTextureData();
    }

    public int getColorDisplayListNum(GL gl, boolean nullsFilled)
    {
        if (colorDisplayListNum != 0)
        {
            gl.glDeleteLists(colorDisplayListNum, 1);
        }
        colorDisplayListNum = gl.glGenLists(1);
        if (colorDisplayListNum == 0)
        {
            StsException.systemError(this, "getColorDisplayListNum", "Failed to allocate a display list");
            return 0;
        }

        //        boolean nullsFilled = ((StsSurfaceClass)getCreateStsClass()).getNullsFilled();
        gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
        createColorList(gl, nullsFilled);
        gl.glEndList();

        return colorDisplayListNum;
    }

    private void createColorList(GL gl, boolean nullsFilled)
    {
        float[][] arrayRGBA = computeRGBAArray(nullsFilled);
        int nColors = arrayRGBA[0].length;
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
        gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
        //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
        //        arrayRGBA = null;
    }

    /*
         public void createColorTLUT(GL gl, boolean nullsFilled)
         {
             float[][] arrayRGBA = computeRGBAArray(nullsFilled);
             StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
             arrayRGBA = null;
         }
    */
    public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
    {
        float[][] arrayRGBA = computeRGBAArray(nullsFilled);
        return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
    }

    /**
     * null values (for now) will be drawn transparent,
     * so we only need to set the single non-null value
     */
    private float[][] computeRGBAArray(boolean nullsFilled)
    {
        float[] rgbaColor = getStsColor().getRGBA();
        float[][] arrayRGBA = new float[4][maxNullIndex + 1];

        // make NOT_NULL index the color of the surface
        // the shader can't handle this sort of ramp, so fill...
        for (int i = 0; i < 4; i++)
        {
            arrayRGBA[i][0] = rgbaColor[i];
            arrayRGBA[i][1] = rgbaColor[i];
        }

        // make all other nulls the null fill color
        if (nullsFilled)
        {
            String nullColorName = ((StsSurfaceClass) getCreateStsClass()).getNullColorName();
            StsColor nullColor = ((StsSurfaceClass) getCreateStsClass()).getNullColor(stsColor);
            if (nullColor == null)
            {
                return arrayRGBA;
            }

            float[] rgbaGray = nullColor.getRGBA();
            for (int i = 0; i < 4; i++)
            {
                for (int n = 2; n <= maxNullIndex; n++)
                {
                    arrayRGBA[i][n] = rgbaGray[i];
                }
            }
        }

        return arrayRGBA;
    }

    protected void drawSurfaceFillWithoutNulls(GL gl)
    {
        int i, j;

        if (timer == null)
        {
            timer = new StsTimer();
        }
        timer.start();

        getComputeNormals();

        float[] point = new float[3];

        float rowStartX = xMin;
        float rowStartY = yMin;

        //        StsColor.setGLColor(gl, this.stsColor);
        if (smoothNormals)
            gl.glShadeModel(GL.GL_SMOOTH);
        else
            gl.glShadeModel(GL.GL_FLAT);

        for (i = 0; i < nRows - 1; i++)
        {
            point[0] = xMin;
            point[1] = rowStartY;

            gl.glBegin(GL.GL_TRIANGLE_STRIP);

            for (j = 0; j < nCols; j++)
            {
                point[2] = pointsZ[i][j];
                if (lighting && smoothNormals)
                {
                    gl.glNormal3fv(normals[i][j], 0);
                }
                gl.glVertex3fv(point, 0);

                point[1] += (float) yInc;
                point[2] = pointsZ[i + 1][j];
                if (lighting)
                {
                    gl.glNormal3fv(normals[i + 1][j], 0);
                }
                gl.glVertex3fv(point, 0);

                point[0] += xInc;
                point[1] -= yInc;
            }

            gl.glEnd();

            rowStartY += yInc;
        }
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {

        float[][] points = getPoints();

        float[][][] normals = null;
        if (lighting)
        {
            normals = getComputeNormals();
        }
        if (hasCoarseGrid)
        {
            drawCoarseTextureTileSurface(tile, gl, points, normals);
        }
        else
        {
            drawFineTextureTileSurface(tile, gl, points, normals);
        }
    }

    public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
    {
    }

    public void drawCoarseTextureTileSurface(StsTextureTile tile, GL gl, float[][] points, float[][][] normals)
    {
        if (pointsZ == null)
        {
            return;
        }

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

        // if tile is narrower than stride, reduce stride to width of tile
        int tileRowStride = Math.min(rowStride, nRows - 1);
        int tileColStride = Math.min(colStride, nCols - 1);

        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        double dRowTexCoorStride = tileRowStride * dRowTexCoor;
        double dColTexCoorStride = tileColStride * dColTexCoor;

        float xIncStride = tileColStride * xInc;
        float yIncStride = tileRowStride * yInc;

        double rowTexCoor = tile.minRowTexCoor;
        double nextRowTexCoor = rowTexCoor + dRowTexCoorStride;

        int row = tile.croppedRowMin;
        int nextRow = row + tileRowStride;
        float xStart = xMin + tile.croppedColMin * xInc;
        float x = xStart;
        float y = yMin + tile.croppedRowMin * yInc;
        float nextY = y + yIncStride;

        if (smoothNormals)
            gl.glShadeModel(GL.GL_SMOOTH);
        else
            gl.glShadeModel(GL.GL_FLAT);

        while (true)
        {
            int col = tile.croppedColMin;
            double colTexCoor = tile.minColTexCoor;
            x = xStart;
            gl.glBegin(GL.GL_QUAD_STRIP);
            while (true)
            {
                float zBot = points[row][col];
                boolean botOk = zBot != nullValue;
                float zTop = points[nextRow][col];
                boolean topOk = zTop != nullValue;
                boolean ok = topOk || botOk;
                if(!botOk)
                    zBot = zTop;
                else if(!topOk)
                    zTop = zBot;

                if(ok)
                {
                    gl.glTexCoord2d(colTexCoor, rowTexCoor);
                    if (lighting && smoothNormals)
                    {
                        gl.glNormal3fv(normals[row][col], 0);
                    }
                    gl.glVertex3f(x, y, zBot);
                    gl.glTexCoord2d(colTexCoor, nextRowTexCoor);
                    if (lighting)
                    {
                        gl.glNormal3fv(normals[nextRow][col], 0);
                    }
                    gl.glVertex3f(x, nextY, zTop);
                }

                if (col == tile.croppedColMax)
                {
                    break;
                }

                int nextCol = col + tileColStride;
                if (nextCol <= tile.croppedColMax)
                {
                    col = nextCol;
                    x += xIncStride;
                    colTexCoor += dColTexCoorStride;
                }
                else // nextCol > tile.colMax
                {
                    int dCol = tile.croppedColMax - col;
                    col = tile.croppedColMax;
                    x += dCol * xInc;
                    colTexCoor += dCol * dColTexCoor;
                }
            }
            gl.glEnd();

            if (nextRow >= tile.croppedRowMax)
            {
                break;
            }

            row = nextRow;
            nextRow = row + tileRowStride;
            if (nextRow <= tile.croppedRowMax)
            {
                y = nextY;
                nextY += yIncStride;
                rowTexCoor = nextRowTexCoor;
                nextRowTexCoor += dRowTexCoorStride;
            }

            else // nextRow > tile.rowMax
            {
                int dRow = tile.croppedRowMax - row;
                nextRow = tile.croppedRowMax;
                y = nextY;
                nextY += dRow * yInc;
                rowTexCoor = nextRowTexCoor;
                nextRowTexCoor += dRow * dRowTexCoor;
            }
        }
    }


    public void drawFineTextureTileSurface(StsTextureTile tile, GL gl, float[][] points, float[][][] normals)
    {
        float xStart = xMin + tile.croppedColMin * xInc;
        float x = xStart;
        float y = yMin + tile.croppedRowMin * yInc;
        if (pointsZ == null)
        {
            return;
        }

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;

        if (smoothNormals)
            gl.glShadeModel(GL.GL_SMOOTH);
        else
            gl.glShadeModel(GL.GL_FLAT);


        for (int row = tile.croppedRowMin; row < tile.croppedRowMax; row++, y += yInc, rowTexCoor += dRowTexCoor)
        {
            double colTexCoor = tile.minColTexCoor;
            x = xStart;
            gl.glBegin(GL.GL_QUAD_STRIP);

            for (int col = tile.croppedColMin; col <= tile.croppedColMax; col++, x += xInc, colTexCoor += dColTexCoor)
            {
                float zBot = points[row][col];
                boolean botOk = zBot != nullValue;
                float zTop = points[row+1][col];
                boolean topOk = zTop != nullValue;
                if(!topOk && !botOk) continue;
                if(!botOk)
                    zBot = zTop;
                else if(!topOk)
                    zTop = zBot;

                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                if (lighting && smoothNormals)
                {
                    gl.glNormal3fv(normals[row][col], 0);
                }
                gl.glVertex3f(x, y, zBot);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                if (lighting)
                {
                    gl.glNormal3fv(normals[row + 1][col], 0);
                }
                gl.glVertex3f(x, y + yInc, zTop);
            }
            gl.glEnd();
        }
    }

    /*
    protected void drawSurfaceFillWithoutNullsGridCentered(GL gl, StsSurfaceDisplayable surfaceDisplayable)
    {
     int i, j;

     float[] point = new float[3];

    Color[][] colors = surfaceDisplayable.get2dColorArray(this, pointsZ, offset);

     float rowStartX = xMin;
     float rowStartY = yMin;

     for (i=0; i<nRows-1; i++)
     {
      point[0] = xMin;
      point[1] = rowStartY;

      gl.glBegin(GL.GL_TRIANGLE_STRIP);

      for (j=0; j<nCols; j++)
      {
       point[2] = pointsZ[i][j];
       StsColor.setGLJavaColor(gl, colors[i][j]);
       if(isLighted) gl.glNormal3fv(normals[i][j]);
       gl.glVertex3fv(point);

       point[1] += (float)yInc;
       point[2] = pointsZ[i+1][j];

       StsColor.setGLJavaColor(gl, colors[i+1][j]);
       if(isLighted) gl.glNormal3fv(normals[i+1][j]);
       gl.glVertex3fv(point);

       point[0] += xInc;
       point[1] -= yInc;
      }

      gl.glEnd();

      rowStartY += yInc;
     }
    }
    */
    public void checkConstructGridNormals()
    {
        computeNormals();
    }

    protected void computeNormals()
    {
        if (smoothNormals)
            normals = StsToolkit.computeSmoothNormals(pointsZ, nRows, nCols, xInc, yInc);
        else
            normals = StsToolkit.computeFaceNormals(pointsZ, nRows, nCols, xInc, yInc);
    }

    /** construct normals from tstrips */
/*
    public void computeTStripNormals() throws StsException
    {
        if (normals != null)
        {
            return;
        }
        normals = TriangleStrip.computeTStripNormals(this, tStrips, pointsZ);
    }
*/
    /*
    public void computeTStripNormals() throws StsException
    {
     int i, j, n;
     int n0, n1;

     if(normals != null) return;
           normals = new float[nRows][nCols][3];

     float[] idif = new float[3];
     float[] jdif = new float[3];
     float[] normal;

     idif[0] = idif[1] = (float)yInc;
     jdif[0] = (float)xInc;

     if (tStrips == null) makeTStrips();
     int nTStrips = tStrips.getSize();
     for (n = 0; n < nTStrips; n++)
     {
      TriangleStrip t = (TriangleStrip)tStrips.getElement(n);
      i = t.rowOrColNumber;
      j = t.firstIndex;

      // first side
      if (t.firstSide == STRIP_TOP)
      {
       idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
       jdif[2] = pointsZ[i+1][j] - pointsZ[i+1][j-1];
       normals[i+1][j-1] = StsGridPoint.leftCrossProduct(idif, jdif);
        }
      else if (t.firstSide == STRIP_BOT)
      {
       idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
       jdif[2] = pointsZ[i][j] - pointsZ[i][j-1];
       normals[i][j-1] = StsGridPoint.leftCrossProduct(idif, jdif);
      }

      // middle sides excluding last middle side
      for (; j<t.lastIndex; j++)
      {
       idif[2] = pointsZ[i+1][j] - pointsZ[i][j];

       jdif[2] = pointsZ[i+1][j+1] - pointsZ[i+1][j];
       normals[i+1][j] = StsGridPoint.leftCrossProduct(idif, jdif);

       jdif[2] = pointsZ[i][j+1] - pointsZ[i][j];
       normals[i][j] = StsGridPoint.leftCrossProduct(idif, jdif);
      }
      // j has been bumped to next j by for loop above
      // last side: last idif is from bottom row if available; otherwise top row

      idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
      normal = StsGridPoint.leftCrossProduct(idif, jdif);
      normals[i][j] = normal;
      normals[i+1][j] = normal;

      if (t.lastSide == STRIP_BOT)
      {
       jdif[2] = pointsZ[i][j+1] -  pointsZ[i][j];
       normals[i][j+1] = StsGridPoint.leftCrossProduct(idif, jdif);
      }
      else if (t.lastSide == STRIP_TOP)
      {
       jdif[2] = pointsZ[i+1][j+1] - pointsZ[i+1][j];
       normals[i+1][j+1] = StsGridPoint.leftCrossProduct(idif, jdif);
      }
     }
    }
    */
    //	public void deleteDisplayLists() { deleteDisplayLists = true; }


    protected void checkDeleteDisplays(GL gl)
    {
        if (geometryChanged || textureChanged) return;
        //        surfaceChanged = true;
        //        if(surfDisplayListNum == 0 && gridDisplayListNum == 0) return;
        deleteDisplayLists(gl);
        deleteTexturesAndDisplayLists(gl);
    }

    public void deleteDisplayLists(GL gl)
    {
        if (textureTiles != null)
        {
            textureTiles.deleteDisplayLists(gl);
        }

        if (surfDisplayListNum > 0)
        {
            gl.glDeleteLists(surfDisplayListNum, 1);
            surfDisplayListNum = 0;
        }
        if (gridDisplayListNum > 0)
        {
            gl.glDeleteLists(gridDisplayListNum, 1);
            gridDisplayListNum = 0;
        }
        geometryChanged = true;
        //        deleteDisplayLists = false;
    }

    private void displayGridLines(StsGLPanel3d glPanel3d)
    {
        //timer.start();

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        if (gridDisplayListNum == 0 && currentModel.useDisplayLists) // build display list
        {
            gridDisplayListNum = gl.glGenLists(1);
            if (gridDisplayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsGrid.displayGrid: " +
                    "Failed to allocate a display list");
                return;
            }

            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glNewList(gridDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
            drawGridLines(glPanel3d, gl);
            gl.glEndList();
            glPanel3d.resetViewShift(gl);

            //timer.stop("display list surface setup: ");

        }
        else if (currentModel.useDisplayLists) // use existing display list
        {
            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glCallList(gridDisplayListNum);
            glPanel3d.resetViewShift(gl);
            //timer.stop("display list surface draw: ");
        }
        else // immediate mode draw
        {
            if (gridDisplayListNum > 0)
            {
                gl.glDeleteLists(gridDisplayListNum, 1);
                gridDisplayListNum = 0;
            }
            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            drawGridLines(glPanel3d, gl);
            glPanel3d.resetViewShift(gl);
        }
    }

    private void drawGridLines(StsGLPanel3d glPanel3d, GL gl)
    {
        gl.glDisable(GL.GL_LIGHTING);
        //        win3d.glPanel3d.setViewShift(StsGraphicParameters.gridShift);

        gl.glLineWidth(StsGraphicParameters.gridLineWidth);

        float[][] points = getPoints();
        drawILines(gl, points);
        drawJLines(gl, points);

        gl.glEnable(GL.GL_LIGHTING);
        //        win3d.glPanel3d.resetViewShift();
    }

    private void drawILines(GL gl, float[][] points)
    {
        int i, j, n;
        boolean notNull;
        boolean drawing = false;

        float[] point = new float[3];

        try
        {
            float rowStartX = xMin;
            float rowStartY = yMin;

            for (i = 0, n = 0; i < nRows; i++)
            {
                point[0] = rowStartX;
                point[1] = rowStartY;

                for (j = 0; j < nCols; j++, n++)
                {
                    notNull = isPointNotNullOrFilled(i, j);
                    if (!drawing)
                    {
                        if (notNull)
                        {
                            drawing = true;
                            gl.glBegin(GL.GL_LINE_STRIP);
                            point[2] = points[i][j];
                            gl.glVertex3fv(point, 0);
                        }
                    }
                    else
                    {
                        if (!notNull)
                        {
                            gl.glEnd();
                            drawing = false;
                        }
                        else
                        {
                            point[2] = points[i][j];
                            gl.glVertex3fv(point, 0);
                        }
                    }
                    point[0] += (float) xInc;
                }

                if (drawing)
                {
                    gl.glEnd();
                    drawing = false;
                }
                rowStartY += yInc;
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.drawILines() failed.", e,
                StsException.WARNING);
        }
    }

    public boolean isPointZOK(int row, int col)
    {
        if (pointsZ == null) return false;
        return pointsZ[row][col] != nullValue;
    }

    public boolean isPointNotNullOrFilled(int row, int col)
    {
        if (pointsNull == null) return true;
        byte pointNull = pointsNull[row][col];
        return pointNull == SURF_PNT || pointNull == SURF_GAP;
    }

    private void drawJLines(GL gl, float[][] points)
    {
        boolean notNull;
        boolean drawing = false;
        int i, j, n;

        float[] point = new float[3];

        try
        {
            float colStartX = xMin;
            float colStartY = yMin;

            for (j = 0; j < nCols; j++)
            {
                point[0] = colStartX;
                point[1] = colStartY;

                for (i = 0; i < nRows; i++)
                {
                    n = getIndex(i, j);
                    notNull = isPointNotNullOrFilled(i, j);

                    if (!drawing)
                    {
                        if (notNull)
                        {
                            drawing = true;
                            gl.glBegin(GL.GL_LINE_STRIP);
                            point[2] = points[i][j];
                            gl.glVertex3fv(point, 0);
                        }
                    }
                    else
                    {
                        if (!notNull)
                        {
                            gl.glEnd();
                            drawing = false;
                        }
                        else
                        {
                            point[2] = points[i][j];
                            gl.glVertex3fv(point, 0);
                        }
                    }
                    point[1] += (float) yInc;
                }

                if (drawing)
                {
                    gl.glEnd();
                    drawing = false;
                }
                colStartX += xInc;
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGrid.drawJLines() failed.", e,
                StsException.WARNING);
        }
    }

    /*
       private void computeTStrips()
       {
     if(pointsZ == null)
     {
      StsException.outputException(new StsException(StsException.WARNING,
       "StsGrid.makeTStrips() failed. StsGrid.pointsZ are null."));
      return;
     }

           tStrips = TriangleStrip.computeTStrips(this, pointsNull);
       }
    */
    /*
    private void makeTStrips() throws StsException
    {
     int row, col;
     int colMin, colMax, nextColMin, nextColMax;
     int nBot, nTop;
     byte pointNull;
     boolean botNotNull, topNotNull;
     TriangleStrip t = null;
     boolean startedTStrip = false;
     int nPoints = 0;

     if(pointsZ == null)
     {
      StsException.outputException(new StsException(StsException.WARNING,
       "StsGrid.makeTStrips() failed. StsGrid.pointsZ are null."));
      return;
     }

     StsGridIterator gridIterator = new StsGridIterator(this);
     int nGridRows = gridIterator.getNRows();
     tStrips = new StsList(nGridRows, nGridRows);

     while((row = gridIterator.getNextRow()) != nullInteger)
     {
      colMin = gridIterator.colMin;
      colMax = gridIterator.colMax;

      nextColMin = gridIterator.nextColMin;
      nextColMax = gridIterator.nextColMax;

      int colStart = Math.min(colMin, nextColMin);
      int colEnd = Math.max(colMax, nextColMax);

      for (col=colStart; col <= colEnd; col++)
      {
       pointNull = pointsNull[row][col];
       botNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED || pointNull == NULL_GAP || pointNull == NULL_GAP_OR_BOUNDARY);
//				botNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED);

       pointNull = pointsNull[row+1][col];
       topNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED || pointNull == NULL_GAP || pointNull == NULL_GAP_OR_BOUNDARY);
//				topNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED);

       if (!startedTStrip)
       {
        nPoints = 0;
        if (botNotNull && topNotNull)
        {
         t = new TriangleStrip();  // create a new tstrip
         t.firstSide = STRIP_BOTH;
         t.firstIndex = col;
         t.lastIndex = col;
         nPoints += 2;
         startedTStrip = true;
        }
        else if (botNotNull) // top point null
        {
         t = new TriangleStrip();  // create a new tstrip
         t.firstSide = STRIP_BOT;
         t.firstIndex = col+1;
         nPoints++;
         startedTStrip = true;
        }
        else if (topNotNull) // bot point null
        {
         t = new TriangleStrip();  // create a new tstrip
         t.firstSide = STRIP_TOP;
         t.firstIndex = col+1;
         nPoints++;
         startedTStrip = true;
        }
       }  // !startedTStrip

       else  // tstrip started
       {
        if (!botNotNull && !topNotNull)  // both null - end tstrip
        {
         if (nPoints > 2)  // need 3 or more points
         {
          // add the triangle strip
          t.rowOrColNumber = row;
          tStrips.add(t);
         }
         startedTStrip = false;  // reset flag
        }
        else if (botNotNull && topNotNull)  // continue tstrip
        {
         t.lastIndex = col;
         t.lastSide = STRIP_BOTH;
         nPoints += 2;
        }
        else if (botNotNull) // top point null  - end this tstrip
        {
         nPoints++;
         if (nPoints > 2)
         {
          t.lastSide = STRIP_BOT;
          t.rowOrColNumber = row;
          tStrips.add(t);
         }
         t = new TriangleStrip();  // start a new tstrip
         t.firstSide = STRIP_BOT;
         t.firstIndex = col+1;
         nPoints = 1;
        }
        else // (topNotNull) top point not null - end this tstrip
        {
         nPoints++;
         if (nPoints > 2)
         {
          t.lastSide = STRIP_TOP;
          t.rowOrColNumber = row;
          tStrips.add(t);
         }
         t = new TriangleStrip();  // start a new tstrip
         t.firstSide = STRIP_TOP;
         t.firstIndex = col+1;
         nPoints = 1;
        }
       }  // tstrip started
      }  // for "col" loop

      // reached end of a row
      if (startedTStrip && nPoints>2)
      {
       t.rowOrColNumber = row;
       tStrips.add(t);
      }
      startedTStrip = false;  // reset flag
     }  // for "i" loop

     tStrips.trimToSize();  // shrink to actual allocation
    }
    */
    /* display surface containing null Z values by drawing triangle strips */
/*
    protected void drawSurfaceFillWithNulls(GL gl) throws StsException
    {
        TriangleStrip.drawSurfaceFillWithNulls(gl, this, tStrips, pointsZ,
            normals, lighting);
    }
*/
    /*
    protected void drawSurfaceFillWithNulls(GL gl) throws StsException
    {
     TriangleStrip t = null;
     float[] point = new float[3];
     int i, j = 0;
     int n;

     float startX = xMin;
     float startY = yMin;

     try
     {
//			StsColor.setGLColor(gl, this.stsColor);

      int nTStrips = tStrips.getSize();
      for (int nt = 0; nt < nTStrips; nt++)
      {
       t = (TriangleStrip)tStrips.getElement(nt);
       i = t.rowOrColNumber;
       j = t.firstIndex;

       gl.glBegin(GL.GL_TRIANGLE_STRIP);

       // first side
       if (t.firstSide == STRIP_BOT)
       {
        point[0] = (float)(startX + (j-1)*xInc);
        point[1] = (float)(startY + i*yInc);
        point[2] = pointsZ[i][j-1];

        if(isLighted) gl.glNormal3fv(normals[i][j-1]);
        gl.glVertex3fv(point);
        if(isLighted) gl.glNormal3fv(normals[i][j-1]);
        gl.glVertex3fv(point);

        point[0] += xInc;
       }
       else if (t.firstSide == STRIP_TOP)
       {
        point[0] = (float)(startX + (j-1)*xInc);
        point[1] = (float)(startY + (i+1)*yInc);
        point[2] = pointsZ[i+1][j-1];

        if(isLighted) gl.glNormal3fv(normals[i+1][j-1]);
        gl.glVertex3fv(point);
        if(isLighted) gl.glNormal3fv(normals[i+1][j-1]);
        gl.glVertex3fv(point);

        point[0] += xInc;
        point[1] -= yInc;
       }
       else
       {
        point[0] = startX + j*xInc;
        point[1] = startY + i*yInc;
       }
       // middle sides

       for (; j<=t.lastIndex; j++)
       {
        point[2] = pointsZ[i][j];

        if(isLighted) gl.glNormal3fv(normals[i][j]);
        gl.glVertex3fv(point);

        point[1] += yInc;
        point[2] = pointsZ[i+1][j];

        if(isLighted) gl.glNormal3fv(normals[i+1][j]);
        gl.glVertex3fv(point);

        point[0] += xInc;
        point[1] -= yInc;
       }

       // last side: j has been bumped to next j in for-loop above
       if (t.lastSide == STRIP_BOT)
       {
        point[2] = pointsZ[i][j];

        if(isLighted) gl.glNormal3fv(normals[i][j]);
        gl.glVertex3fv(point);
       }
       else if (t.lastSide == STRIP_TOP)
       {
        point[1] += yInc;
        point[2] = pointsZ[i+1][j];

        if(isLighted) gl.glNormal3fv(normals[i+1][j]);
        gl.glVertex3fv(point);
       }

       gl.glEnd();
      }
     }
     catch(Exception e)
     {
      StsException.outputException("Exception in display. " +
       "row: " + t.rowOrColNumber + " col: " + j, e, StsException.WARNING);
     }
    }

    protected void drawSurfaceFillWithNullsGridCentered(GL gl, StsSurfaceDisplayable surfaceDisplayable)
    {
     TriangleStrip t = null;
     float[] point = new float[3];
     int i, j = 0;
     int n;

     try
     {

      float startX = xMin;
      float startY = yMin;

      int nTStrips = tStrips.getSize();
      for (int nt = 0; nt < nTStrips; nt++)
      {
       t = (TriangleStrip)tStrips.getElement(nt);
       i = t.rowOrColNumber;
       j = t.firstIndex;

       gl.glBegin(GL.GL_TRIANGLE_STRIP);

       // first side
       if (t.firstSide == STRIP_BOT)
       {
        point[0] = (float)(startX + (j-1)*xInc);
        point[1] = (float)(startY + i*yInc);
        point[2] = pointsZ[i][j-1];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i, j-1, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i][j-1]);
        gl.glVertex3fv(point);
        if(isLighted) gl.glNormal3fv(normals[i][j-1]);
        gl.glVertex3fv(point);

        point[0] += xInc;
       }
       else if (t.firstSide == STRIP_TOP)
       {
        point[0] = (float)(startX + (j-1)*xInc);
        point[1] = (float)(startY + (i+1)*yInc);
        point[2] = pointsZ[i+1][j-1];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i+1, j-1, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i+1][j-1]);
        gl.glVertex3fv(point);
        if(isLighted) gl.glNormal3fv(normals[i+1][j-1]);
        gl.glVertex3fv(point);

        point[0] += xInc;
        point[1] -= yInc;
       }
       else
       {
        point[0] = startX + j*xInc;
        point[1] = startY + i*yInc;
       }
       // middle sides

       for (; j<=t.lastIndex; j++)
       {
        point[2] = pointsZ[i][j];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i, j, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i][j]);
        gl.glVertex3fv(point);

        point[1] += yInc;
        point[2] = pointsZ[i+1][j];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i+1, j, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i+1][j]);
        gl.glVertex3fv(point);

        point[0] += xInc;
        point[1] -= yInc;
       }

       // last side: j has been bumped to next j in for-loop above
       if (t.lastSide == STRIP_BOT)
       {
        point[2] = pointsZ[i][j];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i, j, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i][j]);
        gl.glVertex3fv(point);
       }
       else if (t.lastSide == STRIP_TOP)
       {
        point[1] += yInc;
        point[2] = pointsZ[i+1][j];

        StsColor.setGLJavaColor(gl, surfaceDisplayable.getGridColor(i+1, j, point[2] + offset));
        if(isLighted) gl.glNormal3fv(normals[i+1][j]);
        gl.glVertex3fv(point);
       }

       gl.glEnd();
      }
     }
     catch(Exception e)
     {
      StsException.outputException("Exception in display. " +
       "row: " + t.rowOrColNumber + " col: " + j, e, StsException.WARNING);
     }
    }
    */
    public boolean checkIsLoaded()
    {
        return pointsLoaded(type);
    }

    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        super.setIsVisible(b);

        if (b)
        {
            isVisible = pointsLoaded(type);
        }
        else
        {
            isVisible = b;
        }
        currentModel.viewObjectRepaint(this, this);
    }

    /** add fault polygons */
    public void addFaultPolygons(StsFaultPolygon[] polygons)
    {
        if (polygons == null)
        {
            return;
        }
        if (faultPolygons == null)
        {
            faultPolygons = polygons;
        }
        else
        {
            int nOld = faultPolygons.length;
            int nNew = polygons.length;
            StsFaultPolygon[] oldPolygons = faultPolygons;
            faultPolygons = new StsFaultPolygon[nOld + nNew];
            for (int i = 0; i < nOld; i++)
            {
                faultPolygons[i] = oldPolygons[i];
            }
            for (int i = 0; i < nNew; i++)
            {
                faultPolygons[i + nOld] = polygons[i];
            }
        }

        applyFaultPolygons();
    }

    /** apply fault polygons to grid */
    public void applyFaultPolygons()
    {
        if (!isType(IMPORTED) || faultPolygons == null)
        {
            return;
        }
        applyFaultPolygons(faultPolygons);
    }

    /** remove fault polygons from grid */
    public void removeFaultPolygons()
    {
        if (!isType(IMPORTED) || faultPolygons == null)
        {
            return;
        }
        removeFaultPolygons(faultPolygons);
        faultPolygons = null;
    }

    /** get fault polygon array */
    public StsFaultPolygon[] getFaultPolygons()
    {
        return faultPolygons;
    }

    public void setTransparent(boolean value)
    {
        transparency = value;
        dbFieldChanged("transparency", transparency);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getTransparent() { return transparency; }

    public boolean isType(byte type)
    {
        switch (type)
        {
            case IMPORTED:
                return this.type == IMPORTED;
            case MODEL:
                return this.type == MODEL ||
                    this.type == CHECKED ||
                    this.type == INITIALIZED ||
                    this.type == GAPPED ||
                    this.type == BUILT;
            case CHECKED:
                return this.type == CHECKED ||
                    this.type == INITIALIZED ||
                    this.type == GAPPED;
            case INITIALIZED:
                return this.type == INITIALIZED ||
                    this.type == GAPPED ||
                    this.type == BUILT;
            case GAPPED:
                return this.type == GAPPED ||
                    this.type == BUILT;
            case BUILT:
                return this.type == BUILT;
            default:
            {
                StsException.outputException(new StsException(StsException.
                    WARNING,
                    "StsSurface.isType() failed.", "Undefined type: " + type));
                return false;
            }
        }
    }

    public boolean initialize(StsModel model)
    {
        if (initialized)
        {
            return true;
        }

        try
        {
            //            boolean horizonsBuilt = currentModel.getBooleanStateProperty("horizonsBuilt");
            //            if(type == IMPORTED && horizonsBuilt) isVisible = false;
            //            else                                  isVisible = true;

            setup();
            setIsVisible(isVisible);
            if (faultPolygons != null)
            {
                for (int i = 0; i < faultPolygons.length; i++)
                {
                    faultPolygons[i].initialize();
                }
                applyFaultPolygons();
            }
            setOffsetRange();
            initializeSurfaceTextures();
            initialized = true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "initialize(StsModel)", "Surface: " + getName() + ". ", e);
        }
        return true;
    }

    protected void initializeSurfaceTextures()
    {
        initSurfaceTextureList();
        // Restore the persisted texture. We only persist the name so we need to rebuild the texture.
        textureChanged = true;
        colorListChanged = true;
        if (surfaceTextureName != null)
        {
            for (int i = 0; i < surfaceTextureList.length; i++)
            {
                if (surfaceTextureName.equals(surfaceTextureList[i].getName()))
                {
                    setNewSurfaceTexture(surfaceTextureList[i]);
                    return;
                }
            }
        }
        setNewSurfaceTexture(surfaceTextureList[0]);
    }

    /** find the surface point at a mouse cursor location */
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d)
    {
        return getSurfacePosition(mouse, display, true, glPanel3d);
    }

    // if display, display position in statusArea; if outsideGridOK, allow for a pick
    // off the grid but find the nearest grid boundary point
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, boolean outsideGridOK, StsGLPanel3d glPanel3d)
    {
        StsGridPoint gridPoint;

        StsMousePoint mousePoint = mouse.getMousePoint();

        // get projection points and check for valid z range
        double[] point0 = glPanel3d.getNearProjectedPoint(mousePoint);
        double[] point1 = glPanel3d.getFarProjectedPoint(mousePoint);
        if (Math.abs(point0[2] - point1[2]) < StsParameters.doubleRoundOff)
        {
            return null;
        }

        StsGridDefinition gridDefinition = currentModel.getGridDefinition();

        gridPoint = gridDefinition.getSurfacePosition(point0, point1, this, outsideGridOK);
        if (gridPoint != null && display)
        {
            displaySurfacePosition(gridPoint);
        }
        return gridPoint;
    }

    public void debugSurfacePoint(StsGridPoint gridPoint)
    {
        if (gridPoint == null)
        {
            return;
        }
        StsXYGridable grid = gridPoint.grid;
        if (grid instanceof StsBlockGrid)
        {
            ((StsBlockGrid) grid).debugSurfacePoint(gridPoint);
        }
    }

    public void display(StsGLPanel glPanel)
    {
        if (pointsZ == null)
        {
            return;
        }

        GL gl = glPanel.getGL();

        if (!isVisible)
        {
            checkDeleteDisplays(gl);
            return;
        }

        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
        {
            checkDeleteDisplays(gl);
            return;
        }
        // in displaySurfaceFill we now set textureChanged to true if colorListChanged is true
        // so we don't need to do this operation here as it is a texture operation
        /*
            if(colorListChanged)
            {
                clearTextureTileSurface(glPanel3d);
            }
        */
        byte projectZDomain = currentModel.getProject().getZDomain();
        if (projectZDomain != zDomainDisplayed)
        {
            checkDeleteDisplays(gl);
            zDomainDisplayed = projectZDomain;
        }

        if (!currentModel.useDisplayLists && usingDisplayLists)
        {
            deleteDisplayLists(gl);
            usingDisplayLists = false;
        }

        if (isType(GAPPED))
        {
            displaySurface((StsGLPanel3d) glPanel);
        }
        else
        {
            displaySurface((StsGLPanel3d) glPanel);
        }
    }

    /** Draw any map edges on section */
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, StsPoint[] points)
    {
        if (!currentModel.getProject().supportsZDomain(this.zDomainOriginal))
        {
            return;
        }
        boolean drawGaps = true;
        //       drawOnCursor3d(glPanel3d, points[0], points[3], drawGaps, false);
        drawOnCursor3d(glPanel3d, points[0], points[3], drawGaps, true);
    }

    /** Draw any map edges on section */
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo,
                               float dirCoordinate)
    {
        if (!currentModel.getProject().supportsZDomain(this.zDomainOriginal)) return;
        boolean drawGaps = true;
        drawOnCursor2d(glPanel3d, dirNo, dirCoordinate, drawGaps);
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d,
                                StsSeismicCurtain seismicCurtain)
    {
        if (!currentModel.getProject().supportsZDomain(this.zDomainOriginal)) return;
        boolean drawGaps = true;
        boolean drawDotted = true;
        drawOn3dCurtain(glPanel3d, seismicCurtain, drawGaps, drawDotted);
    }

    /** get center z value to use in ordering surfaces */
    public float getOrderingValue()
    {
        try
        {
            if (!checkIsLoaded())
            {
                return StsParameters.nullValue;
            }
            StsProject project = currentModel.getProject();
            return getInterpolatedZ(project.getXCenter(), project.getYCenter());
        }
        catch (Exception e)
        {
            return StsParameters.nullValue;
        }
    }

    static public Object[] sortSurfaces(StsObject[] surfaces)
    {
        try
        {
            if (surfaces == null)
            {
                return null;
            }
            int nSurfaces = surfaces.length;
            if (nSurfaces == 0)
            {
                return null;
            }

            ArrayList objects = new ArrayList(nSurfaces);
            for (int n = 0; n < nSurfaces; n++)
            {
                objects.add(surfaces[n]);
            }

            Comparator comparator = new SurfaceComparator();
            Collections.sort(objects, comparator);
            Object[] sortedSurfaces = (Object[]) Array.newInstance(surfaces[0].getClass(), nSurfaces);
            Iterator iter = objects.iterator();
            int n = 0;
            while (iter.hasNext())
            {
                sortedSurfaces[n++] = (StsObject) iter.next();
            }
            return sortedSurfaces;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurface.sortSurfaces() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    static public StsSurface[] getTimeSurfaces(StsSurface[] surfaces)
    {
        StsSurface[] timeSurfaces = (StsSurface[])StsMath.arraycopy(surfaces);
        int nSurfaces = surfaces.length;
        int nTimeSurfaces = 0;
        for(int n = 0; n < nSurfaces; n++)
        {
            if(StsProject.supportsTime(surfaces[n].getZDomainOriginal()))
                timeSurfaces[nTimeSurfaces++] = surfaces[n];
        }
        timeSurfaces = (StsModelSurface[])StsMath.trimArray(timeSurfaces, nTimeSurfaces);
        return timeSurfaces;
    }

    static public StsSurface[] getZDomainSurfaces(StsSurface[] surfaces, byte currentZDomain)
    {
        StsSurface[] domainSurfaces = (StsSurface[])StsMath.arraycopy(surfaces);
        int nSurfaces = surfaces.length;
        int nDomainSurfaces = 0;
        for(int n = 0; n < nSurfaces; n++)
        {
            byte surfaceZDomain = surfaces[n].zDomainSupported;
            if(surfaceZDomain == StsProject.TD_TIME_DEPTH || surfaceZDomain == currentZDomain)
                domainSurfaces[nDomainSurfaces++] = surfaces[n];
        }
        domainSurfaces = (StsSurface[])StsMath.trimArray(domainSurfaces, nDomainSurfaces);
        return domainSurfaces;
    }
    
    static StsObject[] getCurrentSurfaces()
    {
        StsObject[] surfaces;
        StsModel model = StsObject.getCurrentModel();
        surfaces = model.getObjectList(StsModelSurface.class);
        if (surfaces.length > 0)
        {
            return surfaces;
        }
        return model.getObjectList(StsSurface.class);
    }

    static public final class SurfaceComparator implements Comparator
    {
        SurfaceComparator()
        {
        }

        public int compare(Object o1, Object o2)
        {
            StsSurface s1 = (StsSurface) o1;
            StsSurface s2 = (StsSurface) o2;

            StsRotatedGridBoundingBox boundingBox = new StsRotatedGridBoundingBox(s1, false);
            boundingBox.originSet = true;
            boundingBox.initializedXY = true;
            boundingBox.intersectBoundingBox(s2);
            boolean isCongruent = s1.congruentWith(s2);
            Iterator quadCellIterator = boundingBox.getQuadCellIterator();
            // get 100 good comparisons between surfaces to determine order
            int iter = 0;
            int aboveSum = 0;
            int belowSum = 0;
            while (quadCellIterator.hasNext() && iter < 100)
            {
                iter++;
                int[] rowCol = (int[]) quadCellIterator.next();
                int row = rowCol[0];
                int col = rowCol[1];
                float x = boundingBox.getXCoor(col);
                float y = boundingBox.getYCoor(row);
                rowCol = s1.getNearestRowColCoors(x, y);
                int row1 = rowCol[0];
                int col1 = rowCol[1];
                rowCol = s2.getNearestRowColCoors(x, y);
                int row2 = rowCol[0];
                int col2 = rowCol[1];
                if (s1.getPointNull(row1, col1) != SURF_PNT || s2.getPointNull(row2, col2) != SURF_PNT) continue;
                float z1 = s1.getPointZ(row1, col1);
                float z2 = s2.getPointZ(row2, col2);
                if (z1 > z2) belowSum++;
                else if (z2 > z1) aboveSum++;
            }
            int totalSum = aboveSum + belowSum;
            if (aboveSum >= belowSum)
            {
                if (aboveSum < 5 * belowSum)
                    StsMessage.printMessage("Surfaces pass thru each other and ordering is not obvious for: " + s1.getName() + " and " + s2.getName() +
                        "\nUse updown buttons to set proper order.");
                return -1;
            }
            else
            {
                if (belowSum < 5 * aboveSum)
                    StsMessage.printMessage("Surfaces pass thru each other and ordering is not obvious for: " + s1.getName() + " and " + s2.getName() +
                        "\nUse updown buttons to set proper order.");
                return 1;
            }
        }
    }

    static public boolean surfacesGridsAreSame(StsSurface[] surfaces)
    {
        if (surfaces == null)
        {
            return false;
        }
        int nSurfaces = surfaces.length;
        if (nSurfaces == 0)
        {
            return false;
        }

        StsSurface surface = surfaces[0];
        for (int i = 1; i < nSurfaces; i++)
        {
            if (!surface.xyGridSameAs(surfaces[i]))
            {
                return false;
            }
        }
        return true;
    }

    public void cropChanged()
    {
        geometryChanged = true;
        setTextureChanged();
    }

    public String getLabel()
    {
        return "surface: " + getName();
    }
    
    public String toString() { return "surface: " + getName(); }

    public String getSurfacePointString(int row, int col)
    {
        return "   row: " + row + " col: " + col + " null: " + getNullName(row, col);
    }

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            surfaceTextureListBean = new StsComboBoxFieldBean(StsSurface.class, "newSurfaceTexture", "Property", "surfaceTextureList");
            propertyMinBean = new StsFloatFieldBean(StsSurface.class, "propertyMin", false, "Property Min:");
            propertyMaxBean = new StsFloatFieldBean(StsSurface.class, "propertyMax", false, "Property Max:");

            offsetBean = new StsFloatFieldBean(StsSurface.class, "offset", true, "Offset", true);
            colorscaleBean = new StsEditableColorscaleFieldBean(StsSurface.class, "colorscale");

            displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(StsSurface.class, "isVisible", "Enable"),
                new StsBooleanFieldBean(this, "transparent", "Make Transparent"),
                bornField, deathField,
                new StsBooleanFieldBean(StsSurface.class, "displayGrid", "Grid"),
                new StsBooleanFieldBean(StsSurface.class, "displayFill", "Fill"),
                new StsBooleanFieldBean(StsSurface.class, "isPixelMode", "Pixel Mode"),
                new StsBooleanFieldBean(StsSurface.class, "lighting", "Lighting"),
                new StsBooleanFieldBean(StsSurface.class, "smoothNormals", "Smooth lighting"),
                new StsBooleanFieldBean(StsSurface.class, "contourColor", "Contoured overlay"),
                //    	new StsBooleanFieldBean(StsSurface.class, "displayProperty", "Property"),
                offsetBean,
                new StsColorComboBoxFieldBean(StsSurface.class, "stsColor", "Color", currentModel.getSpectrum("Basic")),
                surfaceTextureListBean,
                propertyMinBean,
                propertyMaxBean,
                colorscaleBean
            };
        }
        /** When the comboBox is selected, we dynamically rebuild the list since it might have been changed. */
        //        textureDisplayableListBean.setPopupMenuListener(this);
        //        initSurfaceTextureDisplayableList();
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if (propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
                {
                    new StsStringFieldBean(StsSurface.class, "zDomainString", false, "Z Domain"),
                    new StsIntFieldBean(StsSurface.class, "nCols", false, "nX Cols"),
                    new StsIntFieldBean(StsSurface.class, "nRows", false, "nY Rows"),
                    new StsDoubleFieldBean(StsSurface.class, "xOrigin", false, "X Origin"),
                    new StsDoubleFieldBean(StsSurface.class, "yOrigin", false, "Y Origin"),
                    new StsFloatFieldBean(StsSurface.class, "xInc", false, "X Inc"),
                    new StsFloatFieldBean(StsSurface.class, "yInc", false, "Y Inc"),
                    new StsFloatFieldBean(StsSurface.class, "xMin", false, "X Loc Min"),
                    new StsFloatFieldBean(StsSurface.class, "yMin", false, "Y Loc Min"),
                    new StsFloatFieldBean(StsSurface.class, "xSize", false, "X Size"),
                    new StsFloatFieldBean(StsSurface.class, "ySize", false, "Y Size"),
                    new StsFloatFieldBean(StsSurface.class, "angle", false, "XY Rot Angle"),
                    new StsFloatFieldBean(StsSurface.class, "zMin", false, "Min Depth"),
                    new StsFloatFieldBean(StsSurface.class, "zMax", false, "Max Depth"),
                    new StsStringFieldBean(StsSurface.class, "nativeHorizontalUnitsString", false, "Native Horizontal Units:"),
                    new StsStringFieldBean(StsSurface.class, "nativeVerticalUnitsString", false, "Native Vertical Units:"),
                };
        }
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        setOffsetRange();
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        //        initSurfaceTextureDisplayableList();
        //        textureDisplayableListBean.setSelectedItem(newTextureDisplayable);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        getSurfaceClass().selected(this);
    }

    public StsSurfaceClass getSurfaceClass()
    {
        return (StsSurfaceClass) currentModel.getCreateStsClass(this);
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public boolean canExport() { return true; }

    public boolean export()
    {
        progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Surface export", false);
        progressBarDialog.pack();
        progressBarDialog.setVisible(true);
        String exportType = " in time.";
        if (isDepth) exportType = " in depth.";
        progressBarDialog.setLabelText("Exporting surface with name " + getName() + exportType);

        Runnable runExport = new Runnable()
        {
            public void run()
            {
                exportSurface();
            }
        };

        Thread exportThread = new Thread(runExport);

        exportThread.start();
        return true;
    }

    private void exportSurface()
    {
        PrintWriter printWriter = null;
        //        StsCursor cursor = new StsCursor(currentModel.win3d.glPanel3d, Cursor.WAIT_CURSOR);
        int maxProgress = nRows - 1;
        progressBarDialog.setProgressMax(maxProgress);
        try
        {
            // fill in all nulls
            interpolateDistanceTransform();

            // output in grid-seismic.txt format (x, y, t or z, line, xline) if seismic line-xline numbering is available
            // otherwise output in grid.txt format (t or z only)
            StsProject project = currentModel.getProject();
            StsRotatedGridBoundingBox rotatedBoundingBox = project.getRotatedBoundingBox();

            if (!angleSet) setAngle();
            float[][] exportZ;
            if (isDepth && this.adjPointsZ != null)
                exportZ = adjPointsZ;
            else
                exportZ = pointsZ;

            String exportType = "-time";
            if (isDepth) exportType = "-depth";

            float rowNumStart, colNumStart, rowNumInc, colNumInc;
            if (rotatedBoundingBox.rowColLabelsSet())
            {
                rowNumStart = rotatedBoundingBox.rowNumMin;
                colNumStart = rotatedBoundingBox.colNumMin;
                rowNumInc = rotatedBoundingBox.rowNumInc;
                colNumInc = rotatedBoundingBox.colNumInc;
            }
            else
            {
                rowNumStart = 0;
                colNumStart = 0;
                rowNumInc = 1;
                colNumInc = 1;
            }
            String filename = project.getRootDirString() + "grid-seismic.txt." + getName() + exportType;
            File file = new File(filename);
            if (file.exists())
            {
                boolean overWrite = StsYesNoDialog.questionValue(currentModel.win3d, "File " + filename + " already exists. Do you wish to overwrite it?");
                if (!overWrite) return;
                String binaryFilename = project.getBinaryFullDirString() +
                    "grid.bin." + getName();
                File binaryFile = new File(binaryFilename);
                if (binaryFile != null) binaryFile.delete();
            }

            progressBarDialog.setLabelText("Exporting surface to file: " + filename);
            //                progressBarDialog.setLabelText("Exporting surface to file: " + filename + " in seimic format (x, y, t, line xline)");

            printWriter = new PrintWriter(new FileWriter(filename, false));
            double[] origin = getAbsoluteXY(xMin, yMin);
            double rowStartX = origin[0];
            double rowStartY = origin[1];
            float rowNum = rowNumStart;
            for (int row = 0; row < nRows; row++)
            {
                double x = rowStartX;
                double y = rowStartY;
                float colNum = colNumStart;
                for (int col = 0; col < nCols; col++)
                {
                    if (pointsNull[row][col] == SURF_PNT)
                        printWriter.println(x + " " + y + " " + exportZ[row][col] + " " + rowNum + " " + colNum);
                    x += colXInc;
                    y += colYInc;
                    colNum += colNumInc;
                }
                rowStartX += rowXInc;
                rowStartY += rowYInc;
                rowNum += rowNumInc;
                progressBarDialog.setProgress(row);
            }
            /*
                   else
                   {
                       String filename = project.getRootDirString() + "grid.txt." + getName();
                       File file = new File(filename);
                       if (file.exists())
                       {
                           boolean overWrite = StsMessage.questionValue(currentModel.win3d, "File " + filename + " already exists. Do you wish to overwrite it?");
                           if (!overWrite)return;
                           String binaryFilename = project.getBinaryFullDirString() + "grid.bin." + name;
                           File binaryFile = new File(binaryFilename);
                           if (binaryFile != null)binaryFile.delete();
                       }

                       progressBarDialog.setLabelText("Exporting surface to file: " + filename);
       //                progressBarDialog.setLabelText("Exporting surface to file: " + filename + " in surface format (header plus z array)");
                       progressBarDialog.pack();
                       progressBarDialog.setVisible(true);

                       printWriter = new PrintWriter(new FileWriter(filename, false));

                       printWriter.println(xOrigin + " " + xInc + " " + yOrigin + " " + yInc + " " + nCols + " " + nRows + " ROW " + angle);
                       int nValues = 0;
                       for (int row = 0; row < nRows; row++)
                       {
                           for (int col = 0; col < nCols; col++)
                           {
                               if (pointsNull[row][col] != NOT_NULL)
                               {
                                   printWriter.print(nullZValue + " ");
                               }
                               else
                               {
                                   printWriter.print(exportZ[row][col] + " ");
                               }
                               if (++nValues == 5)
                               {
                                   printWriter.println();
                                   nValues = 0;
                               }
                           }
                           if (nValues > 0)
                           {
                               printWriter.println();
                           }
                           progressBarDialog.setProgress(row);
                       }
                       return;
                   }
            */
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurface.export() failed.", e, StsException.WARNING);
            return;
        }
        finally
        {
            if (printWriter != null)
            {
                printWriter.flush();
                printWriter.close();
            }
            // turn off the wait cursor
            //			if(cursor != null) cursor.restoreCursor();
            progressBarDialog.setProgress(maxProgress);
            progressBarDialog.dispose();
        }
    }

    /** Starting with current boundingBox, adjust x, y, and z ranges so they are "nice".
     * @parameters grid a new grid used to adjust the bounding box
     * @see StsRotatedBoundingBox
     * @see addBoundingBox
     * @see adjustBoundingBoxGridLines
     *  */
    /*  moved to StsProject
    public boolean addToProject()
    {
     StsProject project = currentModel.getProject();
     boolean originChanged = project.checkSetOrigin(xOrigin, yOrigin);
     if(!project.checkSetAngle(angle)) return false;
     StsRotatedGridBoundingBox boundingBox = project.getRotatedBoundingBox();

           float[] xy = project.getRotatedRelativeXY(xOrigin, yOrigin);
           xMin = xy[0];
           yMin = xy[1];
           xMax = xMin + (nCols-1)*xInc;
           yMax = yMin + (nRows-1)*yInc;
           xOrigin = project.getXOrigin();
           yOrigin = project.getYOrigin();
           boundingBox.checkMakeCongruent(this);
     float zDif = zMax - zMin;
     boundingBox.adjustZRange(zMin - zDif/2, zMax + zDif/2);
          boundingBox.adjustRowColNumbering();
//        boundingBox.addVolumeBoundingBox(this);
     StsBoundingBox unrotatedBoundingBox = project.getUnrotatedBoundingBox();
     unrotatedBoundingBox.addRotatedBoundingBox(boundingBox);

           setRelativeRotationAngle();
           project.resetCropVolume();
           return true;
    }
    */
    /**
     * surface may be rotated relative to project, so get the relative rotation angle
     * and rotate the projection matrix and draw as unrotated grid
     */
    public void setRelativeRotationAngle()
    {
        float projectAngle = currentModel.getProject().getAngle();
        relativeRotationAngle = angle - projectAngle;
        if (Math.abs(relativeRotationAngle) < 1.0f)
        {
            relativeRotationAngle = 0.0f;
        }
    }

    public void itemStateChanged(ItemEvent e)
    {
        if (e.getItem() instanceof StsColorscale)
        {
            if (displayProperty && newSurfaceTexture != null) textureChanged = true;
            colorListChanged = true;
        }
    }

    private void updateColors()
    {

    }

    public StsObject[] getPropertyList()
    {
        return currentModel.getObjectList(StsSeismicVolume.class);
    }

    public void setPropertyList(StsObject[] objectList)
    {
    }
    /*
        private void addDepthChangeListener()
        {
            depthColorscale.addItemListener
                (
                    new ItemListener()
            {
                public void itemStateChanged(ItemEvent e)
                {
                    if (e.getItem() instanceof StsColorscale)
                    {
                        depthColorscaleChanged();
                    }
                }
            }
            );
            return;
        }

        protected void depthColorscaleChanged()
        {
            float[] editRange = depthColorscale.getEditRange();
            setTextureChanged();
        }
    */


    /*
        private void addGradChangeListener()
        {
            gradColorscale.addItemListener
                (
                    new ItemListener()
            {
                public void itemStateChanged(ItemEvent e)
                {
                    if (e.getItem() instanceof StsColorscale)
                    {
                        gradColorscaleChanged();
                    }
                }
            }
            );
            return;
        }

        protected void gradColorscaleChanged()
        {
            float[] editRange = gradColorscale.getEditRange();
            setTextureChanged();
        }
    */
    public StsColorscale getColorscale()
    {
        return newSurfaceTexture.getColorscale();
    }


    public void initializeColorscaleActionListener()
    {
    }

    public float[] getHistogram()
    {
        return newSurfaceTexture.getHistogram();
    }

    public void setColorscale(StsColorscale colorscale)
    {
        if (colorscale == null)
        {
            return;
        }
        setTextureChanged();
        colorListChanged();
        currentModel.win3dDisplayAll();
    }

    public void colorListChanged()
    {
        colorListChanged = true;
//        currentModel.viewObjectRepaint(this);
    }

    public void colorListChanged(StsSurfaceTexture changedTexture)
    {
        if (changedTexture == currentSurfaceTexture)
            colorListChanged();
    }

    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            textureChanged = true;
            colorListChanged = true;
            currentModel.viewObjectRepaint(this, this);
        }
    }

    public void interpolateNeighbors(int row, int col, float z)
    {
        if(row > 0)
        {
            if(pointsNull[row-1][col] != SURF_PNT) pointsZ[row-1][col] = z;
            if(col > 0 && pointsNull[row-1][col-1] != SURF_PNT) pointsZ[row-1][col-1] = z;
            if(col < nCols-1 && pointsNull[row-1][col+1] != SURF_PNT) pointsZ[row-1][col+1] = z;
        }
        if(row < nRows-1)
        {
            if(pointsNull[row+1][col] != SURF_PNT) pointsZ[row+1][col] = z;
            if(col > 0 && pointsNull[row+1][col-1] != SURF_PNT) pointsZ[row+1][col-1] = z;
            if(col < nCols-1 && pointsNull[row+1][col+1] != SURF_PNT) pointsZ[row+1][col+1] = z;
        }
        if(col > 0 && pointsNull[row][col-1] != SURF_PNT) pointsZ[row][col-1] = z;
        if(col < nCols-1 && pointsNull[row][col+1] != SURF_PNT) pointsZ[row][col+1] = z;
    }

    public void checkClearNeighbors(int row, int col)
    {
        if(row > 0)
        {
            if(pointsNull[row-1][col] != SURF_PNT) pointsZ[row-1][col] = nullValue;
            if(col > 0 && pointsNull[row-1][col-1] != SURF_PNT) pointsZ[row-1][col-1] = nullValue;
            if(col < nCols-1 && pointsNull[row-1][col+1] != SURF_PNT) pointsZ[row-1][col+1] = nullValue;
        }
        if(row < nRows-1)
        {
            if(pointsNull[row+1][col] != SURF_PNT) pointsZ[row+1][col] = nullValue;
            if(col > 0 && pointsNull[row+1][col-1] != SURF_PNT) pointsZ[row+1][col-1] = nullValue;
            if(col < nCols-1 && pointsNull[row+1][col+1] != SURF_PNT) pointsZ[row+1][col+1] = nullValue;
        }
        if(col > 0 && pointsNull[row][col-1] != SURF_PNT) pointsZ[row][col-1] = nullValue;
        if(col < nCols-1 && pointsNull[row][col+1] != SURF_PNT) pointsZ[row][col+1] = nullValue;
    }
    public Class getDisplayableClass() { return StsSurface.class; }

    class SurfaceDefaultTexture extends StsSurfaceTexture
    {
        SurfaceDefaultTexture(StsSurface surface)
        {
            super(surface);
        }

        public byte[] getTextureData()
        {
            byte[] textureData = new byte[nRows * nCols];
            if (pointsNull == null)
            {
                return textureData;
            }

            int n = 0;
            for (int row = 0; row < nRows; row++, n += nCols)
            {
                System.arraycopy(pointsNull[row], 0, textureData, n, nCols);
            }
            return textureData;
        }

        public boolean isDisplayable() { return true; }

        public String getName()
        {
            return surface.getName();
        }

        public String toString()
        {
            return getName();
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            if (colorDisplayListNum != 0)
            {
                gl.glDeleteLists(colorDisplayListNum, 1);
            }
            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsException.systemError(this, "getColorDisplayListNum", "Failed to allocate a display list");
                return 0;
            }

            //        boolean nullsFilled = ((StsSurfaceClass)getCreateStsClass()).getNullsFilled();
            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
            //        arrayRGBA = null;
        }
/*
        public void createColorTLUT(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
            arrayRGBA = null;
        }
*/
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }

        /**
         * null values (for now) will be drawn transparent,
         * so we only need to set the single non-null value
         */
        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            float[] rgbaColor = getStsColor().getRGBA();
            float[][] arrayRGBA = new float[4][maxNullIndex + 1];

            // make NOT_NULL index the color of the surface
            // the shader can't handle this sort of ramp, so fill...
            for (int i = 0; i < 4; i++)
            {
                arrayRGBA[i][0] = rgbaColor[i];
                arrayRGBA[i][1] = rgbaColor[i];
            }

            // make all other nulls the null fill color
            if (nullsFilled)
            {
                String nullColorName = ((StsSurfaceClass) getCreateStsClass()).getNullColorName();
                StsColor nullColor = ((StsSurfaceClass) getCreateStsClass()).getNullColor(stsColor);
                if (nullColor == null)
                {
                    return arrayRGBA;
                }

                float[] rgbaGray = nullColor.getRGBA();
                for (int i = 0; i < 4; i++)
                {
                    for (int n = 2; n <= maxNullIndex; n++)
                    {
                        arrayRGBA[i][n] = rgbaGray[i];
                    }
                }
            }

            return arrayRGBA;
        }

        public StsColorscale getColorscale()
        {
            return null;
        }

        public void initializeColorscaleActionListener()
        {
        }

        public float[] getHistogram()
        {
            return null;
        }
    }

    class SurfaceDepthTexture extends StsSurfaceTexture
    {
        StsColorscale depthColorscale = null;
        int colorDisplayListNum = 0;

        SurfaceDepthTexture(StsSurface surface)
        {
            super(surface);
            initializeColorscaleActionListener();
        }

        public String getName()
        {
            return surface.getName() + "-depth";
        }

        public String toString()
        {
            return getName();
        }

        public StsColorscale getColorscale()
        {
            if (depthColorscale == null)
            {
                createDepthColorscale();
            }
            return depthColorscale;
        }

        private void createDepthColorscale()
        {
            double[] scale = StsMath.niceScale((double) getZMin(), (double) getZMax(), 32, true);
            float zMin = (float) scale[0];
            float zMax = (float) scale[1];

            StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
            depthColorscale = new StsColorscale("Depth", spectrum, zMin, zMax);
            depthColorscale.addActionListener(surface);
        }

        public void initializeColorscaleActionListener()
        {
            StsColorscale colorscale = getColorscale();
            colorscale.addActionListener(surface);
        }

        public float[] getHistogram()
        {
            return StsToolkit.buildHistogram(getTextureData(), getZMin(), getZMax());
        }

        public boolean isDisplayable()
        {
            return surface.isPersistent();
        }

        public byte[] getTextureData()
        {
            float zMin = surface.zMin;
            float zMax = surface.zMax;

            float[][] depths = surface.getPointsZ();
            byte[][] surfaceNulls = surface.getPointsNull();

            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            byte[] textureData = new byte[nRows * nCols];
            int n = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (surfaceNulls[row][col] != StsSurface.SURF_PNT)
                    {
                        textureData[n++] = -1;
                    }
                    else
                    {
                        byte byteIndex = (byte) (255 * ((depths[row][col] - zMin) / (zMax - zMin)));
                        textureData[n++] = byteIndex;
                    }
                }
            }
            //        if(getColorscale() != null)
            //        	StsToolkit.buildHistogram(textureData, surface.getZMin(), surface.getZMax());
            return textureData;
        }

        public void deleteColorDisplayList(GL gl)
        {
            if (colorDisplayListNum == 0)
            {
                return;
            }
            gl.glDeleteLists(colorDisplayListNum, 1);
            colorDisplayListNum = 0;
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            deleteColorDisplayList(gl);

            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsMessageFiles.logMessage(
                    "System Error in StsSurface.SurfaceDepthTexture.getColorListNum(): Failed to allocate a display list");
                return 0;
            }

            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
            arrayRGBA = null;
        }
/*
        public void createColorTLUT(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
        }
*/
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }

        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            Color[] colors = getColorscale().getNewColorsInclTransparency();

            if (nullsFilled)
            {
                String nullColorName = ((StsSurfaceClass) surface.
                    getCreateStsClass()).getNullColorName();
                StsColor nullColor = ((StsSurfaceClass) surface.getCreateStsClass()).
                    getNullColor(surface.getStsColor());
                colors[colors.length - 1] = nullColor.getColor();
            }

            int nColors = colors.length;
            float[][] arrayRGBA = new float[4][nColors];
            float[] rgba = new float[4];
            for (int n = 0; n < nColors; n++)
            {
                colors[n].getComponents(rgba);
                for (int i = 0; i < 4; i++)
                {
                    arrayRGBA[i][n] = rgba[i];
                }
            }
            return arrayRGBA;
        }

        public float getDataMin() { return surface.getZMin(); }
        public float getDataMax() { return surface.getZMax(); }
    }

    class SurfaceGradTexture extends StsSurfaceTexture
    {
        StsColorscale gradColorscale = null;
        int colorDisplayListNum = 0;
        float gradRangeMax;
        byte[] textureData = null;

        SurfaceGradTexture(StsSurface surface)
        {
            super(surface);
            initializeColorscaleActionListener();
        }

        public String getName()
        {
            return surface.getName() + "-grad";
        }

        public String toString()
        {
            return getName();
        }

        public byte[] computeTextureData()
        {
            float[][] grad = new float[nRows][nCols];
            byte[] textureData = new byte[nRows * nCols];

            // put texture nulls in edge rows and cols
            gradRangeMax = -largeFloat;
            for (int col = 0; col < nCols; col++)
            {
                textureData[col] = -1;
            }
            int n = (nRows - 1) * nCols;
            textureData[n++] = -1;

            int n1 = 0;
            int n2 = nCols - 1;
            for (int row = 0; row < nRows; row++, n1 += nCols, n2 += nCols)
            {
                textureData[n1] = -1;
                textureData[n2] = -1;
            }

            for (int row = 1; row < nRows - 1; row++)
            {
                for (int col = 1; col < nCols - 1; col++)
                {
                    if (pointsNull[row][col] == StsSurface.SURF_PNT)
                    {
                        float zc = pointsZ[row][col];
                        float maxGrad = 0.0f;
                        if (pointsNull[row - 1][col - 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(halfSqrt2f *
                                    (pointsZ[row - 1][col - 1] - zc)));
                        }
                        if (pointsNull[row - 1][col + 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(halfSqrt2f *
                                    (pointsZ[row - 1][col + 1] - zc)));
                        }
                        if (pointsNull[row + 1][col - 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(halfSqrt2f *
                                    (pointsZ[row + 1][col - 1] - zc)));
                        }
                        if (pointsNull[row + 1][col + 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(halfSqrt2f *
                                    (pointsZ[row + 1][col + 1] - zc)));
                        }
                        if (pointsNull[row - 1][col] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(pointsZ[row - 1][col] - zc));
                        }
                        if (pointsNull[row + 1][col] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(pointsZ[row + 1][col] - zc));
                        }
                        if (pointsNull[row][col - 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(pointsZ[row][col - 1] - zc));
                        }
                        if (pointsNull[row][col + 1] == StsSurface.SURF_PNT)
                        {
                            maxGrad = Math.max(maxGrad,
                                Math.abs(pointsZ[row][col + 1] - zc));
                        }
                        if (maxGrad != 0.0f)
                        {
                            grad[row][col] = maxGrad;
                            gradRangeMax = Math.max(gradRangeMax, maxGrad);
                        }
                        else
                        {
                            grad[row][col] = -1.0f;
                        }
                    }
                }
            }
            for (int row = 1; row < nRows - 1; row++)
            {
                n = row * nCols + 1;
                for (int col = 1; col < nCols - 1; col++)
                {
                    if (grad[row][col] != 0.0f)
                    {
                        textureData[n++] = (byte) (255 * grad[row][col] /
                            gradRangeMax);
                    }
                    else
                    {
                        textureData[n++] = -1;
                    }
                }
            }
            //        if(getColorscale() != null)
            //        	StsToolkit.buildHistogram(textureData, 0.0f, gradRangeMax);
            return textureData;
        }

        public void initializeColorscaleActionListener()
        {
            StsColorscale colorscale = getColorscale();
            colorscale.addActionListener(surface);
        }

        public StsColorscale getColorscale()
        {
            if (gradColorscale == null)
            {
                createGradColorscale();
            }
            return gradColorscale;
        }

        private void createGradColorscale()
        {
            double[] scale = StsMath.niceScale(0.0, (double) gradRangeMax, 32, true);
            StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
            gradColorscale = new StsColorscale("Gradient", spectrum, 0.0f, (float) scale[1]);
            gradColorscale.addActionListener(surface);
        }

        public float[] getHistogram()
        {
            byte[] textureData = getTextureData();
            return StsToolkit.buildHistogram(textureData, getZMin(), getZMax());
        }

        public boolean isDisplayable()
        {
            return surface.isPersistent();
        }

        public byte[] getTextureData()
        {
            if(textureData == null)
                textureData = computeTextureData();
            return textureData;
        }

        public void deleteColorDisplayList(GL gl)
        {
            if (colorDisplayListNum == 0)
            {
                return;
            }
            gl.glDeleteLists(colorDisplayListNum, 1);
            colorDisplayListNum = 0;
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            deleteColorDisplayList(gl);

            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsMessageFiles.logMessage(
                    "System Error in StsSurface.SurfaceDepthTexture.getColorListNum(): Failed to allocate a display list");
                return 0;
            }

            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
            arrayRGBA = null;
        }
/*
        public void createColorTLUT(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
        }
*/
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }

        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            Color[] colors = getColorscale().getNewColorsInclTransparency();
            if (nullsFilled)
            {
                String nullColorName = ((StsSurfaceClass) surface.
                    getCreateStsClass()).getNullColorName();
                StsColor nullColor = ((StsSurfaceClass) surface.getCreateStsClass()).
                    getNullColor(surface.getStsColor());
                colors[colors.length - 1] = nullColor.getColor();
            }
            int nColors = colors.length;
            float[][] arrayRGBA = new float[4][nColors];
            float[] rgba = new float[4];
            for (int n = 0; n < nColors; n++)
            {
                colors[n].getComponents(rgba);
                for (int i = 0; i < 4; i++)
                {
                    arrayRGBA[i][n] = rgba[i];
                }
            }
            return arrayRGBA;
        }

        public float getDataMin() { return 0.0f; }
        public float getDataMax() { return this.gradRangeMax; }
    }

    public boolean textureChanged()
    {
        textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void geometryChanged()
    {
    }

    public String getZDomainString() { return StsParameters.TD_STRINGS[zDomainOriginal]; }
    public void setZDomainString(String string) { zDomainOriginal = StsParameters.getStringMatchByteIndex(StsParameters.TD_STRINGS, string); }
}
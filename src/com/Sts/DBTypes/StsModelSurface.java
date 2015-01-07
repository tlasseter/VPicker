//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsModelSurface extends StsSurface implements StsSelectable, StsXYSurfaceGridable, StsTreeObjectI
{
    /**
     * surface from which this modelSurface was derived
     */
    protected StsSurface originalSurface = null;
    protected StsZone zoneAbove, zoneBelow;
    protected int subZoneNumber;
    /**
     * Display points which were null but if been filled by extrapolation/interpolation up to faults.
     */
    transient protected boolean displayGaps = true;

    transient protected StsList blockGrids;
    /**
     * horizon patches in each block
     */
    //TODO remove refactor propertyVolume
    transient protected StsPropertyVolumeOld propertyVolume = null;

    transient public byte[][] originalPointsNull;
    transient public float[][] originalPointsZ;
    transient public float[][] originalAdjPointsZ;
    transient public StsBlock[][] blockArray = null; // blockGrid at this ij

    transient public float[][] rowLinkDistances; // scratch array arc length in going from i,j to i,j+1
    transient public float[][] colLinkDistances; // scratch array arc length in going from i+1,j to i,j
    transient public float[][] faultDistances; // current distance at this i,j from point to interpolate (allocated in StsGrid)

    transient protected String propertyName = null;
    transient protected int propLayerNumber;

    /**
     * bidirectional double-link list of StsGridLinks
     */
    transient protected StsEdgeLoopRadialLinkGrid linkedGrid;

    // Used to gap voxelGrids back from faults
    transient StsFaultSurface faultSurface = null;

    static StsObjectPanel modelObjectPanel = null;

    // Class version ID: rerun serialver when the class is changed
    static final long serialVersionUID = -5294107068658877877L;

    static public final byte LINK_NULL = StsParameters.LINK_NULL;
    static public final byte LINK = StsParameters.LINK;
    static public final byte LINK_SPLIT = StsParameters.LINK_SPLIT;

    static public final int LOOP = StsParameters.LOOP;
    static public final int OPEN_LOOP = StsParameters.OPEN_LOOP;

    static StsComboBoxFieldBean modelSurfaceTextureListBean;
    static public StsFieldBean[] modelDisplayFields = null;
    static StsEditableColorscaleFieldBean modelColorscaleBean;
    static public StsFieldBean[] modelPropertyFields = null;

    /**
     * constructor for DB
     */
    public StsModelSurface()
    {
    }

    public StsModelSurface(boolean persistent)
    {
        super(persistent);
    }

    static public StsModelSurface constructModelSurface(String name, StsColor stsColor, byte type)
    {
        try
        {
            return new StsModelSurface(name, stsColor, type);
        }
        catch (Exception e)
        {
            StsException.outputException("StsModelSurface.construct() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /**
     * constructor for surface that may have nulls
     */
    private StsModelSurface(String name, StsColor stsColor, byte type)
    {
        setName(name);
        this.stsColor = stsColor;
        this.type = type;
    }

    static public StsModelSurface constructModelSurfaceFromSurface(String name, StsSurface surface, byte type)
    {
        try
        {
            return new StsModelSurface(name, surface, type);
        }
        catch (Exception e)
        {
            StsException.outputException("StsModelSurface.construct() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public StsModelSurface constructModelSurfaceFromSurface(StsSurface surface, StsModelSurfaceClass modelSurfaceClass, StsGridDefinition gridDef)
    {
        try
        {
            StsModelSurface modelSurface = modelSurfaceClass.getModelSurfaceForOriginal(surface);
            // new surface
            if (modelSurface == null)
            {
                String name = surface.getName() + "-Horizon";
                modelSurface = StsModelSurface.constructModelSurfaceFromSurface(name, surface, StsModelSurface.MODEL);
            }
            boolean hasNulls = surface.getHasNulls();
            byte zDomain = surface.getZDomainOriginal();
            byte zDomainSupported = zDomain;
            modelSurface.initialize(gridDef.getNCols(), gridDef.getNRows(), gridDef.getXOrigin(), gridDef.getYOrigin(),
                    gridDef.getXInc(), gridDef.getYInc(), gridDef.getXMin(), gridDef.getYMin(),
                    gridDef.getAngle(), null, true, StsParameters.nullValue, zDomain, zDomainSupported, null);
            boolean congruent = modelSurface.copyPoints(surface);
            modelSurface.initSurfaceTextureList();
            if (!congruent)
            {
                interpolate(modelSurface);
                modelSurface.constructGrid();
            }
            surface.unload(); // remove original surface from memory and make invisible
            return modelSurface;
        }
        catch (Exception e)
        {
            StsException.outputException("buildHorizonSurfaces failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    static public StsModelSurface constructModelSurfaceFromSurfaceAndOffset(StsModelSurface surface, float offset, byte offsetDomain, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            boolean hasNulls = surface.getHasNulls();
            byte zDomainOriginal = surface.getZDomainOriginal();
            byte zDomainSupported = surface.zDomainSupported;
            if (offset == 0.0f) return surface;
            String newSurfaceName;
            if (offset < 0.0f)
                newSurfaceName = surface.getName() + "-" + offset;
            else
                newSurfaceName = surface.getName() + "+" + offset;
            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            StsModelSurface offsetSurface = new StsModelSurface(newSurfaceName, surface.getStsColor(), INTERPOLATED, nCols, nRows,
                    surface.getXOrigin(), surface.getYOrigin(), surface.getXInc(), surface.getYInc(), surface.getXMin(), surface.getYMin(),
                    surface.getAngle(), null, hasNulls, StsParameters.nullValue, zDomainOriginal, zDomainSupported);
            if (zDomainOriginal == offsetDomain)
                offsetSurface.copyPointsZWithOffset(surface.pointsZ, offset);
            else
                offsetSurface.copyPointsZWithOffset(surface.pointsZ, surface.adjPointsZ, velocityModel, offset);
            offsetSurface.copyPointsNull(surface.pointsNull);
            return offsetSurface;
        }
        catch (Exception e)
        {
            StsException.outputException("buildHorizonSurfaces failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    static private void interpolate(StsModelSurface surface)
    {
        surface.interpolateDistanceTransform();
        //       StsInterpolationRadiusWeighted.constructor(surface).interpolate();
    }


    /**
     * constructor for a model surface built from an existing surface
     */
    private StsModelSurface(String name, StsSurface surface, byte type)
    {
        setName(name);
        this.stsColor = surface.getStsColor();
        this.marker = surface.getMarker();
        this.type = type;
        originalSurface = surface;
        zDomainSupported = surface.zDomainOriginal;
        this.surfaceTextureName = name;
    }

    /*
        static public StsModelSurface constructModelSurface(String name, StsColor stsColor, byte type,
            int nCols, int nRows, double xOrigin, double yOrigin, float xInc, float yInc,
            float xMin, float yMin, float angle, float[][] pointsZ, boolean hasNulls, float nullZValue, byte zDomain)
        {
            try
            {
                return new StsModelSurface(name, stsColor, type, nCols, nRows, xOrigin, yOrigin, xInc, yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue, zDomain);
            }
            catch (Exception e)
            {
                StsException.outputException("StsGrid.construc() failed.", e, StsException.WARNING);
                return null;
            }
        }
    */
    /**
     * constructor for surface that may have nulls
     */

    private StsModelSurface(String name, StsColor stsColor, byte type, int nCols, int nRows,
                            double xOrigin, double yOrigin, float xInc, float yInc, float xMin, float yMin, float angle,
                            float[][] pointsZ, boolean hasNulls, float nullZValue, byte zDomain, byte zDomainSupported)
            throws StsException
    {
        initialize(name, stsColor, type, nCols, nRows, xOrigin, yOrigin, xInc, yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue, zDomain, zDomainSupported);
        initSurfaceTextureList();
    }

    private StsModelSurface(boolean persistent, String name)
    {
        super(persistent, name);
    }

    static public StsModelSurface constructTempModelSurface(String name)
    {
        return new StsModelSurface(false, name);
    }

    public boolean initialize(String name, StsColor stsColor, byte type, int nCols, int nRows,
                              double xOrigin, double yOrigin, float xInc, float yInc, float xMin, float yMin, float angle,
                              float[][] pointsZ, boolean hasNulls, float nullZValue, byte zDomain, byte zDomainSupported)
    {
        // save inputs
        setName(name);
        this.stsColor = stsColor;
        this.type = type;
        return initialize(nCols, nRows, xOrigin, yOrigin, xInc, yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue, zDomain, zDomainSupported, null);
    }
    /*
        public boolean classInitialize(int nCols, int nRows, double xOrigin, double yOrigin,
            float xInc, float yInc, float angle, float[][] pointsZ, boolean hasNulls, float nullZValue)
        {
            this.nCols = nCols;
            this.nRows = nRows;
            this.nPnts = nRows*nCols;
            this.xOrigin = xOrigin;
            this.yOrigin = yOrigin;
            this.xInc = xInc;
            this.yInc = yInc;
            this.angle = angle;
            this.pointsZ = pointsZ;
            this.hasNulls = hasNulls;
            this.nullZValue = nullZValue;

            if (!checkGrid(currentModel.win3d))
                return false;

            checkResetOrigin();
            setAngle();
            initializePointsZ();

            StsProject project = currentModel.getProject();

            if (!project.isOriginSet())
            {
                project.setOriginAndAngle(this.xOrigin, this.yOrigin, this.angle);
                relativeXOrigin = 0.0f;
                relativeYOrigin = 0.0f;
            }
            else
            {
                float[] relativeOrigin = project.computeRelativeOrigin(this.xOrigin, this.yOrigin);
                if (relativeOrigin != null)
                {
                    relativeXOrigin = relativeOrigin[0];
                    relativeYOrigin = relativeOrigin[1];
                }
            }

            if (pointsZ != null)
            {
                StsStatusArea status = StsStatusArea.getStatusArea();
                StsMessageFiles.logMessage("Processing grid...");

                interpolateNullZPoints();
                constructGrid();
            }

    //		project.adjustRange(getRelativeRange());

    //        System.out.println("StsGrid: nulls ? " + hasNulls + " nx " + nX + " ny " + nY + "\n" +
    //        					"Xmin " + xMin + " Xmax " + xMax + " Xinc " + xInc + "\n" +
    //        					"Ymin " + yMin + " Ymax " + yMax + " Yinc " + yInc);
            return true;
       }

       private void setAngle()
       {
           if(angle == 0.0)
           {
               colXInc = xInc;
               colYInc = 0.0;
               rowXInc = 0.0;
               rowYInc = yInc;
           }
           else
           {
               cosAngle = StsMath.cosd(angle);
               sinAngle = StsMath.sind(angle);
               colXInc = cosAngle*xInc;
               colYInc = sinAngle*xInc;
               rowXInc = -sinAngle*yInc;
               rowYInc = cosAngle*yInc;
           }
       }
    */

    public boolean isDisplayable()
    {
        return isPersistent();
    }

    public void initModelSurfaceTextureList()
    {
        getSurfaceTextureList();
        if (modelSurfaceTextureListBean != null)
            modelSurfaceTextureListBean.setListItems(surfaceTextureList);
    }

    protected StsSurfaceTexture getCorrelCoefTexture()
    {
        StsHorpickClass horpickClass = (StsHorpickClass) currentModel.getStsClass(StsHorpick.class);
        if (horpickClass == null) return null;
        StsHorpick horpick = horpickClass.getHorpickWithSurface(originalSurface);
        if (horpick == null) return null;
        return horpick.getCorrelCoefsSurfaceTexture(this);
    }

    public StsList getBlockGrids()
    {
        return blockGrids;
    }

    //    public void setMarkerVisible() { ; }
    //    public boolean getMarkerVisible() { return false; }

    public StsMarker getMarker()
    {
        if (marker != null) return marker;
        if (originalSurface != null)
        {
            marker = originalSurface.getMarker();
        }
        return marker;
    }

    public void setStsColor(StsColor stsColor)
    {
        if (this.stsColor != null && this.stsColor.equals(stsColor)) return;
        if (marker != null) marker.setStsColor(stsColor);
        if (originalSurface != null) originalSurface.setStsColor(stsColor);
        super.setStsColor(stsColor);
        // currentModel.viewObjectRepaint(this, this);
    }

    public boolean getDisplayGaps()
    {
        return displayGaps;
    }

    public void setDisplayGaps(boolean b)
    {
        if (displayGaps == b) return;
        displayGaps = b;
        textureChanged();
        currentModel.win3dDisplayAll();
    }

    public StsSurface getOriginalSurface()
    {
        return originalSurface;
    }

    public void setOriginalSurface(StsSurface surface)
    {
        originalSurface = surface;
    }

    public void setZoneAbove(StsZone zoneAbove)
    {
        this.zoneAbove = zoneAbove;
        dbFieldChanged("zoneAbove", zoneAbove);
    }

    public StsZone getZoneAbove()
    {
        return zoneAbove;
    }

    public void setZoneBelow(StsZone zoneBelow)
    {
        this.zoneBelow = zoneBelow;
        dbFieldChanged("zoneBelow", zoneBelow);
    }

    public StsZone getZoneBelow()
    {
        return zoneBelow;
    }

    public void setZone(StsZone zone)
    {
        this.zoneAbove = zone;
        this.zoneBelow = zone;
    }

    public StsZone getZone()
    {
        return zoneAbove;
    }

    public void setSubZoneNumber(int number)
    {
        subZoneNumber = number;
        dbFieldChanged("subZoneNumber", subZoneNumber);
    }

    public int getSubZoneNumber()
    {
        return subZoneNumber;
    }

    public StsFaultSurface getFaultSurface()
    {
        return faultSurface;
    }

    public StsEdgeLoopRadialLinkGrid getLinkedGrid()
    {
        return linkedGrid;
    }

    public StsBlockGrid getBlockGrid(int row, int col)
    {
        StsBlock block = getBlock(row, col);
        if (block == null) return null;
        return block.getBlockGrid(this);
    }

    public StsBlock getBlock(int row, int col)
    {
        if (blockArray == null) return null;
        if (!insideGrid(row, col)) return null;
        if (blockArray[row] == null) return null;
        return blockArray[row][col];
    }

    public StsModelSurface cloneModelSurface()
    {
        StsModelSurface newSurface = new StsModelSurface();
        newSurface.copyModelSurface(this);
        return newSurface;
    }

    public StsModelSurface copyModelSurface(StsModelSurface surface)
    {
        super.copySurface(surface);
        return this;
    }

    /**
     * set grid points using inverse distance interpolation from well markers
     */
    public boolean setGridPoints(StsMarker marker)
    {
        return setGridPoints(marker, null, StsParameters.largeFloat, 1.0f);
    }

    public boolean setGridPoints(StsMarker marker, StsModelSurface surface, float maxDist, float power)
    {
        if (marker == null) return false;
        StsObjectRefList wellMarkers = marker.getWellMarkers();
        int nWellMarkers = (wellMarkers == null) ? 0 : wellMarkers.getSize();
        if (nWellMarkers == 0) return false;

        if (pointsZ == null) pointsZ = new float[nRows][nCols];
        if (pointsNull == null) pointsNull = new byte[nRows][nCols];
        if (nWellMarkers == 1)    // trivial case
        {
            StsWellMarker wellMarker = (StsWellMarker) wellMarkers.getElement(0);
            float z = wellMarker.getZ();
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    pointsZ[row][col] = z;
                    setPointNull(row, col, SURF_PNT);
                }
            }
            return true;
        }

        // extract marker x-y-z values
        float[] markerX = new float[nWellMarkers];
        float[] markerY = new float[nWellMarkers];
        float[] markerZ = new float[nWellMarkers];
        for (int i = 0; i < nWellMarkers; i++)
        {
            StsWellMarker wellMarker = (StsWellMarker) wellMarkers.getElement(i);
            StsPoint xyz = wellMarker.getLocation();
            markerX[i] = xyz.getX();
            markerY[i] = xyz.getY();
            markerZ[i] = xyz.getZ();
        }

        // go thru all grid points
        setGridPoints(markerX, markerY, markerZ, maxDist, power);
        for (int row = 0; row < nRows; row++)
            for (int col = 0; col < nCols; col++)
                setPointNull(row, col, SURF_PNT);

        // set null point values
        byte[][] gridPointsNull = getPointsNull();
        for (int row = 0; row < nRows; row++)
            for (int col = 0; col < nCols; col++)
                setPointNull(row, col, gridPointsNull[row][col]);

        return true;
    }

    /**
     * create a grid from a well marker isopach above/below a surface
     */
    public boolean setGridPoints(StsModelSurface surface, StsMarker marker,
                                 float maxDist, float power)
    {
        if (surface == null) return false;
        if (marker == null) return false;
        StsObjectRefList wellMarkers = marker.getWellMarkers();
        int nWellMarkers = (wellMarkers == null) ? 0 : wellMarkers.getSize();
        if (nWellMarkers == 0) return false;

        StsMarker surfaceMarker = getMarker();
        if (surfaceMarker == null) return false;
        StsObjectRefList surfaceWellMarkers = surfaceMarker.getWellMarkers();
        int nSurfaceWellMarkers = (surfaceWellMarkers == null) ? 0
                : surfaceWellMarkers.getSize();
        if (nSurfaceWellMarkers == 0) return false;

        if (pointsZ == null) pointsZ = new float[nRows][nCols];
        if (pointsNull == null) pointsNull = new byte[nRows][nCols];

        // extract x-y-dz values
        StsList markerValues = new StsList(Math.min(nWellMarkers, nSurfaceWellMarkers));
        for (int i = 0; i < nWellMarkers; i++)
        {
            StsWellMarker wellMarker = (StsWellMarker) wellMarkers.getElement(i);
            StsWell well = wellMarker.getWell();
            for (int j = 0; j < nSurfaceWellMarkers; j++)
            {
                StsWellMarker surfaceWellMarker = (StsWellMarker) surfaceWellMarkers.getElement(j);
                StsPoint p = null;
                if (surfaceWellMarker.getWell() == well)
                {
                    p = new StsPoint(wellMarker.getLocation());
                    p.setZ(p.getZ() - surfaceWellMarker.getZ());
                    markerValues.add(p);
                }
            }
        }
        int nMarkerValues = markerValues.getSize();
        if (nMarkerValues == 0) return false;


        if (nMarkerValues == 1)    // trivial case
        {
            StsPoint p = (StsPoint) markerValues.getElement(0);
            float z = p.getZ();
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    pointsZ[row][col] = z;
                    setPointNull(row, col, SURF_PNT);
                }
            }
            return true;
        }

        // put points in arrays
        float[] x = new float[nMarkerValues];
        float[] y = new float[nMarkerValues];
        float[] dz = new float[nMarkerValues];
        for (int i = 0; i < nMarkerValues; i++)
        {
            StsPoint p = (StsPoint) markerValues.getElement(i);
            x[i] = p.getX();
            y[i] = p.getY();
            dz[i] = p.getZ();
        }

        // go thru all grid points
        setGridPoints(x, y, dz, maxDist, power);

        // add surface point values
        float[][] surfacePoints = getPointsZ();
        byte[][] surfacePointsNull = getPointsNull();

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                pointsZ[row][col] += surfacePoints[row][col];
                setPointNull(row, col, surfacePointsNull[row][col]);
            }
        }

        return true;
    }


    /**
     * we can unload imported surfaces only: too dangerous for model surfaces
     */
    public void unload()
    {
        setIsVisible(false);
        if (getType() != StsModelSurface.IMPORTED) return;

        // these are the only arrays currently used by imported surfaces
        pointsZ = null;
        pointsNull = null;
        normals = null;
    }

    /*
        public void setIsVisible(boolean b)
        {
            if(b == isVisible) return;
            super.setIsVisible(b);

            if(b) isVisible = pointsLoaded(type);
            else isVisible = b;

            if(b)
            {
                StsModelSurfaceClass surfaceClass = getModelSurfaceClass();
                surfaceClass.setDisplaySurfaces(true);
            }
    //		currentModel.win3dDisplayAll();
        }
    */
    /*
        public void copyPoints(StsModelSurface surface)
        {
            super.copyPoints(surface);
        }
    */
    public void displaySurfacePosition(StsGridPoint gridPoint)
    {
        float[] xyz = gridPoint.getXYZorT();

        double[] absXY = currentModel.getProject().getAbsoluteXYCoordinates(xyz);

        // Get ij location info
        //        StsGridPoint gridPoint = new StsGridPoint(p, this);
        //        getUnboundedGridFloatCoordinates(gridPoint);

        // display x-y-z location

        int i = gridPoint.row;
        int j = gridPoint.col;
        int n = getIndex(i, j);

        String message;
        message = getName() + "  - True X: " + (int) (absXY[0]) +
                "  True Y: " + (int) (absXY[1]) +
                "  X: " + (int) xyz[0] +
                "  Y: " + (int) xyz[1] + "  Z: " + (int) xyz[2] + " row: " + i + " col: " + j;
        StsMessageFiles.infoMessage(message);

       if (getPropertyVolume() != null)
       {
            String propertyName = getPropertyName();
            String name = propertyName == null ? "  Property: " : "  " + propertyName + ": ";
            float value = getPropValue(xyz);
            message = "    " + name + value + "\n";
            StsMessageFiles.infoMessage(message);
        }
        if(!hasBlockGrids())
        {
            message = "    " + getSurfacePointString(i, j) + "\n";
            StsMessageFiles.infoMessage(message);
        }
        else
        {
            int nBlockGrids = blockGrids.getSize();
            for(int b = 0; b < nBlockGrids; b++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid)blockGrids.getElement(b);
                message = "    " + blockGrid.getPointString(i, j);
                StsMessageFiles.infoMessage(message);
            }
        }
    }

    public boolean hasBlockGrids() { return blockGrids != null && blockGrids.getSize() > 0; }

    public String getBlockGridPointsString(int row, int col)
    {
        String string =  "   row: " + row + " col: " + col;
        int nBlockGrids = blockGrids.getSize();
        for(int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid)blockGrids.getElement(n);
            string = string.concat(blockGrid.getPointString(row, col));
        }
        return string;
    }
    public void initializeBlockGridArray()
    {
        try
        {
            if (blockArray != null) return;
            blockArray = new StsBlock[nRows][];

            int rowPlusLinkIndex = StsEdgeLoopRadialGridLink.ROW_PLUS;
            StsList blockGrids = getBlockGrids();
            if (blockGrids == null) return;
            int nBlockGrids = blockGrids.getSize();
            for (int g = 0; g < nBlockGrids; g++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(g);
                StsList links = blockGrid.getEdgeLoop().getLoopLinks();
                if (links == null) continue;
                int nLinks = links.getSize();
                for (int n = 0; n < nLinks; n++)
                {
                    StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) links.getElement(n);
                    StsEdgeLoopRadialGridLink nextLink = link.getConnectedLink(rowPlusLinkIndex);
                    if (nextLink != null) initializeBlockArray(link, nextLink, blockGrid);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsFaultSurface.initializeBlockGridArray() failed.",
                    e, StsException.WARNING);
        }
    }

    public void initializePointArrays()
    {
        if(originalPointsNull == null)
        {
            originalPointsNull = StsMath.copyByteArray(pointsNull);
            originalPointsZ = StsMath.copyFloatArray(pointsZ);
            if(adjPointsZ != null)
                originalAdjPointsZ = StsMath.copyFloatArray(adjPointsZ);
        }
        else
        {
            pointsNull = StsMath.copyByteArray(originalPointsNull);
            pointsZ = StsMath.copyFloatArray(originalPointsZ);
            if(originalAdjPointsZ != null)
                adjPointsZ = StsMath.copyFloatArray(originalAdjPointsZ);
        }
    }

    private void initializeBlockArray(StsEdgeLoopRadialGridLink link, StsEdgeLoopRadialGridLink nextLink, StsBlockGrid blockGrid)
    {
        if (link == null || nextLink == null) return;
        if (blockGrid == null) return;

        int row = (int) link.getRowF();
        int colStart = StsMath.ceiling(link.getColF());
        int colEnd = StsMath.floor(nextLink.getColF());

        for (int col = colStart; col <= colEnd; col++)
            setBlock(row, col, blockGrid.block);

        blockGrid.adjustRowColRange(row, colStart);
        blockGrid.adjustRowColRange(row, colEnd);
    }

    public void setBlock(int row, int col, StsBlock block)
    {
        if (blockArray == null) return;
        StsBlock[] blockRow = blockArray[row];
        if (blockRow == null)
        {
            blockRow = new StsBlock[nCols];
            blockArray[row] = blockRow;
        }
        blockRow[col] = block;
    }

    public void computeFaultArcLengthWeights()
    {
        StsList blockGrids = getBlockGrids();

        rowLinkDistances = new float[nRows][];
        colLinkDistances = new float[nCols][];
        faultDistances = new float[nRows][nCols];

        int nBlockGrids = blockGrids.getSize();
        for (int g = 0; g < nBlockGrids; g++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(g);
            StsList edges = blockGrid.getEdgeLoop().getEdges();

            if (!blockGrid.hasAuxiliaryEdges())
                setClosedLoopArcLengthWeights(edges);
            else
            {
                int nEdges = edges.getSize();

                // now for each auxiliaryEdge, find the next one and build point lengths
                // for edges in between

                StsSurfaceEdge edge, nextEdge;

                for (int e = 0; e < nEdges; e++)
                {
                    edge = (StsSurfaceEdge) edges.getElement(e);
                    nextEdge = (StsSurfaceEdge) (edge.getNextConnectedSurfaceEdge());

                    // find break between a non-fault edge and the next fault edge
                    if (!edge.isFaulted() && nextEdge != null && nextEdge.isFaulted())
                    {
                        StsList faultEdges = new StsList(4, 4);
                        while (nextEdge != null && nextEdge.isFaulted())
                        {
                            faultEdges.add(nextEdge);
                            nextEdge = (StsSurfaceEdge) (nextEdge.getNextConnectedSurfaceEdge());
                        }
                        computeFaultArcLengthWeights(faultEdges);
                    }
                }
            }
        }
    }

    private void computeFaultArcLengthWeights(StsList edges)
    {
        int e;

        int nEdges = edges.getSize();
        if (nEdges == 0) return;

        StsSurfaceEdge firstEdge = (StsSurfaceEdge) edges.getFirst();
        StsSurfaceEdge prevEdge = (StsSurfaceEdge) firstEdge.getPrevConnectedSurfaceEdge();
        boolean firstIsDyingFault = prevEdge == null || prevEdge.isAuxiliary();

        StsSurfaceEdge lastEdge = (StsSurfaceEdge) edges.getLast();
        StsSurfaceEdge nextEdge = (StsSurfaceEdge) lastEdge.getNextConnectedSurfaceEdge();
        boolean lastIsDyingFault = nextEdge == null || nextEdge.isAuxiliary();

        // first compute arc lengths for each edge individually and get sum of all edge lengths
        float sumLengths = 0.0f;
        for (e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(e);
            edge.computeArcLengths();
            float edgeLength = edge.getTotalArcLength();
            sumLengths = sumLengths + edgeLength;
        }

        int zeroEnd;
        if (firstIsDyingFault && !lastIsDyingFault)
            zeroEnd = MINUS;
        else if (!firstIsDyingFault && lastIsDyingFault)
            zeroEnd = PLUS;
        else
            zeroEnd = StsParameters.PLUS_AND_MINUS;

        float runningSumLengths = 0.0f;
        for (e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(e);
            runningSumLengths = edge.adjustArcLengths(zeroEnd, runningSumLengths, sumLengths);
        }

        // now for each edge, put arc lengths in row and col faultArcLengths arrays depending on
        // position of each edge point

        for (e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(e);
            StsList edgePoints = edge.getGridEdgePointsList();
            float[] arcLengths = edge.getArcLengths();
            if (arcLengths == null)
            {
                StsException.outputException(new StsException(StsException.WARNING,
                        "StsGrid.computeFaultArcLengthWeights() failed.",
                        "arcLengths have not been computed"));
                continue;
            }

            int nEdgePoints = edgePoints.getSize();

            if (arcLengths.length != nEdgePoints)
            {
                StsException.outputException(new StsException(StsException.WARNING,
                        "StsGrid.computeFaultArcLengthWeights() failed.",
                        "weights length: " + arcLengths.length +
                                " not equal to number of edgePoints: " + nEdgePoints));
                continue;
            }

            for (int p = 0; p < nEdgePoints; p++)
            {
                StsGridSectionPoint point = (StsGridSectionPoint) edgePoints.getElement(p);
                int row = (int) point.getRowF(null);
                int col = (int) point.getColF(null);
                int rowOrCol = point.getRowOrCol(null);
                if (rowOrCol == ROW)
                    setPointLinkDistance(row, col, rowLinkDistances, arcLengths[p]);
                else
                    setPointLinkDistance( col, row, colLinkDistances, arcLengths[p]);
            }
        }
    }

    private void setPointLinkDistance(int row, int col, float[][] linkDistances, float arcLength)
    {
        if (linkDistances[row] == null) linkDistances[row] = new float[nCols];
        float currentDistance = linkDistances[row][col];

        if (currentDistance == 0.0f)
            linkDistances[row][col] = 2 * arcLength;
        else
            linkDistances[row][col] = Math.max(currentDistance, 2 * arcLength);
    }

    private void setClosedLoopArcLengthWeights(StsList edges)
    {
        int nEdges = edges.getSize();
        if (nEdges == 0) return;

        for (int e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(e);
            StsList edgePoints = edge.getGridEdgePointsList();
            int nEdgePoints = edgePoints.getSize();

            for (int p = 0; p < nEdgePoints; p++)
            {
                StsGridSectionPoint point = (StsGridSectionPoint) edgePoints.getElement(p);
                int row = (int) point.getRowF(null);
                int col = (int) point.getColF(null);
                int rowOrCol = point.getRowOrCol(null);
                if (rowOrCol == ROW)
                    setPointLinkDistance(row, col, rowLinkDistances, StsParameters.largeFloat);
                else
                    setPointLinkDistance(row, col, colLinkDistances, StsParameters.largeFloat);
            }
        }
    }

    public void debugDisplayFaultArcLengths(StsGLPanel3d glPanel3d)
    {
        int row, col;
        StsColor stsColor;

        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        float[] xyz = new float[3];

        try
        {
            if (rowLinkDistances == null) return;

            gl.glDisable(GL.GL_LIGHTING);
            glPanel3d.setViewShift(gl, 4.0 * StsGraphicParameters.gridShift);
            StsSpectrum spectrum = currentModel.getSpectrum("Basic");

            for (row = 0; row < nRows; row++)
            {
                if (rowLinkDistances[row] == null) continue;
                for (col = 0; col < nCols; col++)
                {
                    float arcLength = rowLinkDistances[row][col];
                    if (arcLength > 0.0f)
                    {
                        if (arcLength == StsParameters.largeFloat)
                            stsColor = StsColor.GREY;
                        else
                        {
                            int index = (int) (arcLength / 4);
                            index = index % 8;
                            stsColor = spectrum.getColor(index);
                        }
                        xyz[0] = getXCoor(row, (col + 0.25f));
                        xyz[1] = getYCoor(row, (col + 0.25f));
                        xyz[2] = getZorT(row, col);
                        StsGLDraw.drawPoint(xyz, stsColor, glPanel3d, 4, 0.0);
                    }
                }
            }

            for (row = 0; row < nRows; row++)
            {
                if (colLinkDistances[row] == null) continue;
                for (col = 0; col < nCols; col++)
                {
                    float arcLength = colLinkDistances[row][col];
                    if (arcLength > 0.0f)
                    {
                        int index = (int) (arcLength / 4);
                        index = index % 8;
                        stsColor = spectrum.getColor(index);
                        xyz[0] = getXCoor(row + 0.25f, (float) col);
                        xyz[1] = getYCoor(row + 0.25f, (float) col);
                        xyz[2] = getZorT(row, col);
                        StsGLDraw.drawPoint(xyz, stsColor, glPanel3d, 4, 0.0);
                    }
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("RCGrid.debugDisplayFaultArcLengths() failed.",
                    e, StsException.WARNING);
        }
        finally
        {
            glPanel3d.resetViewShift(gl);
        }
    }


    public void debugSurfacePoint(StsGridPoint gridPoint)
    {
        if (gridPoint == null)
        {
            return;
        }
        StsXYGridable grid = gridPoint.grid;
        if (grid == this)
        {
            for (int n = 0; n < blockGrids.getSize(); n++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
                if (blockGrid.insideGrid(gridPoint.row, gridPoint.col))
                    blockGrid.debugSurfacePoint(gridPoint);
            }
        }
    }

    public boolean deleteModel()
    {
        if (blockGrids != null)
        {
            blockGrids = null;
            blockArray = null;
        }
        deleteTransientArrays();
        setType(MODEL);

        this.changeNullTypes(SURF_GAP_SET, SURF_PNT);
        return true;
    }

    public void resetNullTypes()
    {
        changeNullTypes(StsGridPoint.SURF_GAP, StsGridPoint.SURF_PNT);
    }

    public boolean delete()
    {
        if (!super.delete()) return false;
        deleteSectionEdges();
        return true;
    }

    private void deleteSectionEdges()
    {
        StsObject[] sectionObjects = currentModel.getObjectList(StsSection.class);
        for (StsObject sectionObject : sectionObjects)
        {
            StsSection section = (StsSection) sectionObject;
            section.deleteSectionEdge(this);
        }
    }

    public void deleteTransientArrays()
    {
        faultSurface = null;
        if (blockGrids != null) blockGrids.forEach("deleteTransientArrays");

        rowLinkDistances = null;
        colLinkDistances = null;
        faultDistances = null;
        super.deleteTransientArrays();
    }

    public void deleteDisplayLists(GL gl)
    {
        super.deleteDisplayLists(gl);

        if (blockGrids != null)
        {
            int nBlockGrids = blockGrids.getSize();
            for (int n = 0; n < nBlockGrids; n++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
                StsEdgeLoop edgeLoop = blockGrid.getEdgeLoop();
                if (edgeLoop == null) continue;
                edgeLoop.deleteDisplayLists(gl);
            }
        }
    }

    public void debugDisplayBlockGrids(StsGLPanel3d glPanel3d)
    {
        if (blockArray == null) return;

        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        int nColors = 8;

        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        gl.glDisable(GL.GL_LIGHTING);

        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                StsBlock block = getBlock(i, j);
                if (block != null)
                {
                    int blockIndex = block.getIndex();
                    StsColor color = spectrum.getColor(blockIndex % nColors);
                    StsGLDraw.drawPoint(getXYZorT(i, j), color, glPanel3d, 4, 4, 0.0);
                }
            }
        }
        gl.glEnable(GL.GL_LIGHTING);
    }

    public StsModelSurface getModelSurfaceBelow()
    {
        if (zoneBelow == null) return null;
        return zoneBelow.getBaseModelSurface();
    }


    public void setProperty(StsPropertyVolumeOld pv, int nPropLayer)
    {
        propLayerNumber = nPropLayer;

        propertyVolume = pv;
        propertyName = (pv == null) ? null : pv.getPropertyName() + "_" + nPropLayer;

        deleteSurfDisplayList();
    }

    public StsPropertyVolumeOld getPropertyVolume()
    {
        return propertyVolume;
    }

    public String getPropertyName()
    {
        return propertyName;
    }

    public int getPropertyLayerNumber()
    {
        return propLayerNumber;
    }


    public StsColor getPropColor(int row, int col)
    {
        return propertyVolume.getScaledColor(row, col, propLayerNumber);
    }

    public float getPropValue(int row, int col)
    {
        return propertyVolume.getValue(row, col, propLayerNumber);
    }

    public float getPropValue(float[] xyz)
    {
        return propertyVolume.getValue(xyz[0], xyz[1], propLayerNumber);
    }

    public boolean initialize(StsModel model)
    {
        if (initialized) return true;

        try
        {
            if (originalSurface != null) originalSurface.unload();

            setup();
            setIsVisible(isVisible);
            if (faultPolygons != null)
            {
                for (int i = 0; i < faultPolygons.length; i++) faultPolygons[i].initialize();
                applyFaultPolygons();
            }
            //            if(blockGrids != null) constructLinkedGrid();
            //            if(originalSurface != null) originalSurface.unload();
            initializeSurfaceTextures();
            initialized = true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsModelSurface.classInitialize() failed " +
                    "for surface: " + getName() + ". ", e, StsException.WARNING);
        }
        return true;
    }

    public StsObjectList getSurfaceEdges()
    {
        StsObjectList edges = new StsObjectList(50, 50);

        if (blockGrids == null) return new StsObjectList(1);
        int nBlockGrids = blockGrids.getSize();
        for (int g = 0; g < nBlockGrids; g++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(g);
            StsList groupEdges = blockGrid.getEdges();
            int nEdges = groupEdges.getSize();
            for (int e = 0; e < nEdges; e++)
            {
                StsSurfaceEdge edge = (StsSurfaceEdge) groupEdges.getElement(e);
                edges.add(edge);
            }
        }
        return edges;
    }


    // if display, display position in statusArea; if outsideGridOK, allow for a pick
    // off the grid but find the nearest grid boundary point
    public StsGridPoint getSurfacePosition(StsMousePoint mousePoint, boolean display, boolean outsideGridOK, StsGLPanel3d glPanel3d)
    {
        //        StsPoint point0, point1;
        StsGridPoint nearestGridPoint = null;
        StsBlockGrid nearestBlockGrid = null;
        StsGridPoint gridPoint;
        double nearestGridPointScreenZ = StsParameters.largeDouble;

        boolean debug = false; //StsTrace.getTrace();

        // get projection points and check for valid z range
        double[] point0 = currentModel.getGlPanel3d().getNearProjectedPoint(mousePoint);
        double[] point1 = currentModel.getGlPanel3d().getFarProjectedPoint(mousePoint);
        if (Math.abs(point0[2] - point1[2]) < StsParameters.doubleRoundOff) return null;

        StsGridDefinition gridDefinition = currentModel.getGridDefinition();

        if (type == BUILT && blockGrids != null)
        {
            StsSurfaceEdge currentEdge = StsSurfaceEdge.getCurrentEdge();
            if (currentEdge != null)
            {
                nearestBlockGrid = currentEdge.getBlockGrid();
                if (nearestBlockGrid != null)
                    nearestGridPoint = gridDefinition.getSurfacePosition(point0, point1, nearestBlockGrid, outsideGridOK);
            } else
            {
                int nBlockGrids = blockGrids.getSize();
                for (int g = 0; g < nBlockGrids; g++)
                {
                    StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(g);
                    gridPoint = gridDefinition.getSurfacePosition(point0, point1, blockGrid, outsideGridOK);
                    if (gridPoint == null) continue;
                    if (nearestGridPoint == null)
                    {
                        nearestGridPoint = gridPoint;
                        nearestBlockGrid = blockGrid;
                    } else
                    {
                        double[] gridPointScreenXYZ = currentModel.getGlPanel3d().getScreenCoordinates(gridPoint.getPoint());
                        if (gridPointScreenXYZ[2] < nearestGridPointScreenZ)
                        {
                            nearestGridPoint = gridPoint;
                            nearestBlockGrid = blockGrid;
                            nearestGridPointScreenZ = gridPointScreenXYZ[2];
                        }
                    }
                }
            }
            if (display && nearestGridPoint != null) nearestBlockGrid.displaySurfacePosition(nearestGridPoint);
            return nearestGridPoint;
        } else
        {
            gridPoint = gridDefinition.getSurfacePosition(point0, point1, this, outsideGridOK);
            if (gridPoint != null && display) displaySurfacePosition(gridPoint);
            return gridPoint;
        }
    }


    private void deleteSurfDisplayList()
    {
        if (blockGrids == null) return;
        GL gl = currentModel.getWin3dGL();
        if (gl == null) return;

        int nBlockGrids = blockGrids.getSize();
        for (int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            StsEdgeLoop edgeLoop = blockGrid.getEdgeLoop();
            if (edgeLoop != null) edgeLoop.propertyChanged = true;

        }
    }

    public void display(StsGLPanel glPanel)
    {
        boolean debug;

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
        byte projectZDomain = currentModel.getProject().getZDomain();
        if (projectZDomain != zDomainDisplayed)
        {
            checkDeleteDisplays(gl);
            /*
               if(projectZDomain != zDomainOriginal && isType(BUILT))
               {
                   if(!buildNewDomainModelGeometry(projectZDomain)) return;
               }
            */
            zDomainDisplayed = projectZDomain;
        }
        boolean isDepth = (projectZDomain == StsProject.TD_DEPTH);
        //        boolean isDepth = zDomainDisplayed != zDomainOriginal;
        //        useDisplayLists = currentModel.useDisplayLists;

        StsModelSurfaceClass surfaceClass = getModelSurfaceClass();
        StsGLPanel3d glPanel3d = (StsGLPanel3d) glPanel;
        if (isType(BUILT))
        {
            if (surfaceClass.getDisplayHorizonEdges()) blockGrids.forEach("displayEdges", glPanel3d);

            debug = currentModel.getBooleanProperty("debugDisplayBlockGrids");
            if (debug && faultSurface != null) debugDisplayBlockGrids(glPanel3d);
            String displayMode = StsBuiltModelClass.staticGetDisplayModeString();

            if (blockGrids == null)
            {
                StsException.systemDebug(this, "display", "blockGrids seems to have been nulled asynchronously...debug problem.");
                return;
            }
            int nBlockGrids = blockGrids.getSize();
            for (int n = 0; n < nBlockGrids; n++)
            {
                if (blockGrids == null)
                {
                    StsException.systemDebug(this, "display", "blockGrids seems to have been nulled asynchronously...debug problem.");
                    return;
                }
                StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
                if (!blockGrid.block.isVisible) continue;
                StsColor displayColor = null;
                if (displayMode == StsBuiltModelClass.displayZonesString)
                {
                    if (zoneBelow != null)
                        displayColor = zoneBelow.getStsColor();
                    else
                        displayColor = zoneAbove.getStsColor();
                } else if (displayMode == StsBuiltModelClass.displayBlocksString)
                {
                    StsZoneBlock zoneBlockBelow = getZoneBlockBelow(blockGrid);
                    if (zoneBlockBelow != null)
                        displayColor = zoneBlockBelow.getStsColor();
                    else
                        displayColor = getZoneBlockAbove(blockGrid).getStsColor();
                } else // user layer color if displayLayers or properties don't exist (this color will be used instead)
                {
                    if (zoneBelow != null)
                    {
                        int nLayer = zoneBelow.getTopLayerNumber();
                        displayColor = getLayerColor(nLayer);
                    } else
                        displayColor = zoneAbove.getStsColor();
                }

                boolean displayProperties = (displayMode == StsBuiltModelClass.staticGetDisplayModeString());
                blockGrid.display(currentModel, glPanel3d, displayFill, displayGrid, displayMode, displayProperties, displayColor);
            }

            debug = currentModel.getBooleanProperty("Grid Weights");
            if (debug) debugDisplayFaultArcLengths(glPanel3d);

            if (linkedGrid != null) linkedGrid.debugDisplay(currentModel, glPanel3d);
        } else if (isType(GAPPED))
        {
            if (surfaceClass.getDisplayHorizonEdges()) blockGrids.forEach("displayEdges", glPanel3d);
            displaySurface(glPanel3d);
        } else
            displaySurface(glPanel3d);
    }

    private StsZoneBlock getZoneBlockBelow(StsBlockGrid blockGrid)
    {
        StsZoneBlock[] zoneBlocks = getZoneBlocks(blockGrid);
        if (zoneBlocks == null) return null;
        for (StsZoneBlock zoneBlock : zoneBlocks)
        {
            if (zoneBlock.topGrid == blockGrid)
                return zoneBlock;
        }
        return null;
    }


    private StsZoneBlock getZoneBlockAbove(StsBlockGrid blockGrid)
    {
        StsZoneBlock[] zoneBlocks = getZoneBlocks(blockGrid);
        for (StsZoneBlock zoneBlock : zoneBlocks)
        {
            if (zoneBlock.botGrid == blockGrid)
                return zoneBlock;
        }
        return null;
    }

    private StsZoneBlock[] getZoneBlocks(StsBlockGrid blockGrid)
    {
        StsBlock block = blockGrid.block;
        return blockGrid.block.getZoneBlocks();
    }

    private StsColor getLayerColor(int layer)
    {
        int colorIndex = 0;
        if (layer > 0)
            colorIndex = layer % 32;
        return StsColor.colors32[colorIndex];
    }

    protected StsColor getPointTypeColor(byte pointType)
    {
        return StsParameters.getGapPointTypeColor(pointType);
    }

    /*
        private boolean buildNewDomainModelGeometry(byte zDomain)
        {
            StsSeismicVelocityModel velocityModel = currentModel.project.getSeismicVelocityModel();
            int nBlockGrids = blockGrids.getSize();
            for(int n = 0; n < nBlockGrids; n++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
                blockGrid.buildNewDomainGeometry(zDomain, velocityModel);
            }
            return true;
        }
    */
    private void drawBlockGridsCursorEdge(StsGLPanel3d glPanel3d, StsPoint point0, StsPoint point1)
    {
        boolean drawGaps = false;

        if (blockGrids == null)
        {
            StsException.outputException(new StsException(StsException.WARNING,
                    "StsModelSurface.drawBlockGridsCursorEdge() failed.",
                    "blockGrids is null."));
            return;
        }
    }

    /**
     * We will build a temporary linkedGrid for the gapping operations;
     * a second linkedGrid will be constructed when the fault edges have
     * been computed.
     */
    public boolean gapSurfaceGrid(StsBuiltModel builtModel, StsStatusUI status, float progress, float incProgress)
    {
        try
        {

            changeNullTypes(StsGridPoint.SURF_GAP, StsGridPoint.SURF_PNT);
            // changeNullTypes(StsGridPoint.SURF_GAP_SET, StsGridPoint.SURF_GAP);
            if (!isType(INITIALIZED))
            {
                status.setText("Initializing block voxelGrids for surface: " + getName());
                if (!initializeBlockGrids(builtModel)) return false;
                status.setProgress(progress += 0.3f * incProgress);
                //                StsStatusArea.staticSetProgress(progress += 0.3f*incProgress);
                /*
                                logMessage("Constructing initial section intersections for surface: " + name);
                                if(!intersectSections()) return false;
                                StsStatusArea.staticSetProgress(progress += 0.3f*incProgress);
                */
                status.setText("Constructing linked grid for surface: " + getName());
                constructLinkedGrid();
                status.setProgress(progress += 0.2f * incProgress);
                //                StsStatusArea.staticSetProgress(progress += 0.2f*incProgress);
            }
            setType(INITIALIZED);

            status.setText("Gapping grid for surface: " + getName());
            faultSurface = new StsFaultSurface(currentModel, this, linkedGrid);
            status.setProgress(progress += 0.5f * incProgress);
            //            StsStatusArea.staticSetProgress(progress += 0.5f*incProgress);

            setType(GAPPED);

            return true;
        }
        catch (Exception e)
        {
            status.setText("gapping the surface grid failed for surface: " + getName());
            StsException.outputException("StsModelSurface.GapSurfaceGrid() failed.",
                    e, StsException.WARNING);
            status.setProgress(progress += incProgress);
            //            StsStatusArea.staticSetProgress(progress += incProgress);

            return false;
        }
    }

    /**
     * blockGrids and extrapolate/intersect against sections
     */
    public boolean constructBlockGrids(StsStatusUI status, float progress, float incProgress)
    {

        if (blockGrids == null || blockGrids.getSize() == 0)
        {
            status.setProgress(progress += incProgress);
            //            StsStatusArea.staticSetProgress(progress += incProgress);
            return false;
        }

        //        if(!computeBlockBoundingRectangles()) return false;

        status.setText("Gapping individual block voxelGrids for surface: " + getName());
        if (!gapBlockGrids(status, progress, 0.5f * incProgress)) return false;
        status.setProgress(progress += 0.35f * incProgress);
        //        StsStatusArea.staticSetProgress(progress += 0.35f*incProgress);

        status.setText("Compute block voxelGrids for surface: " + getName());
        constructBlockGrids1(status, progress, 0.4f * incProgress);
        status.setProgress(progress += 0.35f * incProgress);
        //        StsStatusArea.staticSetProgress(progress += 0.35f*incProgress);

        status.setText("Construct linkedGrid for surface: " + getName());
        constructLinkedGrid();
        status.setProgress(progress += 0.05f * incProgress);
        //        StsStatusArea.staticSetProgress(progress += 0.05*incProgress);

        status.setText("Complete blockGrids for surface: " + getName());
        blockGrids.forEach("fillInsidePoints");
        blockGrids.forEach("adjustAuxiliaryEdges");
        progress += 0.05 * incProgress;

        status.setText("Construct block voxelGrids for surface: " + getName());
        blockGrids.forEach("construct");
        progress += 0.2 * incProgress;

        type = BUILT;
        return true;
    }

    private boolean gapBlockGrids(StsStatusUI status, float progress, float incProgress)
    {
        initializeBlockGridArray();
        // initializePointArrays();
        blockGrids.forEach("initializeArrays");
        if (!faultSurface.gapBlockGrids()) return false;
        status.setProgress(progress += incProgress);
        return true;
    }

    private void constructBlockGrids1(StsStatusUI status, float progress, float incProgress)
    {
        int nBlockGrids = blockGrids.getSize();
        incProgress = incProgress / nBlockGrids;
        for (int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            status.setText("Constructing gap grid " + n + " for surface: " + getName());
            blockGrid.construct2dGrid(currentModel);
            status.setProgress(progress += incProgress);
        }
    }

    /**
     * Construct edges along both sides of fault sections and right
     * side of boundary sections.  Build edges between each pair of
     * wells along the side. Edges are connected at vertices so that
     * we can subsequently walk edge loops and form edge groups.
     */
    private boolean initializeBlockGrids(StsBuiltModel builtModel)
    {
        boolean constructionOK = true;
        try
        {
            StsBlock[] blocks = builtModel.getBlocks();
            int nBlocks = blocks.length;
            blockGrids = new StsObjectList(nBlocks, 1, StsBlockGrid.class);

            for (StsBlock block : blocks)
            {
                StsBlockGrid blockGrid = new StsBlockGrid(this, block);
                if (blockGrid.constructInitialEdges())
                {
                    blockGrids.add(blockGrid);
                    block.addBlockGrid(blockGrid);
                } else constructionOK = false;
            }
            return constructionOK;
        }
        catch (Exception e)
        {
            StsException.outputException("StsModelSurface.constructBlockGrids() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }

    private void constructLinkedGrid()
    {
        linkedGrid = new StsEdgeLoopRadialLinkGrid(this);
        for (int n = 0; n < blockGrids.getSize(); n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            blockGrid.addToLinkedGrid(linkedGrid);
        }
        // blockGrids.forEach("addToLinkedGrid", linkedGrid);
        blockGrids.forEach("constructGridBoundingBox");
        linkedGrid.resetMinMaxLimits();
        linkedGrid.orderLinks();
    }

    /*
        protected void drawSurfaceFill(GL gl)
        {
            try
            {
                if (hasNulls && !fillGaps)
                {
                    computeTStripNormals();

                    if(surfaceDisplayable == null)
                    {
                        drawSurfaceFillWithNulls(gl);
                        return;
                    }
                    byte displayType = surfaceDisplayable.getDisplayType();
                    if(displayType == StsSurfaceDisplayable.STRAT_CELL_CENTERED)
                        drawSurfaceFillWithNullsCellCentered(gl, surfaceDisplayable);
                    else if(displayType == StsSurfaceDisplayable.STRAT_APPROX_CELL_CENTERED)
                        drawSurfaceFillWithNullsApproxCellCentered(gl, surfaceDisplayable);
                    else if(displayType == StsSurfaceDisplayable.ORTHO_GRID_CENTERED)
                        drawSurfaceFillWithNullsGridCentered(gl, surfaceDisplayable);
                }
                 else
                 {
                    computeNormals();

                    if(surfaceDisplayable == null)
                    {
                        drawSurfaceFillWithNulls(gl);
                        return;
                    }
                    byte displayType = surfaceDisplayable.getDisplayType();
                    if(displayType == StsSurfaceDisplayable.STRAT_CELL_CENTERED)
                        drawSurfaceFillWithoutNullsCellCentered(gl, surfaceDisplayable);
                    else if(displayType == StsSurfaceDisplayable.STRAT_APPROX_CELL_CENTERED)
                        drawSurfaceFillWithoutNullsCellCentered(gl, surfaceDisplayable);
                    else if(displayType == StsSurfaceDisplayable.ORTHO_GRID_CENTERED)
                        drawSurfaceFillWithoutNullsGridCentered(gl, surfaceDisplayable);
                 }
            }
            catch(Exception e)
            {
    //            StsException.outputException("StsGrid.drawSurface() failed.",
    //                e, StsException.WARNING);
            }
        }

        private void drawSurfaceFillWithoutNullsCellCentered(GL gl, StsSurfaceDisplayable surfaceDisplayable)
        {
            int i, j;

            float[] point = new float[3];

            double diagonalX = -rowXInc + colXInc;
            double diagonalY = -rowYInc + colYInc;

            float rowStartX = xMin;
            float rowStartY = yMin;

            for (i=0; i<nRows-1; i++)
            {
                point[0] = rowStartX;
                point[1] = rowStartY;

                gl.glBegin(GL.GL_TRIANGLE_STRIP);

                for (j=0; j<nCols; j++)
                {
                    point[2] = pointsZ[i][j];
                    StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));
                    gl.glNormal3fv(point);
                    gl.glVertex3fv(normals[i][j]);

                    point[0] += (float)rowXInc;
                    point[1] += (float)rowYInc;
                    point[2] = pointsZ[i+1][j];

                    gl.glNormal3fv(point);
                    gl.glVertex3fv(normals[i+1][j]);

                    point[0] += (float)diagonalX;
                    point[1] += (float)diagonalY;
                }

                gl.glEnd();

                rowStartX += rowXInc;
                rowStartY += rowYInc;
            }
        }

        private void drawSurfaceFillWithNullsCellCentered(GL gl, StsSurfaceDisplayable surfaceDisplayable)
        {
            float startX = xMin;
            float startY = yMin;

            TriangleStrip t = null;
            float[] point = new float[3];
            int i, j = 0;
            int n;

            try
            {
                int nTStrips = tStrips.getSize();
                for (int nt = 0; nt < nTStrips; nt++)
                {
                    t = (TriangleStrip)tStrips.getElement(nt);
                    i = t.rowOrColNumber;
                    j = t.firstIndex;

                    gl.glBegin(GL.GL_TRIANGLES);

                    // first side
                    if (t.firstSide == STRIP_BOT)
                    {
                        point[0] = startX + (j-1)*xInc;
                        point[1] = startY + i*yInc;
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j-1, propLayerNumber));

                        point[2] = pointsZ[i][j-1];
                        gl.glNormal3fv(normals[i][j-1]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[1] -= yInc;
                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);
                    }
                    else if (t.firstSide == STRIP_TOP)
                    {
                        point[0] = startX + (j-1)*xInc;
                        point[1] = startY + i*yInc;
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j-1, propLayerNumber));

                        point[1] += yInc;
                        n = index(i+1,j-1);
                        point[2] = pointsZ[i+1][j-1];
                        gl.glNormal3fv(normals[i+1][j-1]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[1] -= yInc;
                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);
                    }
                    else
                    {
                        point[0] = startX + j*xInc;
                        point[1] = startY + i*yInc;
                    }

                    // middle sides
                    for (; j<t.lastIndex; j++)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);

                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[2] = pointsZ[i+1][j+1];
                        gl.glNormal3fv(normals[i+1][j+1]);
                        gl.glVertex3fv(point);

                        point[0] -= xInc;
                        point[1] -= yInc;
                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j+1];
                        gl.glNormal3fv(normals[i+1][j+1]);
                        gl.glVertex3fv(point);

                        point[1] -= yInc;
                        point[2] = pointsZ[i][j+1];
                        gl.glNormal3fv(normals[i][j+1]);
                        gl.glVertex3fv(point);
                    }

                    // last side: j has been bumped to next j in for loop above
                    if (t.lastSide == STRIP_BOT)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);

                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[1] -= yInc;
                        n = index(i,j+1);
                        point[2] = pointsZ[i][j+1];
                        gl.glNormal3fv(normals[i][j+1]);
                        gl.glVertex3fv(point);
                    }
                    else if (t.lastSide == STRIP_TOP)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);

                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[2] = pointsZ[i+1][j+1];
                        gl.glNormal3fv(normals[i+1][j+1]);
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
    /**
     * uses fast/fuzzy cell-centered drawing
     */
    /*
        private void drawSurfaceFillWithNullsApproxCellCentered(GL gl, StsSurfaceDisplayable surfaceDisplayable)
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
                        point[0] = startX + (j-1)*xInc;
                        point[1] = startY + i*yInc;
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j-1, propLayerNumber));

                        point[2] = pointsZ[i][j-1];
                        gl.glNormal3fv(normals[i][j-1]);
                        gl.glVertex3fv(point);
                        gl.glNormal3fv(normals[i][j-1]);
                        gl.glVertex3fv(point);
                    }
                    else if (t.firstSide == STRIP_TOP)
                    {
                        point[0] = startX + (j-1)*xInc;
                        point[1] = startY + i*yInc;
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j-1, propLayerNumber));

                        point[2] = pointsZ[i+1][j-1];
                        gl.glNormal3fv(normals[i+1][j-1]);
                        gl.glVertex3fv(point);
                        gl.glNormal3fv(normals[i+1][j-1]);
                        gl.glVertex3fv(point);
                    }
                    else
                    {
                        point[0] = startX + j*xInc;
                        point[1] = startY + i*yInc;
                    }

                    // middle sides
                    for (; j<=t.lastIndex; j++)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);

                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
                        gl.glVertex3fv(point);

                        point[0] += xInc;
                        point[1] -= yInc;
                    }

                    // last side: j has been bumped to next j in for loop above
                    if (t.lastSide == STRIP_BOT)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[2] = pointsZ[i][j];
                        gl.glNormal3fv(normals[i][j]);
                        gl.glVertex3fv(point);
                    }
                    else if (t.lastSide == STRIP_TOP)
                    {
                        StsColor.setGLColor(gl, surfaceDisplayable.getCellColor(i, j, propLayerNumber));

                        point[1] += yInc;
                        point[2] = pointsZ[i+1][j];
                        gl.glNormal3fv(normals[i+1][j]);
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
    public void applyMinCorrelFilter(boolean applyMinCorrelFilter, StsSurfaceTexture correlCoefTexture)
    {
        super.applyMinCorrelFilter(applyMinCorrelFilter, correlCoefTexture);
        interpolateDistanceTransform();
        geometryChanged = true;
    }

    public float[][] initializeDistances()
    {
        if (pointsNull == null) return null;

        float[][] distances = new float[nRows][nCols];
        StsSurfaceTexture correlCoefTexture = originalSurface.getCorrelCoefTexture();
        if (applyMinCorrelFilter && correlCoefTexture != null)
        {
            byte[] data = correlCoefTexture.getTextureData();
            int n = 0;
            for (int row = 0; row < nRows; row++)
                for (int col = 0; col < nCols; col++, n++)
                {
                    pointsZ[row][col] = originalSurface.getZNotNull(row, col);
                    if (pointsZ[row][col] != nullZValue && data[n] != -1)
                    {
                        distances[row][col] = 0.0f;
                        setPointNull(row, col, SURF_PNT);
                    } else
                    {
                        distances[row][col] = largeFloat;
                        setPointNull(row, col, SURF_GAP_SET);
                    }
                }
        } else
        {
            for (int row = 0; row < nRows; row++)
                for (int col = 0; col < nCols; col++)
                    if (getZNotNull(row, col) != nullZValue)
                        distances[row][col] = 0.0f;
                    else
                        distances[row][col] = largeFloat;
        }
        return distances;
    }

    public float distanceTransformInterpolation(int row, int col, StsDistanceTransformPoint[] points, float maxInterpolationDistance)
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
                System.out.println("    row: " + nearestRow + " col: " + nearestCol + " distance: " + distance + " z: " + pointsZ[nearestRow][nearestCol]);
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                z = pointsZ[nearestRow][nearestCol];
            }
        }

        float d0 = points[0].distance;

        if (nearestDistance == largeFloat)
        {
            if (debug) System.out.println("    no available points to interpolate");
            return d0;
        }

        if (d0 != largeFloat)
        {
            float z0 = pointsZ[row][col];
            float w0 = 1.0f / d0;
            float w1 = 1.0f / nearestDistance;
            float f = w0 / (w0 + w1);
            if (debug)
                System.out.print("  interpolation between old z: " + z0 + " and new z: " + z + " interpolation factor: " + f + " nearest distance: " + nearestDistance);
            z = f * z0 + (1 - f) * z;
            if (debug) System.out.println(" resulting z: " + z);
            if (d0 < nearestDistance) nearestDistance = d0;
        } else
        {
            if (debug)
                System.out.println("    new z value at this point: " + z + " nearest distance: " + nearestDistance);
        }
        //        if(wtSum == 0.0f) return largeFloat;

        //        float z;
        //        if(wtSum == 1.0f)
        //            z = zSum;
        //        else
        //            z = zSum/wtSum;

        pointsZ[row][col] = z;
        return nearestDistance;
    }


    public void adjustBlockGridThickness(StsModelSurface baseSurface)
    {

    }

    public void interpolateBlockGrids(StsModelSurface baseSurface, StsModelSurface topSurface, float f)
    {
    }

    public void interpolateSurfaces(StsModelSurface base, StsModelSurface top, float f)
    {
    }
/*
    public final void setPointNull(int row, int col, byte nullType)
    {
        System.out.println(pointsNull[row][col] + " changed to " + nullType);
        super.setPointNull(row, col, nullType);
    }
*/
    public StsFieldBean[] getDisplayFields()
    {
        if (modelDisplayFields == null)
        {
            modelSurfaceTextureListBean = new StsComboBoxFieldBean(StsModelSurface.class, "newSurfaceTexture", "Property", "surfaceTextureList");
            offsetBean = new StsFloatFieldBean(StsModelSurface.class, "offset", true, "Offset", true);
            colorscaleBean = new StsEditableColorscaleFieldBean(StsModelSurface.class, "colorscale");

            modelDisplayFields = new StsFieldBean[]
                    {
                            //    	new StsBooleanFieldBean(StsModelSurface.class, "markerVisible", "Marker"),
                            new StsBooleanFieldBean(StsModelSurface.class, "transparent", "Transparent:"),
                            new StsBooleanFieldBean(StsModelSurface.class, "isVisible", "Surface"),
                            new StsBooleanFieldBean(StsModelSurface.class, "displayGrid", "Grid"),
                            new StsBooleanFieldBean(StsModelSurface.class, "displayFill", "Fill"),
                            new StsBooleanFieldBean(StsModelSurface.class, "displayGaps", "Gaps"),
                            //    	new StsBooleanFieldBean(StsModelSurface.class, "displayDomains", "Blocks"),
                            new StsColorComboBoxFieldBean(StsModelSurface.class, "stsColor", "Color", currentModel.getSpectrum("Basic")),
                            modelSurfaceTextureListBean,
                            offsetBean,
                            colorscaleBean                            
                    };
        }
        // initModelSurfaceTextureList();
        return modelDisplayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if (modelPropertyFields == null)
        {
            modelPropertyFields = new StsFieldBean[]
                    {
                            new StsStringFieldBean(StsModelSurface.class, "markerName", "Marker"),
                            new StsIntFieldBean(StsModelSurface.class, "nCols", false, "nX Cols"),
                            new StsIntFieldBean(StsModelSurface.class, "nRows", false, "nY Rows"),
                            new StsDoubleFieldBean(StsModelSurface.class, "xOrigin", false, "X Origin"),
                            new StsDoubleFieldBean(StsModelSurface.class, "yOrigin", false, "Y Origin"),
                            new StsFloatFieldBean(StsModelSurface.class, "xInc", false, "X Inc"),
                            new StsFloatFieldBean(StsModelSurface.class, "yInc", false, "Y Inc"),
                            new StsFloatFieldBean(StsSurface.class, "xMin", false, "X Loc Min"),
                            new StsFloatFieldBean(StsSurface.class, "yMin", false, "Y Loc Min"),
                            new StsFloatFieldBean(StsModelSurface.class, "xSize", false, "X Size"),
                            new StsFloatFieldBean(StsModelSurface.class, "ySize", false, "Y Size"),
                            new StsFloatFieldBean(StsModelSurface.class, "angle", false, "XY Rot Angle"),
                            new StsFloatFieldBean(StsModelSurface.class, "zMin", false, "Min Depth"),
                            new StsFloatFieldBean(StsModelSurface.class, "zMax", false, "Max Depth"),
                    };
        }
        return modelPropertyFields;
    }

    //    public StsFieldBean[] getDisplayFields() { return modelDisplayFields; }
    //    public StsFieldBean[] getPropertyFields() { return modelPropertyFields; }

    public StsObjectPanel getObjectPanel()
    {
        modelObjectPanel = StsObjectPanel.checkConstruct(modelObjectPanel, this, true);
        return modelObjectPanel;
    }

    public void treeObjectSelected()
    {
        getModelSurfaceClass().selected(this);
    }

    static public StsModelSurfaceClass getModelSurfaceClass()
    {
        return (StsModelSurfaceClass) currentModel.getCreateStsClass(StsModelSurface.class);
    }

    public boolean anyDependencies()
    {
        return false;
    }
}


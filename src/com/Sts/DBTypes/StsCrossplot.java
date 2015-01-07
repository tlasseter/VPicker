package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.text.*;

public class StsCrossplot extends StsMainObject implements StsTreeObjectI, StsVolumeDisplayable, ActionListener, StsViewable
{
    protected StsObjectRefList volumes;
    protected StsColorscale colorscale;
    protected float[][] totalAxisRanges = null;
    protected StsObjectRefList polygons;
    protected StsTypeLibrary polygonTypeLibrary;
    protected boolean cropIt = true;
//    protected StsObjectRefList crossplotPoints = null;
    protected float[] dataHist = new float[255]; // Histogram of the data distribution
    protected float[] attHist = new float[255]; // Histogram of the attribute Data

    transient int currentDirNo = -1;
    transient int currentPlaneIndex = -1;
    transient float currentDirCoordinate = StsParameters.largeFloat;
    transient public byte[] planeDensityData;
    transient public byte[] planeAttributeData;
    transient public byte[] volumeDensityData;
    transient public byte[] volumeAttributeData;
    transient int densityMax = 0;
    transient float maxAtt = Float.MIN_VALUE, minAtt = Float.MAX_VALUE;
//    transient public boolean dataDisplayed = false;
    transient int nRows = 255, nCols = 255;

    transient int colorListNum = 0;
    transient boolean colorListChanged = true;
    transient float[][] arrayRGBA = null;
//    transient int texture = 0;
    boolean backgroundTransparent = true;

    transient boolean runTimer = false;
    transient StsTimer timer;
    transient public StsToolbar toolbar = null;
//    transient private boolean reDraw = false;
    transient private boolean filledPolygons = false;
    transient public byte[] polygonData = null;
    transient public boolean polygonsProcessed = false;
//    transient String titleLabel = null;

    transient private float planeCoor;
    transient byte[][] axisData;
    transient byte[] attData;
    transient int nCursorRows, nCursorCols;
    transient boolean debug = false;
    protected boolean displayDensity = false;

    static protected String defaultLibraryName = new String("crossplot");
    static StsObjectPanel objectPanel = null;

    static String[] directionLabels = new String[] { "Crossline Number: ", "Inline Number: ", "Slice Number: " };

    static public final StsFieldBean[] displayFields = null;

    static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsCrossplot.class, "colorscale");
    static public final StsFieldBean[] propertyFields = new StsFieldBean[]
    {
        new StsStringFieldBean(StsCrossplot.class, "name", true, "Name:"),
        new StsBooleanFieldBean(StsCrossplot.class, "cropIt", "Crop Volumes:"),
        colorscaleBean
    };

    public StsCrossplot()
    {
    }


    public StsCrossplot(boolean persistent)
    {
        super(persistent);
    }

    static public StsCrossplot constructor(String name, StsSeismicVolume[] volumes)
    {
        try
        {
            StsCrossplot crossplot = new StsCrossplot(false);
            if(!crossplot.initialize(name, volumes)) return null;
            return crossplot;
        }
        catch(Exception e)
        {
            StsException.outputException("StsCrossplot.constructor(win3d) failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private boolean initialize(String name, StsSeismicVolume[] volumes)
    {
        return initialize(name, volumes, currentModel.getSpectrum(getCrossplotClass().getXplotSpectrumName()));
    }

    private boolean initialize(String name, StsSeismicVolume[] volumes, StsSpectrum spectrum)
    {
        try
        {
            if(volumes == null || volumes.length < 2)
            {
                new StsMessage(currentModel.win3d, StsMessage.ERROR, "Minimum of 2 seismic volumes must be specified for a crossplot.");
                return false;
            }
            setName(name);
            int nVolumes = volumes.length;
            this.volumes = StsObjectRefList.constructor(nVolumes, 1, "volumes", this);
            this.volumes.add(volumes);
            this.colorscale = new StsColorscale("Crossplot", spectrum , 0.0f, 255.0f);
            if(runTimer) timer = new StsTimer();

            totalAxisRanges = new float[2][2];

            for(int n = 0; n < 2; n++)
            {
                totalAxisRanges[n][0] = volumes[n].dataMin;
                totalAxisRanges[n][1] = volumes[n].dataMax;
            }

            polygonTypeLibrary = StsTypeLibrary.getCreateGenericLibrary();
            if(polygons == null) polygons = StsObjectRefList.constructor(5, 5, "polygons", this);

            initialize(currentModel);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsCrossplot.classInitialize() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static public StsCrossplotClass getCrossplotClass()
    {
        return (StsCrossplotClass)currentModel.getCreateStsClass(StsCrossplot.class);
    }

    public boolean initialize(StsModel model)
    {
        StsCrossplotClass crossplotClass = getCrossplotClass();

        if(crossplotClass.getCurrentCrossplot() == null) return true;
        if(colorscale != null) colorscale.addActionListener(this);
        colorscaleBean.setHistogram(dataHist);

        StsSeismicVolume[] sv = getVolumes();
        for(int i=0; i<sv.length; i++)
            sv[i].addActionListener(this);

        processPolygons();

        return true;
    }    

    public void addToModel()
    {
		super.addToModel();
        currentModel.checkAddToCursor3d(this);
    }

    public boolean hasSameSeismicVolumes(StsSeismicVolume[] otherVolumes)
    {
        int nOtherVolumes = otherVolumes.length;
        if(nOtherVolumes != volumes.getSize()) return false;
        for(int n = 0; n < nOtherVolumes; n++)
            if(otherVolumes[n] != volumes.getElement(n)) return false;
        return true;
    }
    
    public boolean hasAttributeData()
    {
    	if(volumes == null) return false;
    	if(volumes.getSize() > 2)
    		return true;
    	return false;
    }
/*
    public void itemStateChanged(ItemEvent e)
    {
        StsMessageFiles.infoMessage("Crossplot itemStateChanged Performed...");
        if(e.getItem() instanceof StsColorscale)
        {
            updateCrossplotColors();
        }
        else if(!(e.getItem() instanceof StsSeismicVolume)) return;
    }
*/
	public void actionPerformed(ActionEvent e)
	{
//		StsMessageFiles.infoMessage("Crossplot actionPerformed called...");
		if(e.getSource() instanceof StsColorscale)
		{
			updateCrossplotColors();
		}
	}

    private void updateCrossplotColors()
    {
        colorListChanged = true;
		if (currentModel == null) return;
        if (currentModel.viewPersistManager.families == null) return;

        for(int i = 0; i<currentModel.viewPersistManager.families.length; i++)
        {
            StsWin3dBase[] windows = currentModel.getWindows(i);
            for (int w = 0; w < windows.length; w++)
            {
                StsWin3dBase window = windows[w];
                StsView[] views = window.getDisplayedViews();
                for(int v = 0; v < views.length; v++)
                    if (views[v] instanceof StsViewXP)
                        ((StsViewXP)views[v]).checkClearTexture(this);
            }
        }

//        currentModel.win3dDisplay();
    }

    public void clearCrossplotTextureDisplays()
    {
        clearData();
        currentModel.toggleOnCursor3dObject(this);
        StsComboBoxToolbar toolbar = (StsComboBoxToolbar)currentModel.win3d.getToolbarNamed(StsComboBoxToolbar.NAME);
        if(toolbar != null) toolbar.comboBoxSetItem(this);
    }

    public void subVolumeVisibilityChanged(Boolean isVisibleObject)
    {
        clearData();
    }

    public void clearData()
    {
        planeDensityData = null;
        volumeDensityData = null;
    }

    public StsObjectRefList getPolygons() { return polygons; }

    public StsObject[] getPolygonArray()
    {
        if(polygons == null) return new StsObject[0];
        else return polygons.getElements();
    }

    public StsSeismicVolume[] getVolumes()
    {
        if(volumes == null) return null;
        return (StsSeismicVolume[]) volumes.getCastList(StsSeismicVolume.class);
    }

    public boolean getIsVisibleOnCursor()
    {
        return isVisible && getCrossplotClass().getIsVisibleOnCursor();
    }

    public boolean delete()
    {
		currentModel.win3d.displayPreviousObject(StsCrossplot.class);
		super.delete();
		return true;
	/*
        StsCrossplotClass volumeClass = getCrossplotClass();
        StsCrossplot currentCrossplot = volumeClass.getCurrentCrossplot();
        super.delete(); // sets this.index to -1
        if(currentCrossplot == this)
        {
            int nVolumes = volumeClass.getSize();
            if(nVolumes > 0)
                currentCrossplot = (StsCrossplot)volumeClass.getElement(nVolumes-1);
            else
            {
                currentCrossplot = null;
                // Remove XP Button
                currentModel.validateToolbarStates();
            }

            volumeClass.setCurrentObject(currentCrossplot);
        }
        return true;
    */
    }

    public boolean displayOK()
    {
        if(volumes.getSize() < 2)
            return false;
        if(polygons == null || polygons.getSize() == 0)
           return false;
        return true;
    }


    public byte[] getDensityData(int dirNo, float dirCoordinate, StsCursor3d cursor3d)
    {
        if(getCrossplotClass().getDisplayEntireVolume())
        {
            if(checkLoadVolumeData(cursor3d))
            {
                colorscale.setRange(0.0f, densityMax);
                currentModel.clearDisplayTextured3dCursors(this);
                return volumeDensityData;
            }
            else
                return null;
        }
        else
        {
            if(checkLoadPlaneData(dirNo, dirCoordinate, cursor3d))
            {
                colorscale.setRange(0.0f, densityMax);
                currentModel.clearDisplayTextured3dCursors(this);
                return planeDensityData;
            }
            else
                return null;
        }
    }

    public byte[] getAttributeData(int dirNo, float dirCoordinate, StsCursor3d cursor3d)
    {
        minAtt = Integer.MAX_VALUE;
        maxAtt = Integer.MIN_VALUE;
        if(getCrossplotClass().getDisplayEntireVolume())
        {
            if(checkLoadVolumeData(cursor3d))
            {
                colorscale.setRange(minAtt, maxAtt);
                currentModel.clearDisplayTextured3dCursors(this);
                return volumeAttributeData;
            }
            else
                return null;
        }
        else
        {
            if(checkLoadPlaneData(dirNo, dirCoordinate, cursor3d))
            {
                colorscale.setRange(minAtt, maxAtt);
                currentModel.clearDisplayTextured3dCursors(this);
                return planeAttributeData;
            }
            else
                return null;
        }
    }

    private boolean checkLoadPlaneData(int dirNo, float dirCoordinate, StsCursor3d cursor3d)
    {
        if(debug) System.out.println("currentDirNo: " + currentDirNo + " newDirNo: " + dirNo +
                                     "currentDirCoordinate: " + currentDirCoordinate + " newDirCoordinate: " + dirCoordinate);
        if(planeDensityData != null && currentDirNo == dirNo && currentDirCoordinate == dirCoordinate)
            return true;

        currentDirNo = dirNo;
        currentDirCoordinate = dirCoordinate;

//        dataDisplayed = false;

        // Is the Cursor inside of the Crop Limits
        StsCropVolume cropVolume = currentModel.getProject().getCropVolume();
        if(cropIt && cropVolume.isDirCoordinateCropped(currentDirNo, currentDirCoordinate))
        {
//            planeDensityData = new byte[nRows*nCols];
            return false;
        }

        if(volumes == null)
            return false;
        int nVolumes = volumes.getSize();
        if(nVolumes < 2)
            return false;
        StsSeismicVolume volume = null;
        StsSeismicVolume attVolume = null;

        axisData = new byte[2][];
        boolean[][] skip = new boolean[2][];
        for(int n = 0; n < 2; n++)
        {
            volume = (StsSeismicVolume)volumes.getElement(n);
            axisData[n] = volume.readBytePlaneData(currentDirNo, dirCoordinate);

            if(axisData[n] == null) return false;

            // Set the Combined Crop Limits
//            if(cropIt) setCropLimits(volume);
        }

        if(nVolumes > 2)
        {
            attVolume = (StsSeismicVolume)volumes.getElement(2);
            attData = volume.readBytePlaneData(currentDirNo, dirCoordinate);
        }
        if(debug) System.out.println("checkLoadData axisData[0][0]: " + axisData[0][0] + " axisData[1][0]: " + axisData[1][0]);

        // compute pointsPerBin (255 x 255 bins)
        int[][] pointsPerBin = new int[nRows][nCols];
        float[][] attribute = new float[nRows][nCols];

        int cursorRowMin = cropVolume.getCursorRowMin(currentDirNo);
        int cursorRowMax = cropVolume.getCursorRowMax(currentDirNo);
        int cursorColMin = cropVolume.getCursorColMin(currentDirNo);
        int cursorColMax = cropVolume.getCursorColMax(currentDirNo);
        int nCursorCols = currentModel.getProject().getRotatedBoundingBox().getNCursorCols(currentDirNo);

        int volRowMin = volume.getSubVolumeCursorRowMin(cropVolume, currentDirNo);
        int volRowMax = volume.getSubVolumeCursorRowMax(cropVolume, currentDirNo);
        int volColMin = volume.getSubVolumeCursorColMin(cropVolume, currentDirNo);
        int volColMax = volume.getSubVolumeCursorColMax(cropVolume, currentDirNo);
        int nVolumeCols = volume.getNCursorCols(currentDirNo);

        byte[] subVolumePlane = null;
        boolean subVolumeApplied = isSubVolumeApplied();
        if(subVolumeApplied)
        {
            if(debug) System.out.println("Crossplot subvolume applied.");
            if(!cursor3d.hasSubVolumes()) subVolumeApplied = false;
        }
        if(volRowMin >= 0 && volRowMax >= 0 && volColMin >= 0 && volColMax >= 0)
        {
            int cursorRow = cursorRowMin;
            for (int volRow = volRowMin; volRow <= volRowMax; volRow++, cursorRow++)
            {
                int cursorIndex = cursorRow*nCursorCols + cursorColMin;
                int volumeIndex = volRow*nVolumeCols + volColMin;
                for (int volCol = volColMin; volCol <= volColMax; volCol++, cursorIndex++, volumeIndex++)
                {
                    if (subVolumeApplied && subVolumePlane[cursorIndex] == 0)continue;
                    int dataX = StsMath.signedByteToUnsignedInt(axisData[0][volumeIndex]);
                    int dataY = StsMath.signedByteToUnsignedInt(axisData[1][volumeIndex]);
                    float attValue = 0.0f;
                    if(nVolumes > 2)
                        attValue = attVolume.getScaledValue(attData[volumeIndex]);

                    if (dataX < 255 && dataY < 255) // 255 value represents null
                    {
                        pointsPerBin[dataY][dataX]++;
                        attribute[dataY][dataX] += attValue;
                        densityMax = Math.max(densityMax, pointsPerBin[dataY][dataX]);
                    }
                }
            }
        }
        // Calculate the average attribute
        if(nVolumes > 2)
        {
            for (int n = 0; n < nRows; n++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if(pointsPerBin[n][j] == 0)
                        attribute[n][j] = 0.0f;
                    else
                    {
                        attribute[n][j] = attribute[n][j] / (float)pointsPerBin[n][j];
                    }
                    maxAtt = Math.max(maxAtt, attribute[n][j]);
                    minAtt = Math.min(minAtt, attribute[n][j]);
                }
            }
        }
        // byte 255 (signed byte -1) is the transparent color
        if(planeDensityData == null)
            planeDensityData = new byte[nRows*nCols];
        for(int n = 0; n < nRows*nCols; n++)
            planeDensityData[n] = (byte)255;

        if(nVolumes > 2)
        {
            if(planeAttributeData == null)
                planeAttributeData = new byte[nRows * nCols];
            for (int n = 0; n < nRows * nCols; n++)
                planeAttributeData[n] = (byte) 255;
        }
        // scale to 0 to 254
        float scale = 254.0f/densityMax;
        // hack so we have a reasonable amount of stuff to look at
        scale = 3.0f*scale;

        int n = 0;
        int idx = 0;
        for(int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++, n++)
            {
                int density = 0;
                int attValue = 0;
                if(pointsPerBin[i][j] > 0)
                {
                    density = StsMath.ceiling(pointsPerBin[i][j]*scale);
                    planeDensityData[n] = unsignedIntToSignedByte254(density);
                    if(nVolumes > 2)
                    {
                        attValue = (int) (254 * (attribute[i][j] - attVolume.getDataMin()) / (attVolume.getDataMax() - attVolume.getDataMin()));
                        planeAttributeData[n] = (byte)attValue;
                        if(attValue > 254) attValue = 254;
                    }
                }
                if(density > 254) density = 254;
                dataHist[density]++;
                attHist[attValue]++;
            }
        }

        int numValidValues = (nRows*nCols) - (int)dataHist[0];
        dataHist[0] = 0;
        attHist[0] = 0;
        for(int i=0; i<255; i++)
        {
            dataHist[i] = (float) ((float) dataHist[i] / ((float) numValidValues) * 100.0f);
            if (nVolumes > 2)
                attHist[i] = (float) ((float) attHist[i] / ((float) numValidValues) * 100.0f);
        }
        return true;
    }

    private boolean isSubVolumeApplied()
    {
        // Need to verify that subvolumes are activated and that the crossplot wants to
        // filter by them.
        boolean svApplied = ((StsSubVolumeClass)currentModel.getStsClass(StsSubVolume.class)).getIsApplied();
        boolean filterCrossplotBySubVolume = getCrossplotClass().getDisplayOnSubVolumes();
        if(svApplied && filterCrossplotBySubVolume)
            return true;
        else
            return false;
//        return currentModel.getClassBooleanBeanValue("com.Sts.DBTypes.StsSubVolume", "isApplied", false);
    }

    private boolean checkLoadVolumeData(StsCursor3d cursor3d)
    {
        if(volumeDensityData != null)
            return true;

        if(volumes == null)
            return false;

        int nVolumes = volumes.getSize();
        if(nVolumes < 2)
            return false;

        StsSeismicVolume volumeX = (StsSeismicVolume)volumes.getElement(0);
        StsSeismicVolume volumeY = (StsSeismicVolume)volumes.getElement(1);
        StsSeismicVolume attVolume = null;
        float[][] attribute = null;
        if(nVolumes > 2)
        {
            attVolume = (StsSeismicVolume) volumes.getElement(2);
            attribute = new float[nRows][nCols];
        }
        // compute pointsPerBin (255 x 255 bins)
        int[][] pointsPerBin = new int[nRows][nCols];


        StsProject project = currentModel.getProject();
        StsRotatedGridBoundingBox projectBoundingBox = project.getRotatedBoundingBox();
        StsCropVolume cropVolume = project.getCropVolume();

        StsRotatedGridBoundingBox boundingBox = new StsRotatedGridBoundingBox(false);
        if(!cropVolume.isCropped())
            boundingBox.initialize(projectBoundingBox);
        else
            boundingBox.initialize(cropVolume);

        boolean subVolumeApplied = isSubVolumeApplied();
        if(debug) System.out.println("subVolumeApplied: " + subVolumeApplied);
        StsGridBoundingBox subVolumeBoundingBox = null;
        int nSubVolumeCols = -1;
        if(subVolumeApplied)
        {
            StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)currentModel.getStsClass(StsSubVolume.class);
            subVolumeBoundingBox = subVolumeClass.getBoundingBox();
            if(subVolumeBoundingBox != null)
            {
                boundingBox.intersectBoundingBox(subVolumeBoundingBox);
                nSubVolumeCols = subVolumeBoundingBox.getNCols();
            }
            else
                subVolumeApplied = false;
        }

//        volumeX.setupSliceBlocks();
//        volumeY.setupSliceBlocks();

        ByteBuffer sliceXBuffer, sliceYBuffer, sliceAttBuffer = null;

        int rowMin = projectBoundingBox.getNearestBoundedRowCoor(boundingBox.yMin);
        int rowMax = projectBoundingBox.getNearestBoundedRowCoor(boundingBox.yMax);
        int colMin = projectBoundingBox.getNearestBoundedColCoor(boundingBox.xMin);
        int colMax = projectBoundingBox.getNearestBoundedColCoor(boundingBox.xMax);
        int sliceMin = projectBoundingBox.getNearestSliceCoor(boundingBox.zMin);
        int sliceMax = projectBoundingBox.getNearestSliceCoor(boundingBox.zMax);
        int nTotalCols = projectBoundingBox.nCols;

        float z = projectBoundingBox.getZCoor(sliceMin);
        float zInc = projectBoundingBox.zInc;

        if(debug) System.out.println("Computing volume crossplot density data for volume defined by bounding box: " + boundingBox.toDetailString());

        if(subVolumeApplied)
        {
            byte[] subVolumePlane = null;

            for (int slice = sliceMin; slice <= sliceMax; slice++, z += zInc)
            {
                int volumeSlice = getVolumeSlice(volumeX, z);
                if(volumeSlice == -1)
                    continue;
                sliceXBuffer = volumeX.getSliceBuffer(volumeSlice);
                if (sliceXBuffer == null)
                    continue;
                sliceYBuffer = volumeY.getSliceBuffer(volumeSlice);
                if (sliceYBuffer == null)
                    continue;
                if(attVolume != null)
                    sliceAttBuffer = attVolume.getSliceBuffer(volumeSlice);

                subVolumePlane = cursor3d.getSubVolumePlane(StsCursor3d.ZDIR, z, attVolume.getZDomain());
                if (subVolumePlane != null)
                {
                    for (int row = rowMin; row <= rowMax; row++)
                    {
                        volumeX.setSliceBufferPosition(sliceXBuffer, row, colMin);
                        volumeY.setSliceBufferPosition(sliceYBuffer, row, colMin);
                        if(attVolume != null)
                            attVolume.setSliceBufferPosition(sliceAttBuffer, row, colMin);

                        int n = row*nTotalCols + colMin;
                        for (int col = colMin; col <= colMax; col++, n++)
                        {
                            if (subVolumePlane[n] != 0)
                            {
                                int dataX = StsMath.signedByteToUnsignedInt(sliceXBuffer.get());
                                int dataY = StsMath.signedByteToUnsignedInt(sliceYBuffer.get());
                                float attValue = 0.0f;
                                if(nVolumes > 2)
                                    attValue = attVolume.getScaledValue(sliceAttBuffer.get());

                                if (dataX < 255 && dataY < 255) // 255 value represents null
                                {
                                    if(attVolume != null)
                                        attribute[dataY][dataX] += attValue;
                                    pointsPerBin[dataY][dataX]++;
                                    densityMax = Math.max(densityMax, pointsPerBin[dataY][dataX]);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (int slice = sliceMin; slice <= sliceMax; slice++, z += zInc)
            {
                int volumeSlice = getVolumeSlice(volumeX, z);
                if(volumeSlice == -1)
                    continue;
                sliceXBuffer = volumeX.getSliceBuffer(volumeSlice);
                if (sliceXBuffer == null)
                    continue;
                sliceYBuffer = volumeY.getSliceBuffer(volumeSlice);
                if (sliceYBuffer == null)
                    continue;
                if(attVolume != null)
                    sliceAttBuffer = attVolume.getSliceBuffer(volumeSlice);

                for (int row = rowMin; row <= rowMax; row++)
                {
                    volumeX.setSliceBufferPosition(sliceXBuffer, row, colMin);
                    volumeY.setSliceBufferPosition(sliceYBuffer, row, colMin);
                    if(attVolume != null)
                            attVolume.setSliceBufferPosition(sliceAttBuffer, row, colMin);
                    for (int col = colMin; col <= colMax; col++)
                    {
                        int dataX = StsMath.signedByteToUnsignedInt(sliceXBuffer.get());
                        int dataY = StsMath.signedByteToUnsignedInt(sliceYBuffer.get());
                        float attValue = 0.0f;
                        if(nVolumes > 2)
                            attValue = attVolume.getScaledValue(sliceAttBuffer.get());

                        if (dataX < 255 && dataY < 255) // 255 value represents null
                        {
                            if(attVolume != null)
                                attribute[dataY][dataX] += attValue;

                            pointsPerBin[dataY][dataX]++;
                            densityMax = Math.max(densityMax, pointsPerBin[dataY][dataX]);
                        }
                    }
                }
            }
        }

        // Calculate the average attribute
        if(attVolume != null)
        {
            for (int n = 0; n < 255; n++)
            {
                for (int j = 0; j < 255; j++)
                {
                    if(pointsPerBin[n][j] == 0)
                       attribute[n][j] = 0.0f;
                   else
                       attribute[n][j] = attribute[n][j]/(float)pointsPerBin[n][j];
                    maxAtt = Math.max(maxAtt, attribute[n][j]);
                    minAtt = Math.min(minAtt, attribute[n][j]);
                }
            }
        }

        // byte 255 (signed byte -1) is the transparent color
        if(volumeDensityData == null)
            volumeDensityData = new byte[nRows*nCols];
        for(int n = 0; n < nRows*nCols; n++)
            volumeDensityData[n] = (byte) 255;

        if(nVolumes > 2)
        {
            volumeAttributeData = new byte[nRows * nCols];
            for (int n = 0; n < nRows * nCols; n++)
                volumeAttributeData[n] = (byte) 255;
        }

        // scale to 0 to 254
        float scale = 254.0f/densityMax;
        // hack so we have a reasonable amount of stuff to look at
        scale = 3.0f*scale;

        int n = 0;
        int idx = 0;
        for(int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++, n++)
            {
                int density = 0;
                int attValue = 0;
                if(pointsPerBin[i][j] > 0)
                {
                    density = StsMath.ceiling(pointsPerBin[i][j]*scale);
                    volumeDensityData[n] = unsignedIntToSignedByte254(density);
                    if(nVolumes > 2)
                    {
                        attValue = (int) (254 * (attribute[i][j] - attVolume.getDataMin()) / (attVolume.getDataMax() - attVolume.getDataMin()));
                        volumeAttributeData[n] = (byte)attValue;
                        if(attValue > 254) attValue = 254;
                    }
                    if(density > 254) density = 254;
                }
//                densityData[n] = unsignedIntToSignedByte254((int)(pointsPerBin[i][j] * scale));
//                idx = StsMath.signedByteToUnsignedInt(densityData[n++]);
                dataHist[density]++;
                attHist[attValue]++;
            }
        }

        int numValidValues = (nRows*nCols) - (int)dataHist[0];
        dataHist[0] = 0;
        attHist[0] = 0;
        for(int i=0; i<255; i++)
        {
            dataHist[i] = (float) ((float) dataHist[i] / ((float) numValidValues) * 100.0f);
            if(nVolumes > 2)
                attHist[i] = (float) ((float) attHist[i] / ((float) numValidValues) * 100.0f);
        }

//        volumeX.deleteSliceBlocks();
//        volumeY.deleteSliceBlocks();
        return true;
    }

    private int getVolumeSlice(StsSeismicVolume volume, float z)
    {
        if(z < volume.zMin || z > volume.zMax) return -1;
        return volume.getNearestSliceCoor(z);
    }
/*
    private byte[] checkLoadVolumeData(StsCursor3d cursor3d)
    {
        if(volumeDensityData != null) return volumeDensityData;

        if(volumes == null) return null;
        int nVolumes = volumes.getSize();
        if(nVolumes < 2) return null;

        StsSeismicVolume volumeX = (StsSeismicVolume)volumes.getElement(0);
        StsSeismicVolume volumeY = (StsSeismicVolume)volumes.getElement(1);

        // compute pointsPerBin (255 x 255 bins)
        int[][] pointsPerBin = new int[nRows][nCols];

        int max = 0;

        StsProject project = currentModel.getProject();
        StsRotatedGridBoundingBox cursorBoundingBox = project.getRotatedBoundingBox();
        StsCropVolume cropVolume = project.getCropVolume();

        StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
        if(!cropVolume.isCropped())
            boundingBox.classInitialize(cursorBoundingBox);
        else
            boundingBox.classInitialize(cropVolume);

        boolean subVolumeApplied = isSubVolumeApplied();
        if(mainDebug) System.out.println("subVolumeApplied: " + subVolumeApplied);
        if(subVolumeApplied)
        {
            StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)currentModel.getStsClass(StsSubVolume.class);
                StsGridBoundingBox subVolumeBoundingBox = subVolumeClass.getBoundingBox();
                if(subVolumeBoundingBox != null) boundingBox.intersectBoundingBox(subVolumeBoundingBox);
        }

        volumeX.setupSliceBlocks();
        volumeY.setupSliceBlocks();

        int nTotalCols = cursorBoundingBox.getNCols();
        int subVolRowStart = boundingBox.getRowMin();
        int subVolRowEnd = boundingBox.getRowMax();
        int subVolColStart = boundingBox.getColMin();
        int subVolColEnd = boundingBox.getColMax();

        int subVolIndex = 0;
        for (int subVolRow = subVolRowStart; subVolRow <= subVolRowEnd; subVolRow++)
        {
            int cursorIndex = subVolRow * nTotalCols + subVolColStart;
            for (int subVolCol = subVolColStart; subVolCol <= subVolColEnd; subVolCol++, cursorIndex++, n++)
            {
                if (subVolumePlane[subVolIndex] == 0) subVolumeData[subVolIndex] = -1;
                else subVolumeData[subVolIndex] = planeData[cursorIndex];
            }
        }
    }



        int nTotalRows = volumeX.getNRows();
        int nTotalCols = volumeY.getNCols();
//        int nTotalRows = projectBoundingBox.getNRows();
//        int nTotalCols = projectBoundingBox.getNCols();
        int nSlicePoints = nTotalRows*nTotalCols;

        MappedByteBuffer sliceXBuffer, sliceYBuffer;

        int subVolSliceMin = boundingBox.sliceMin;
        int subVolSliceMax = boundingBox.sliceMax;

        if(mainDebug) System.out.println("Computing volume crossplot density data for volume defined by bounding box: " + boundingBox.toString());

        byte[] subVolumePlane = null;
        for(int slice = sliceMin; slice <= sliceMax; slice++)
        {
            sliceXBuffer = volumeX.getSliceBuffer(slice);
            if(sliceXBuffer == null) continue;
            sliceYBuffer = volumeY.getSliceBuffer(slice);
            if(sliceYBuffer == null) continue;
            int slicePosition = volumeX.getSlicePosition(slice);

            boolean subVolumeOK = false;
            if(subVolumeApplied)
            {
                subVolumePlane = cursor3d.getSubVolumePlane(StsCursor3d.ZDIR, slice);
                subVolumeOK = (subVolumePlane != null);
            }
            for(int row = rowMin; row <= rowMax; row++)
            {
                sliceXBuffer.position(slicePosition + row*nTotalCols + colMin);
                sliceYBuffer.position(slicePosition + row*nTotalCols + colMin);
                int n = row*nTotalCols + colMin;
                for(int col = colMin; col <= colMax; col++, n++)
                {
                    if(subVolumeOK && subVolumePlane[n] == 0) continue;
                    int dataX = StsMath.signedByteToUnsignedInt(sliceXBuffer.get());
                    int dataY = StsMath.signedByteToUnsignedInt(sliceYBuffer.get());
                    if(dataX < 255 && dataY < 255)  // 255 value represents null
                    {
                        pointsPerBin[dataY][dataX]++;
                        max = Math.max(max, pointsPerBin[dataY][dataX]);
                    }
                }
            }
        }

        // byte 255 (signed byte -1) is the transparent color
        if(volumeDensityData == null) volumeDensityData = new byte[nRows*nCols];
        for(int n = 0; n < nRows*nCols; n++) volumeDensityData[n] = (byte)255;

        // scale to 0 to 254
        float scale = 254.0f/max;
        // hack so we have a reasonable amount of stuff to look at
        scale = 3.0f*scale;

        int n = 0;
        int idx = 0;
        for(int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++, n++)
            {
                int density = 0;
                if(pointsPerBin[i][j] > 0)
                {
                    density = StsMath.ceiling(pointsPerBin[i][j]*scale);
                    volumeDensityData[n] = unsignedIntToSignedByte254(density);
                }
//                densityData[n] = unsignedIntToSignedByte254((int)(pointsPerBin[i][j] * scale));
//                idx = StsMath.signedByteToUnsignedInt(densityData[n++]);
                if(density > 254) density = 254;
                dataHist[density]++;
            }
        }

        int numValidValues = (nRows*nCols) - (int)dataHist[0];
        dataHist[0] = 0;
        for(int i=0; i<255; i++)
            dataHist[i] = (float)((float)dataHist[i]/((float)numValidValues)*100.0f);

        if(max != colorscale.range()[1]) colorscale.setRange(0.0f, max);

        volumeX.deleteSliceBlocks();
        volumeY.deleteSliceBlocks();

        return volumeDensityData;
    }
*/
    public boolean isVolumeDisplayOK()
    {
        return getCrossplotClass().getDisplayEntireVolume() && volumeDensityData != null;
    }

    /** converts an unsigned int to a signedByte value between/including 0 to 254 */
    final public static byte unsignedIntToSignedByte254(int i)
    {
        if(i >= 255) i = 254;
        if(i < 0) i = 0;
        return (byte)i;
    }

/*
    public void checkSetCrop()
    {
        if(!cropIt) return null;
        // Crop Variables
        cropped = false;
        iCropMin = -100000.0f;
        iCropMax = 100000.0f;
        xCropMin = -100000.0f;
        xCropMax = 100000.0f;
        zCropMin = -100000.0f;
        zCropMax = 100000.0f;

        StsSeismicVolume volume = (StsSeismicVolume)volumes.getElement(0);
        setCropLimits(volume);
    }
*/
    private void setCropLimits(StsSeismicVolume sv)
    {
/*
        if(iCropMin < sv.getCropRowLabelMin())
            iCropMin = sv.getCropRowLabelMin();
        if(iCropMax > sv.getCropRowLabelMax())
            iCropMax = sv.getCropRowLabelMax();
        if(xCropMin < sv.getCropColLabelMin())
            xCropMin = sv.getCropColLabelMin();
        if(xCropMax > sv.getCropColLabelMax())
            xCropMax = sv.getCropColLabelMax();
        if(zCropMin < sv.getCropZMin())
            zCropMin = sv.getCropZMin();
        if(zCropMax > sv.getCropZMax())
            zCropMax = sv.getCropZMax();

        if(iCropMin > sv.getRowLabelMin() || iCropMax < sv.getRowLabelMax() || xCropMin > sv.getColLabelMin() ||
           xCropMax < sv.getColLabelMax() || zCropMin > sv.getZMin() || zCropMax < sv.getZMax())
            cropped = true;
*/
        return;
    }

    public boolean isCursorInsideCroppedLimits(StsSeismicVolume sv, int dirNo)
    {
//        return sv.isCoorInsideCropLimits(dirNo, planeCoor);
        return true;
    }

    public String getName() { return super.getName(); }
    public boolean getCropIt() { return cropIt; }
    public void setCropIt(boolean crop)
    {
        cropIt = crop;
//        reDraw = true;
        currentModel.win3dDisplayAll();
    }

    public float[][] getTotalAxisRanges() { return totalAxisRanges; }

    public StsToolbar getToolbar() { return toolbar; }
    public void setToolbar(StsToolbar tb)
    {
        toolbar = tb;
        return;
    }

    public StsColorscale getColorscale() { return colorscale; }
    public void setColorscale(StsColorscale colorscale) { this.colorscale = colorscale; }

    public boolean getIsPixelMode() { return getCrossplotClass().getIsPixelMode(); }
//    public boolean getCropMode() { return getCrossplotClass().getCropMode(); }

    public void addPolygon(StsXPolygon polygon)
    {
        polygons.add(polygon);
        polygonsChanged();
    }

    public void deletePolygon(StsXPolygon polygon)
    {
        polygons.delete(polygon);
        polygon.delete();
//        resetPolygonIDs();
        polygonsChanged();
    }

    public void deletePolygonPoint(StsXPolygon polygon, int index)
    {
        polygon.deletePoint(index);
        polygonsChanged();
    }

    public void polygonsChanged()
    {
        polygonsProcessed = false;
        polygonData = null;
    }

    public void drawPolygons(GL gl, GLU glu, boolean filled)
    {
        if(polygons == null) return;
        int nPolygons = polygons.getSize();
        for(int n = 0; n < nPolygons; n++)
        {
            StsXPolygon polygon = (StsXPolygon)polygons.getElement(n);
            polygon.draw(gl, glu, filled);
        }
    }

    public boolean getDisplayDensity()
    {
        StsToggleButton attributeOrDensityBtn =
            (StsToggleButton)currentModel.win3d.getToolbarComponentNamed(StsCrossplotToolbar.NAME, StsCrossplotToolbar.ATTRIBUTE_DENSITY);
        return !attributeOrDensityBtn.isSelected();
//        return displayDensity;
    }

    public void checkPolygons()
    {
        if(polygons == null || polygons.getSize() == 0) return;
        polygons.forEach("checkDeletePolygon", polygons);
    }

    public boolean processPolygons()
    {
        try
        {
            if(polygonsProcessed) return true;
            if(polygons == null) return false;
            int nPolygons = polygons.getSize();
            if(nPolygons == 0) return false;
			if(totalAxisRanges == null) return false;
            StsGLOffscreenPolygon offscreen = new StsGLOffscreenPolygon(255, 255,
                totalAxisRanges[0][0], totalAxisRanges[0][1], totalAxisRanges[1][0], totalAxisRanges[1][1], false);

            offscreen.setPolygons(polygons);
			offscreen.startGL(); // jbw
			offscreen.repaint(); // jbw
            offscreen.close();
            polygonData = offscreen.getData();
            polygonsProcessed = (polygonData != null);
            return polygonsProcessed;
        }
        catch(Exception e)
        {
            StsException.outputException("StsCrossplot.processPolygons() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
/*
    public void addCrossplotPoint(StsCursorPoint cursorPoint, StsColor color, byte colValue, byte rowValue)
    {
        StsCrossplotPoint crossplotPoint = new StsCrossplotPoint(cursorPoint, color, colValue, rowValue);
        crossplotPoints.add(crossplotPoint);
    }

    public void clearCrossplotPoints()
    {
        crossplotPoints = null;
    }
*/
    /** polygon ID is the sequence number in the ref list it belongs to;
     *  used to define color list for color-indexed texture map
     */
/*
    public void resetPolygonIDs()
    {
        int nPolygons = polygons.getSize();
        for(int n = 0; n < nPolygons; n++)
        {
            StsXPolygon polygon = (StsXPolygon)polygons.getElement(n);
            polygon.setID(n);
        }
    }
*/
/*
    public boolean processPolygons()
    {
        try
        {
            if(polygonsProcessed) return true;
            if(polygons == null || polygons.getSize() == 0) return false;
            if(offscreen == null)
            {
                offscreen = new StsGLOffscreenXPData(256, 256, 100.0f, this);
//                offscreen.startGL();
            }
//            else
            offscreen.repaint();

            polygonData = offscreen.getData();
            polygonsProcessed = (polygonData != null);
            return polygonsProcessed;
        }
        catch(Exception e)
        {
            StsException.outputException("StsCrossplot.processPolygons() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
*/
    /** Find all unique color-names used in this crossplot and
     *  classInitialize the colorListSelector with them.
     */
/*
    public void setCurrentCrossplotColors()
    {
        int nItems;

        StsColorListComboBox comboBox = getColorListComboBox();
        comboBox.removeAll();

        StsObjectRefList types = polygonTypeLibrary.getTypes();
        int nTypes = types.getSize();
        if(nTypes == 0) return;
        for(int n = 0; n < nTypes; n++)
        {
            StsType type = (StsType)types.getElement(n);
            String name = type.getName();
            Color color = type.getColor();
            comboBox.addItem(name, color);
        }
        comboBox.setSelectedIndex(1);
    }
*/
/*
    public StsColorListComboBox getColorListComboBox()
    {
        Component component = currentModel.win3d.getToolbarComponentNamed("Cross Plot Toolbar", "editColors");
        if(component == null) return null;
        StsEditColorListPanel panel = (StsEditColorListPanel)component;
        panel.toggleDisplay(true);
        StsColorListComboBox comboBox = panel.getComboBox();
        return comboBox;
    }
*/
    public StsType getPolygonType(String name, Color color)
    {
        return polygonTypeLibrary.getType(name, new StsColor(color));
    }

    public StsType getPolygonType(String name, StsColor color)
    {
        return polygonTypeLibrary.getType(name, color);
    }

    public StsTypeLibrary getTypeLibrary() { return polygonTypeLibrary; }
    public void setTypeLibrary(StsTypeLibrary lib)
    {
        polygonTypeLibrary = lib;
    }

//    public byte[] getPolygonData() { return polygonData; }

    // displays the crossplot on a seismic section in 3d
/*
    public void displayOnCursor(int dirNo, float dirCoordinate, GL gl, GLU glu)
      {
          StsSeismicVolume volume;
          StsSeismicCursorSection[] seismicSections;

  //        if(!processPolygons()) return;
          if(polygonData == null) return;

          if(volumes.getSize() < 2) return;

          seismicSections = new StsSeismicCursorSection[2];
          for(int n = 0; n < 2; n++)
          {
              volume = (StsSeismicVolume)volumes.getElement(n);
              seismicSections[n] = volume.getCursorSection(dirNo);
              if(seismicSections[n] == null)
              {
                  StsMessageFiles.infoMessage("Must select a seismic slice first.");
                  return;
              }
          }

          int nPlane = seismicSections[0].getCursorPlaneIndex(dirCoordinate);
          if(nPlane < 0) return;

          int nRows = seismicSections[0].nRows;
          int nCols = seismicSections[0].nCols;

          StsCursor3dTexture cursorSection = cursorSectionDisplayables[dirNo];
          if(cursorSection == null)
          {
              cursorSection = new StsCursor3dTexture(seismicSections[0]);
              cursorSectionDisplayables[dirNo] = cursorSection;
          }

          StsGLDraw.enableTransparentOverlay(gl);

          gl.glDisable(GL.GL_LIGHTING);
          gl.glEnable(GL.GL_TEXTURE_2D);
          gl.glEnable(GL.GL_BLEND);
          gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
          gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
          gl.glShadeModel(GL.GL_FLAT);

          createPolygonColorList(gl);

          cursorSection.bindTexture(gl);

          if(checkAddData(cursorSection, nPlane, seismicSections, gl))
              cursorSection.display(gl, seismicSections[0].planePoints);

          gl.glDisable(GL.GL_TEXTURE_2D);
          gl.glDisable(GL.GL_BLEND);
          gl.glEnable(GL.GL_LIGHTING);

          StsGLDraw.disableTransparentOverlay(gl);
      }
  */
    // displays the crossplot on a seismic section in 2d
    public void displayOnCursor2d(int dirNo, float dirCoordinate, GL gl, GLU glu)
	{
/*
        StsSeismicVolume volume;
        StsSeismicCursorSection[] seismicSections;
        byte[][] seismicData;

//        if(!processPolygons()) return;
        if(polygonData == null) return;

        if(volumes.getSize() < 2) return;

        seismicSections = new StsSeismicCursorSection[2];
        for(int n = 0; n < 2; n++)
        {
            volume = (StsSeismicVolume)volumes.getElement(n);
            seismicSections[n] = volume.getCursorSection(dirNo);
            if(seismicSections[n] == null)
            {
                StsMessageFiles.infoMessage("Must select a seismic slice first, or not initialized yet.");
                return;
            }
        }

        int nPlane = seismicSections[0].getCursorPlaneIndex(dirCoordinate);
        if(nPlane < 0) return;

        StsCursor3dTexture cursorSection = cursorSectionDisplayables[dirNo];
        if(cursorSection == null)
        {
            cursorSection = new StsCursor3dTexture(seismicSections[0]);
            cursorSectionDisplayables[dirNo] = cursorSection;
        }

        StsGLDraw.enableTransparentOverlay(gl);

        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);

        createPolygonColorList(gl);

        cursorSection.bindTexture(gl);

        if(checkAddData(cursorSection, nPlane, seismicSections, gl))
            cursorSection.display2d(gl, seismicSections[0].planePoints);

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

        StsGLDraw.disableTransparentOverlay(gl);
*/
    }

    public byte[] readBytePlaneData(int dir, float dirCoordinate)
    {
        if (polygonData == null)return null;

        StsSeismicVolume[] volumes = getVolumes();

        // check that all volumes are available; if not, delete this volume
        for (int n = 0; n < 2; n++)
        {
            StsObject volume = (StsObject) volumes[n];
            StsClass volumeClass = currentModel.getStsClass(volume.getClass());
            if (!volumeClass.contains(volume))
            {
                getStsClass().delete(this);
                return null;
            }
        }

        byte[][] seismicData = new byte[2][];
        for (int n = 0; n < 2; n++)
        {
            seismicData[n] = volumes[n].readBytePlaneData(dir, dirCoordinate);
            if (seismicData[n] == null)
            {
                StsMessageFiles.infoMessage("StsCrossplot.computeCursorDisplayData() failed." +
                                            " seismic data unavailable from volume: " + volumes[n].getName());
                return null;
            }
        }
        if (seismicData[0].length != seismicData[1].length)
        {
            StsException.systemError("StsCrossplot.displayOnCursor() failed." +
                                     " seismic data from cursor not the same size for 2 volumes.");
            return null;
        }

        int nRows = volumes[0].getNCursorRows(dir);
        int nCols = volumes[0].getNCursorCols(dir);

        byte[] displayData = new byte[nRows * nCols];
        int n = 0;
        int nDataPoints = 0;
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                int v0 = StsMath.signedByteToUnsignedInt(seismicData[0][n]);
                int v1 = StsMath.signedByteToUnsignedInt(seismicData[1][n]);
                if (v0 < 255 && v1 < 255) // we assume unsigned 255 is a null
                {
                    byte d = polygonData[v1 * 255 + v0];
                    displayData[n] = d;
                    if (d > 0)
                        nDataPoints++;
                }
                n++;
            }
        }
        if(debug) System.out.println("Number of data points on seismic section " + nDataPoints);

        return displayData;
    }

    public final float getScaledValue(byte byteValue)
    {
        return (float)StsMath.signedByteToUnsignedInt(byteValue);
    }

    public final boolean isByteValueNull(byte byteValue)
    {
         return byteValue == 0;
    }

    public int getCursorPlaneIndex(int dirNo, float dirCoordinate)
    {
        StsSeismicVolume volume = (StsSeismicVolume)volumes.getElement(0);
        return volume.getCursorPlaneIndex(dirNo, dirCoordinate);
    }

/*
    public byte[] compute3dDisplayData()
    {
        if(polygonData == null) return null;

        byte[][] seismicData = new byte[2][];
        StsSeismicVolume[] volumes = volumes();
        for(int n = 0; n < 2; n++)
        {
            seismicData[n] = volumes[n].readIndexPlaneData(dirNo, nPlane);
            if(seismicData[n] == null)
            {
                StsMessageFiles.infoMessage("StsCrossplot.computeCursorDisplayData() failed." +
                    " seismic data unavailable from volume: " + volumes[n].getName());
                return null;
            }
        }
        if(seismicData[0].length != seismicData[1].length)
        {
            StsException.systemError("StsCrossplot.displayOnCursor() failed." +
                " seismic data from cursor not the same size for 2 volumes.");
            return null;
        }

        int nRows = volumes[0].getCursorNRows(dirNo);
        int nCols = volumes[0].getCursorNCols(dirNo);

        int nPoints = nRows*nCols;
        byte[] isDisplayData = new byte[nRows*nCols];
        int n = 0;
        int nDataPoints = 0;
        for(int row = 0; row < nRows; row++)
        {
            for(int col = 0; col < nCols; col++)
            {
                int v0 = StsMath.signedByteToUnsignedInt(seismicData[0][n]);
                int v1 = StsMath.signedByteToUnsignedInt(seismicData[1][n]);
                byte d = polygonData[v1*256 + v0];
                isDisplayData[n++] = d;
                if(d > 0) nDataPoints++;
            }
        }
//        System.out.println("Number of data points on seismic section " + cursorSection.dirNo + ": " + nDataPoints);

        return isDisplayData;
    }
*/
    private void createPolygonColorList(GL gl)
    {
        StsColor[] colors = polygonTypeLibrary.getStsColors();
        int nColors = colors.length + 1;
        int nAllocatedColors = StsMath.nextBaseTwoInt(nColors);
        float[][] arrayRGBA = new float[4][nAllocatedColors];

        float[] rgba = new float[4];
        for(int n = 1; n < nColors; n++)
        {
            rgba = colors[n-1].getRGBA();
            for(int c = 0; c < 4; c++) arrayRGBA[c][n] = rgba[c];
        }
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nAllocatedColors, arrayRGBA[0], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nAllocatedColors, arrayRGBA[1], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nAllocatedColors, arrayRGBA[2], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nAllocatedColors, arrayRGBA[3], 0);
        gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
    }

    public StsFieldBean[] getDisplayFields() { return displayFields; }
    public StsFieldBean[] getPropertyFields() { return propertyFields; }
    public Object[] getChildren() { return new Object[0]; }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.setCurrentObject(this);
        currentModel.getGlPanel3d().checkAddView(StsViewXP.class);
//        currentModel.glPanel3d.cursor3d.setCurrentCrossplot(this);
        currentModel.win3dDisplayAll();
    }
/*
    public void resetAxisRanges()
    {
        axisRanges = new float[2][2];
        axisLabels = new String[2];

        for(int n = 0; n < 2; n++)
        {
            StsSeismicVolume volume = (StsSeismicVolume)this.volumes.getElement(n);
            axisRanges[n][0] = volume.dataMin;
            axisRanges[n][1] = volume.dataMax;
            axisLabels[n] = volume.getName();
        }
        totalAxisRanges = StsMath.copyFloatArray(axisRanges);
    }
*/
    public boolean anyDependencies()
    {
        return false;
    }
    public boolean canExport() { return true; }
    public boolean export()
    {
        PrintWriter printWriter;
        double[] point = new double[2];

        if (polygons == null || polygons.getSize() == 0)
        {
            StsMessageFiles.infoMessage("No polygons defined for this crossplot.");
            return false;
        }

        StsModel model = StsObject.getCurrentModel();
        String name =  StsStringUtils.cleanString(getName());
        String filename = model.getProject().getRootDirString() + "polygons.txt." + name;

        try
        {
            printWriter = new PrintWriter(new FileWriter(filename, false));
            printWriter.println("<Header>");
            printWriter.println("Cross Plotting Polygons");

            StsSeismicVolume[] volumes = getVolumes();
            printWriter.println("<Data>");
            printWriter.println("Seismic Volumes   " + volumes.length);
            for(int j=0; j<volumes.length; j++)
                printWriter.println("       " + volumes[j].getName());

            NumberFormat number = NumberFormat.getInstance();
            number.setGroupingUsed(false);
            number.setMaximumFractionDigits(4);
            number.setMinimumFractionDigits(2);

            for(int i=0; i<polygons.getSize(); i++)
            {
                StsXPolygon polygon = (StsXPolygon) polygons.getElement(i);
                printWriter.println("<Polygon>");
                printWriter.println(polygon.libraryType.getName() + "   " + polygon.libraryType.stsColor.red + "  " + polygon.libraryType.stsColor.green
                                    + "  " + polygon.libraryType.stsColor.blue + "  " + polygon.getNPoints());

                for(int j=0; j < polygon.getNPoints(); j++)
                {
                    point = polygon.getPoint(j);
                    printWriter.println("       " + number.format(point[0]) + "        " + number.format(point[1]));
                }
            }
            printWriter.close();

            StsMessageFiles.infoMessage("Polygon exported to file: " + model.getProject().getRootDirString()
                                    + "polygons.txt." + name);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrum.export() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
//
    public int getColorListNum(GL gl, StsColorscale colorscale)
    {
        if(colorListChanged)
        {
            gl.glDeleteLists(colorListNum, 1);
            colorListNum = 0;
//            deleteTexture(gl);
            colorListChanged = false;
        }

        if(colorListNum != 0) return colorListNum;

        colorListNum = gl.glGenLists(1);
        if(colorListNum == 0)
        {
            StsMessageFiles.logMessage("System Error in StsCrossplot.createColorListNum(): Failed to allocate a display list");
            return 0;
        }

        gl.glNewList(colorListNum, GL.GL_COMPILE);
        createColorList(gl, colorscale.getNewColorsInclTransparency());
        gl.glEndList();

        return colorListNum;
    }
    private void createColorList(GL gl, Color[] colors)
    {
        if(arrayRGBA == null) setRGBAArray(colors);
        createColorList(gl);
        arrayRGBA = null;
    }

    private void createColorList(GL gl)
    {
        if(arrayRGBA == null) return;
        int nColors = arrayRGBA[0].length;
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
        gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
    }

    /**
     * Reset the colors for dynamic update of the crossplot colors
     * @param colors the new colors
     */
    private void setRGBAArray(Color[] colors)
    {
        int nColors = colors.length;

        float norm = 1.0f/255.0f;

        if(arrayRGBA == null || arrayRGBA[0].length != nColors)
            arrayRGBA = new float[4][nColors];

        for(int n = 0; n < nColors; n++)
        {
            arrayRGBA[0][n] = colors[n].getRed()*norm;
            arrayRGBA[1][n] = colors[n].getGreen()*norm;
            arrayRGBA[2][n] = colors[n].getBlue()*norm;
            float alpha = colors[n].getAlpha()*norm;
            if(alpha < 0.0f) alpha = 1.0f;
            arrayRGBA[3][n] = alpha;
        }

        if(backgroundTransparent)
            for(int i = 0; i < 4; i++) arrayRGBA[i][0] = 0.0f;
    }

    /** If the crossplot has changed, then clear the colors and display and force redraw */
    public void clearColors()
    {
        colorListChanged = true;
    }

    public StsColorListItem getTypeItem() {  return polygonTypeLibrary.getTypeItem(); }
    public void setTypeItem(StsColorListItem colorListItem) { polygonTypeLibrary.setTypeItem(colorListItem); }

    public Class[] getViewClasses() { return new Class[] { StsViewXP.class }; }

/*
    public void clearTextureDisplay()
    {
        clearDataDisplay();
    }

    public void clearDataDisplay()
    {
        clearData();
    }
*/
}

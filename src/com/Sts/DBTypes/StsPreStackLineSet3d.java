package com.Sts.DBTypes;

import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsPreStackLineSet3d extends StsPreStackLineSet implements ItemListener, StsTreeObjectI, StsViewable
{
    transient public StsPreStackVolume[] outputVolumes;
    
    final static public String volumeFilePrefix = "seis.vol.stack.";

    static StsComboBoxFieldBean previousVelocityModelBean;
 //	static StsComboBoxFieldBean currentVelocityModelBean = new StsComboBoxFieldBean(StsPreStackLineSet3d.class, "velocityModel", "Current Velocity:", null);
	static public StsFieldBean[] displayFields3d = null;
    static public StsFieldBean[] propertyFields3d = null;

    static public final String group = "prestackVol"; // not used currently; prestack data saved in velocity, stack, and velocity instances of StsPreStackVolume

    static final byte DISPLAY_VELOCITY = StsPreStackLineSetClass.DISPLAY_VELOCITY;
    static final byte DISPLAY_STACKED = StsPreStackLineSetClass.DISPLAY_STACKED;
    static final byte DISPLAY_SEMBLANCE = StsPreStackLineSetClass.DISPLAY_SEMBLANCE;
    static final byte DISPLAY_NONE = StsPreStackLineSetClass.DISPLAY_NONE;
    static final byte DISPLAY_ATTRIBUTE = StsPreStackLineSetClass.DISPLAY_ATTRIBUTE;

    static final int nVolumeTypes = StsPreStackLineSetClass.nVolumeTypes;

    public String getGroupname()
    {
        return group3dPrestack;
    }

    public StsPreStackLineSet3d()
	{
	}

	public StsPreStackLineSet3d(boolean persistent)
	{
		super(persistent);
	}

	private StsPreStackLineSet3d(String name, StsFile[] files, StsModel model, StsProgressPanel panel) throws FileNotFoundException
	{
		super(false);

		setName(name);
//		toggleZDomain(StsProject.TD_TIME);
		int nFiles = files.length;
        panel.initialize(nFiles);
        stsDirectory = files[0].getDirectory();
        for(int n = 0; n < nFiles; n++)
		{
			if(!files[n].exists())
			{
				panel.setDescriptionAndLevel("File not found for file " + files[n].getPathname(), StsProgressBar.ERROR);
				new StsMessage(model.win3d, StsMessage.WARNING, "Failed to find file " + files[n].getPathname() + ".");
				throw new FileNotFoundException();
			}
		}
		ArrayList preStackLines = new ArrayList(100);
		boolean originSet = false;
		for(int n = 0; n < nFiles; n++)
        {
//			checkMemoryStatus("Before line "+n);
			StsPreStackLine3d seismicLine3d = StsPreStackLine3d.constructor(model, files[n], this);
			if(seismicLine3d == null)continue;
			checkSetOrigin(seismicLine3d);

			// TODO:
			// dataMin & dataMax could be different for volume if lines were processed separately.
			// We really can't allow this as it makes the texture mapping enormously more expensive
			// to translate each unsigned byte on a byte-by-byte basis to the correct index into the
			// colorscale.  In this case, we are better off reading in floats, scaling them to bytes and
			// putting them in the texture array.
			dataMin = Math.min(dataMin, seismicLine3d.dataMin);
			dataMax = Math.max(dataMax, seismicLine3d.dataMax);
			preStackLines.add(seismicLine3d);

			// haque until we figure this out....
			if(seismicLine3d.yInc < 0.0f) seismicLine3d.yInc = -seismicLine3d.yInc;

			addLineRotatedBoundingBox(seismicLine3d);
			setTraceOffsetRange(seismicLine3d);
//            addAttributes(seismicLine3d);
            panel.setValue(n+1);
            panel.setDescription("Loaded line : " + (n+1) + " of " + nFiles);
            panel.appendLine("Completed loading line index:" + seismicLine3d.getIndex() + " from file " + seismicLine3d.getName());
		}
        panel.finished();
        stsDirectory = files[0].getDirectory();
		isVisible = true;
		//       setVelStat(lineSet.getShowVelStat());

		setRowColNumbering(preStackLines);

		sortLines(preStackLines);
    /*
        if(!allocateVolumes("rw"))
        {
            throw new StsException(StsException.WARNING, "Failed to allocate stack/semblance/velocity volumes.");
        }
    */
//		sortLines2d();
//		initializeStsClass();
//		setCurrentLineAndGather2d(currentModel.glPanel3d, 0, 0);
//		allocateVolumes();
		addToolbar();
		currentModel.getProject().setCursorDisplayXYAndGridCheckbox(false);

//		stackedVolume = StsPreStackVolume.constructor(this, StsPreStackSeismicClass.DISPLAY_MODE_STACKED);
//		semblanceVolume = StsPreStackVolume.constructor(this, StsPreStackSeismicClass.DISPLAY_MODE_SEMBLANCE);
	}

	static public StsPreStackLineSet3d constructor(String name, StsFile[] files, StsModel model, StsProgressPanel panel)
	{
		try
		{
			StsPreStackLineSet3d volume = new StsPreStackLineSet3d(name, files, model, panel);
			return volume;
		}
		catch(FileNotFoundException fnfe)
		{
			return null;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackLineSet3d.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}

    private void setRowColNumbering(ArrayList preStackLines)
	{
		int nLines = preStackLines.size();
		if(nLines == 0)return;
		StsPreStackLine3d line3d = (StsPreStackLine3d)preStackLines.get(0);
		rowNumMin = line3d.rowNumMin;
		rowNumMax = line3d.rowNumMax;
		colNumMin = line3d.colNumMin;
		colNumMax = line3d.colNumMax;
		rowNumInc = line3d.rowNumInc;
		colNumInc = line3d.colNumInc;
        isXLineCCW = line3d.isXLineCCW;

        for(int n = 1; n < nLines; n++)
		{
			line3d = (StsPreStackLine3d)preStackLines.get(n);
			rowNumMin = Math.min(line3d.rowNumMin, rowNumMin);
			rowNumMax = Math.max(line3d.rowNumMax, rowNumMax);
			colNumMin = Math.min(line3d.colNumMin, colNumMin);
			colNumMax = Math.max(line3d.colNumMax, colNumMax);
		}
		for(int n = 0; n < nLines; n++)
		{
			line3d = (StsPreStackLine3d)preStackLines.get(n);
			line3d.initializeLineIndex(this);
			line3d.initializeGatherIndexes();
		}
        if(!isXLineCCW)
        {
            float temp = rowNumMin;
            rowNumMin = rowNumMax;
            rowNumMax = temp;
            rowNumInc = -rowNumInc;
        }
    }

    public void applyEditedBox(StsEditedBoundingBox editedBox)
    {
        initializeToBoundingBox(editedBox);
    }

    public void completeLoading()
    {
		if (!addToProject())
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to add volume " + name + " to project.");
		}
        initialize();
        setCurrentDataRowCol(currentModel.win3d, 0, 0);
    }

    public StsPreStackVelocityModel constructVelocityModel()
    {
        return new StsPreStackVelocityModel3d(this);
    }

    public int getNColsForRow(int row)
    {
        return nCols;
    }
/*
	private void setRowColNumbering2d()
	{
		int nLines = preStackLines2d.size();
		if(nLines == 0)return;
		StsPreStackLine2d line = (StsPreStackLine2d)preStackLines2d.get(0);

		rowNumMin = line.rowNumMin;
		rowNumMax = line.rowNumMax;
		colNumMin = line.colNumMin;
		colNumMax = line.colNumMax;
		rowNumInc = line.rowNumInc;
		colNumInc = line.colNumInc;
		for(int n = 1; n < nLines; n++)
		{
			line = (StsPreStackLine2d)preStackLines2d.get(n);
			rowNumMin = Math.min(line.rowNumMin, rowNumMin);
			rowNumMax = Math.max(line.rowNumMax, rowNumMax);
			colNumMin = Math.min(line.colNumMin, colNumMin);
			colNumMax = Math.max(line.colNumMax, colNumMax);
		}
		for(int n = 0; n < nLines; n++)
		{
			line = (StsPreStackLine2d)preStackLines2d.get(n);
			line.initializeLineIndex(this);
			line.initializeGatherIndexes();
		}
		rowNumMax = Math.max(line.rowNumMax, 10);
		rowNumInc = Math.max(line.rowNumInc, 1);
		colNumMax = Math.max(line.colNumMax, 10);
		colNumInc = Math.max(line.colNumInc, 1);
	}
*/
	public boolean addToProject()
	{
		if(!currentModel.getProject().addToProject(this, true))return false;
		return true;
	}
/*
	public boolean initialize(StsFile file, StsModel model)
	{
		try
		{
			String pathname = file.getDirectory() + file.getFilename();
			StsParameterFile.readObjectFields(pathname, this, StsSeismicBoundingBox.class, StsBoundingBox.class);
			setName(getStemname());
			stsDirectory = file.getDirectory();
			if(!initialize(model))
			{
				return false;
			}
			isVisible = true;
			((StsPreStackLineSet3dClass)lineSetClass).setIsVisibleOnCursor(true);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackLineSet3d.loadFile() failed.", e, StsException.WARNING);
			return false;
		}
	}
*/
	// if crossline+ direction is 90 degrees CCW from inline+, this is isColCCW; otherwise not
	// angle is from X+ direction to inline+ direction (0 to 360 degrees)
	public boolean initialize(StsModel model)
	{
		if (!super.initialize(model)) return false;
//		initializeVolumeLists(); // now done when displayFields are fetched
		return true;
	}

	private void sortLines(ArrayList preStackLines)
	{
		int nLines = preStackLines.size();

		Collections.sort
        (
            preStackLines,
            new Comparator()
            {
                public int compare(Object a, Object b)
                {
                    StsPreStackLine3d line3dA = (StsPreStackLine3d)a;
                    StsPreStackLine3d line3dB = (StsPreStackLine3d)b;
                    return line3dA.compareTo(line3dB);
                }
            }
		);
        lines = new StsPreStackLine3d[nLines];
        for(int n = 0; n < nLines; n++)
		{
			StsPreStackLine3d line3d = (StsPreStackLine3d)preStackLines.get(n);
//			initializeLineIndex(line3d);
			lines[n] = line3d;
		}
	}
/*
    public void initializeLineIndex(StsPreStackLine line)
    {
        if (isInline)
        {
            line.lineNum = rowNumMin;
            line.lineIndex = getNearestRowCoor(yMin);
        }
        else
        {
            line.lineNum = colNumMin;
            line.lineIndex = getNearestColCoor(xMin);
        }
    }
*/
	public StsPreStackLine getDataLine(int row, int col)
	{
		if (isInline)
		{
			return getLineFromIndex(row, col, nRows);
		}
		else
		{
			return getLineFromIndex(col, row, nCols);
		}
	}

	private StsPreStackLine getLineFromIndex(int lineIndex, int gatherIndex, int nLinesMax)
	{
		if (lines == null)
		{
			return null;
		}
		int nLines = lines.length;
		if (lineIndex < 0 || lineIndex >= nLinesMax)
		{
			return null;
		}
		int nLine;
		// if we don't have a line for every row, then search from beginning;
		// otherwise search from this lineIndex
		if (lineIndex >= nLines)
		{
			nLine = 0;
		}
		else
		{
			nLine = lineIndex;
		}

		StsPreStackLine line = lines[nLine];
		int compare = line.compareTo(lineIndex, gatherIndex);
		if (compare == 0)
		{
			return line;
		}
		while (true)
		{
			nLine -= compare;
			if (nLine < 0 || nLine >= nLines)
			{
				return null;
			}
			line = lines[nLine];
			int newCompare = line.compareTo(lineIndex, gatherIndex);
			if (newCompare == 0)
			{
				return line;
			}
			if (newCompare != compare)
			{
				return null; // compare has passed it, so bail
			}
		}
	}
/*
	private void setCurrentVelocityProfile(int row, int col)
	{
		if(velocityModel == null)return;
		lineSetClass.setDisplayStack(false);
		velocityModel.setCurrentVelocityProfile(row, col);
	}
*/
/*
    public StsVelocityProfile getComputeVelocityProfile(int row, int col)
	{
		if(velocityModel == null)return null;
		return velocityModel.getExistingVelocityProfile(row, col);
	}
*/   
    public void setVelocityModel(StsPreStackVelocityModel newVelocityModel)
	{
        if (newVelocityModel == null) return;
        super.setVelocityModel(newVelocityModel);
        allocateVolumes("rw", false);
    }

    public void setInputVelocityModel(Object object)
	{
		if (object instanceof StsPreStackVelocityModel3d)
		{
			StsPreStackVelocityModel newInputVelocityModel = (StsPreStackVelocityModel3d)object;
			dbFieldChanged("inputVelocityModel", inputVelocityModel, newInputVelocityModel);
		}
		else if (object instanceof StsSeismicVolume)
		{
			// Create a Velocity Model
			StsPreStackVelocityModel newInputVelocityModel = new StsPreStackVelocityModel3d(this, this.getName() + "-Input");
			newInputVelocityModel.addToModel();
			dbFieldChanged("inputVelocityModel", inputVelocityModel, newInputVelocityModel);
			inputVelocityModel.setPreStackLineSet((StsPreStackLineSet)object);
		}
		else
		{
			dbFieldChanged("inputVelocityModel", inputVelocityModel, null);
		}
	}
	/** Allocate a "volume" which is a plane for each gather and there are nRows*nCols gathers.
	 *  A gather plane has a maximum size of nOffsetsMax by nCroppedSlices. Make this a type "2" (ZDIR)
	 *  volume so that nRows corresponds to nOffsetsMax, nCols corresponds to nCroppedSlices, and nPlanes
	 *  corresponds to nRows*nCols.
	 */
	public boolean allocateVolumes(String mode, boolean loadFromFiles)
	{
        if(velocityModel == null) return true;
        if(outputVolumes != null) return true;
        outputVolumes = new StsPreStackVolume[nVolumeTypes];
        for(byte n = 0; n < nVolumeTypes; n++)
        {
            outputVolumes[n] = StsPreStackVolume.loadConstruct(currentModel, this, n, stemname, loadFromFiles);
            if(outputVolumes[n] == null) return false;
        }
		return true;
	}
/*
    public void setOutputVolumeTime(int type, long time)
    {
        outputVolumeTimes[type] = time;
        dbFieldChanged("outputVolumeTimes", outputVolumeTimes);
    }
*/
    public StsPreStackVolume getVelocityVolume()
    {
        return outputVolumes[DISPLAY_VELOCITY];
    }

    public StsPreStackVolume getStackedVolume()
    {
        return outputVolumes[DISPLAY_STACKED];
    }

    public StsPreStackVolume getSemblanceVolume()
    {
        return outputVolumes[DISPLAY_SEMBLANCE];
    }

    public String getGroup()
    {
        return group;
    }

    public void jumpToRowCol(int[] rowCol, StsGLPanel3d glPanel3d)
	{
		if (rowCol == null) return;
		int row = rowCol[0];
		int col = rowCol[1];

        if(velocityModel != null) velocityModel.checkCurrentVelocityProfile(glPanel3d.window);

        StsSuperGather gather = getSuperGather(glPanel3d.window);
		if(gather.superGatherRow == row && gather.superGatherCol == col) return;
//		int dirNo = currentModel.getCursor3d().getCurrentDirNo();
		float currentX = xMin;
		float currentY = yMin;
        StsWindowFamily windowFamily = currentModel.getWindowFamily(glPanel3d.window);
        if (gather.superGatherRow != row)
		{
			currentY = getYCoor(row);
            windowFamily.adjustCursorAndSlider(YDIR, currentY);
//            currentModel.adjustCursorAndSlider(YDIR, currentY, glPanel3d.window);
		}
		if (gather.superGatherCol != col)
		{
			currentX = getXCoor(col);
            windowFamily.adjustCursorAndSlider(XDIR, currentX);
//            currentModel.adjustCursorAndSlider(XDIR, currentX, glPanel3d.window);
		}
    /*
        if (isInline)
		{
			currentModel.getCursor3d(glPanel3d.window).setCurrentDirNo(YDIR);
		}
		else
		{
			currentModel.getCursor3d(glPanel3d.window).setCurrentDirNo(XDIR);
		}
    */
		if (debug)
		{
			System.out.println(" row " + row + " col " + col + " x " + currentX + " y " + currentY);
		}
		super.jumpToRowCol(rowCol, glPanel3d.window);
	 }

     public int[] adjustLimitRowCol(int nRow, int nCol)
     {
         nRow = StsMath.minMax(nRow, 0, nRows-1);
         nCol = StsMath.minMax(nCol, 0, nCols-1);
         return new int[] { nRow, nCol };
     }

    public int[] getRowColFromCoors(float x, float y)
    {
        int row = this.getNearestBoundedRowCoor(y);
        int col = getNearestBoundedColCoor(x);
        return new int[] { row, col } ;
    }
/*
	public int[] advanceRowCol(int row , int col , int direction)
	{
		int nGridRows = (int) ( (nRows - 1.0f) / analysisRowInc) + 1;
		int nGridCols = (int) ( (nCols - 1.0f) / analysisColInc) + 1;
		int gridRow = row / analysisRowInc;
		int gridCol = col / analysisColInc;
		if (mainDebug)
		{
			System.out.println("maxRows=" + nGridRows + " maxCols=" + nGridCols + " currentRow=" + gridRow + " currentCol=" + gridCol);
		}

		int maxGridIndex = nGridRows * nGridCols - 1;
		if (isInline)
		{
			int gridIndex = gridRow * nGridCols + gridCol;
			if (direction > 0 || col % analysisColInc == 0)
			{
				gridIndex += direction;
			}
			if (gridIndex > maxGridIndex)
			{
				gridIndex = 0;
			}
			else if (gridIndex < 0)
			{
				gridIndex = maxGridIndex;
			}
			gridRow = gridIndex / nGridCols;
			gridCol = gridIndex - gridRow * nGridCols;
		}
		else
		{
			int gridIndex = gridCol * nGridRows + gridRow;
			if (direction > 0 || row % analysisRowInc == 0)
			{
				gridIndex += direction;
			}
			if (gridIndex > maxGridIndex)
			{
				gridIndex = 0;
			}
			else if (gridIndex < 0)
			{
				gridIndex = maxGridIndex;
			}
			gridCol = gridIndex / nGridRows;
			gridRow = gridIndex - gridCol * nGridRows;
		}
		row = gridRow * analysisRowInc;
		col = gridCol * analysisColInc;
		return new int[]
			{row , col};
	}
*/
    /*
    private void initializePlaneStatus()
	{
		StsDBMethodCmd cmd = new StsDBMethodCmd(stackedVolume, "initialisePlaneOKFlags");
		currentModel.addTransactionCmd("processVolume initialisePlaneOKFlags", cmd);
		stackedVolume.initialisePlaneOKFlags();
		cmd = new StsDBMethodCmd(semblanceVolume, "initialisePlaneOKFlags");
		currentModel.addTransactionCmd("processVolume initialisePlaneOKFlags", cmd);
		semblanceVolume.initialisePlaneOKFlags();
	}
    */
    public void resetColors()
	{
		for(int i = 0; i < currentModel.viewPersistManager.families.length; i++)
		{
			StsWin3dBase[] windows = currentModel.getWindows(i);
			for(int n = 0; n < windows.length; n++)
			{
				StsWin3dBase window = windows[n];
				window.getCursor3d().objectChanged(this);
			}
		}
	}

	public boolean getIsPixelMode()
	{
		return lineSetClass.getIsPixelMode();
	}

    public StsFieldBean[] getDisplayFields()
	{
        if(displayFields3d == null)
        {
            // colorscalesBean = new StsComboBoxFieldBean(StsPreStackLineSet3d.class, "currentColorscale", "Colorscales");
            colorscaleBean = new StsEditableColorscaleFieldBean(StsPreStackLineSet3d.class, "currentColorscale");
            displayCultureBean = new StsComboBoxFieldBean(StsPreStackLineSet3d.class, "displayCulture", "S/R Attribute");

            StsIntFieldBean traceThresholdBean = new StsIntFieldBean(StsPreStackLineSet3d.class, "traceThreshold", 0, 100, "Minimum Traces Required:", true);
            traceThresholdBean.setToolTipText("Minimum number of traces required in a gather for velocity analysis");

            displayFields3d =  new StsFieldBean[]
            {
                new StsBooleanFieldBean(StsPreStackLineSet3d.class, "isVisible", "Enable"),
                new StsBooleanFieldBean(StsPreStackLineSet3d.class, "readoutEnabled", "Mouse Readout"),
                new StsButtonFieldBean("Wiggle Display Properties", "Edit Wiggle Display Properties.", this, "displayWiggleProperties"),
                new StsButtonFieldBean("Semblance Display Properties", "Edit Semblance Display Properties.", this, "displaySemblanceProperties"),
                new StsButtonFieldBean("CVS/VVS Display Properties", "Edit CVS Display Properties.", this, "displayCVSProperties"),
                new StsComboBoxFieldBean(StsPreStackLineSet3d.class, "displayAttribute", "Display Attribute", "displayAttributes"),
                new StsBooleanFieldBean(StsPreStackLineSet3d.class, "showSource", "Show Sources"),
                new StsBooleanFieldBean(StsPreStackLineSet3d.class, "showReceiver", "Show Receivers"),
                new StsIntFieldBean(StsPreStackLineSet3d.class, "analysisRowInc", 0, 100, "Row Analysis Increment:", true),
                new StsIntFieldBean(StsPreStackLineSet3d.class, "analysisColInc", 0, 100, "Col Analysis Increment:", true),
                traceThresholdBean,
                displayCultureBean,
                colorscalesBean = new StsComboBoxFieldBean(StsPreStackLineSet.class, "currentColorscale", "Colorscales", "colorscaleList"),
                colorscaleBean
            };
        }


        // colorscalesBean.setListItems(colorscales.getElements());
        // colorscalesBean.setSelectedItem(currentColorscale);
        colorscaleBean.setValueObject(currentColorscale);
        colorscaleBean.setHistogram(dataHist);
        displayCultureBean.setListItems(displayCultures);
        return displayFields3d;
	}

	public StsFieldBean[] getPropertyFields()
	{
        if(propertyFields3d == null)
        {
            previousVelocityModelBean = new StsComboBoxFieldBean(StsPreStackLineSet3d.class, "inputVelocityModel", "Previous Velocity:");
            propertyFields3d = new StsFieldBean[]
            {
                new StsStringFieldBean(StsPreStackLineSet3d.class, "name", true, "Name"),
                new StsStringFieldBean(StsPreStackLineSet3d.class, "zDomainString", false, "Z Domain"),
                new StsBooleanFieldBean(StsPreStackLineSet3d.class, "isNMOed", "Is NMOed"),
                previousVelocityModelBean,
                new StsButtonFieldBean("AGC Properties", "Edit Automatic Gain Control Properties.", this, "displayAGCProperties"),
                new StsButtonFieldBean("Filter Properties", "Edit Filter properties.", this, "displayFilterProperties"),
                new StsButtonFieldBean("Datum Properties", "Edit Datim properties.", this, "displayDatumProperties"),
                new StsStringFieldBean(StsPreStackLineSet3d.class, "segyFilename", false, "SEGY Filename"),
                new StsStringFieldBean(StsPreStackLineSet3d.class, "segyFileDate", false, "SEGY creation date"),
                new StsIntFieldBean(StsPreStackLineSet3d.class, "nRows", false, "Number of Lines"),
                new StsIntFieldBean(StsPreStackLineSet3d.class, "nCols", false, "Number of Crosslines"),
                new StsDoubleFieldBean(StsPreStackLineSet3d.class, "xOrigin", false, "X Origin"),
                new StsDoubleFieldBean(StsPreStackLineSet3d.class, "yOrigin", false, "Y Origin"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "xInc", false, "X Inc"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "yInc", false, "Y Inc"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "zInc", false, "Z Inc"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "xMin", false, "X Loc Min"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "yMin", false, "Y Loc Min"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "zMin", false, "Min Z or T"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "zMax", false, "Max Z or T"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "angle", false, "Angle to Line Direction"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "rowNumMin", false, "Min Line"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "rowNumMax", false, "Max Line"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "colNumMin", false, "Min Crossline"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "colNumMax", false, "Max Crossline"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "dataMin", false, "Data Min"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "dataMax", false, "Data Max"),
                new StsFloatFieldBean(StsPreStackLineSet3d.class, "dataAvg", false, "Data Avg")
            };
        }
//		previousVelocityModelBean.removeAll(); // do we need this?  TJL 12/2/06
		previousVelocityModelBean.setListItems(getAvailableModelsAndVolumesList());
		if(inputVelocityModel == null)
			previousVelocityModelBean.setSelectedIndex(0);
		else
			previousVelocityModelBean.setSelectedItem(inputVelocityModel);

        return propertyFields3d;
	}

	public String getSegyFileDate()
	{
		if(segyLastModified == 0)
		{
			File segyFile = new File(segyDirectory + segyFilename);
			if(segyFile != null)
			{
				segyLastModified = segyFile.lastModified();
			}
		}
		DateFormat dateFormat = DateFormat.getDateTimeInstance();
		return dateFormat.format(new Date(segyLastModified));
	}

	public boolean anyDependencies()
	{
		StsCrossplot[] cp = (StsCrossplot[])currentModel.getCastObjectList(StsCrossplot.class);
		for(int n = 0; n < cp.length; n++)
		{
			StsSeismicBoundingBox[] volumes = cp[n].getVolumes();
			for(int j = 0; j < cp[n].volumes.getSize(); j++)
			{
				if(this == cp[n].volumes.getElement(j))
				{
					StsMessageFiles.infoMessage("Seismic PostStack3d " + getName() + " used by Crossplot " + cp[n].getName());
					return true;
				}
			}
		}
		return false;
	}

	/*
	  public void setDirCoordinate(int dirNo, float dirCoordinate)
	  {
	   if(dirNo == StsParameters.XDIR)
	   {
	 int col = getNearestColCoor(dirCoordinate);
	 if(col == currentCol) return;
	 currentCol = col;
	   }
	   else if(dirNo == StsParameters.YDIR)
	 currentRow = getNearestRowCoor(dirCoordinate);
	  }
	 */
	public boolean setDirCoordinate(int dirNo, float dirCoordinate, StsWin3dBase window)
	{
		if (lines != null)
		{
			if(dirNo == StsParameters.XDIR)
			{
				int col = getNearestColCoor(dirCoordinate);
			    setCurrentDataRowCol(window, -1, col);
			    return false;
		    }
		    else if(dirNo == StsParameters.YDIR)
            {
                int row = getNearestRowCoor(dirCoordinate);
                setCurrentDataRowCol(window, row, -1);
                return false;
			}
		}
		return true;
	}

	public int getPlaneValue(float[] xyz)
	{
		return -1;
	}
/*
	public void setCursorDisplayChanged()
	{
		cursorDisplayChanged[0] = true;
		cursorDisplayChanged[1] = true;
	}

	public boolean isCursorDisplayChanged(int dirNo)
	{
		if(!cursorDisplayChanged[dirNo])return false;
		cursorDisplayChanged[dirNo] = false;
		return true;
	}
*/

	public boolean setGLColorList(GL gl, boolean nullsFilled, byte displayType, int shader)
	{
        if(displayType == StsPreStackLineSetClass.DISPLAY_ATTRIBUTE)
            return setDisplayAttributeColorList(gl, false, shader);
        else
            return super.setGLColorList(gl, false, displayType, shader);
    }
    /** colorList is either for the seismic data or a basemap (z-slice) attribute.
	 *  If this is a z-slice (dir == ZDIR) and displayAttribute is not _SEIMIC,
	 *  then we want the displayAttribute colorList, otherwise we want the seismic colorList.
	 */
    public boolean setDisplayAttributeColorList(GL gl, boolean nullsFilled, int shader)
	{
        StsDisplaySeismicAttribute currentDisplayAttribute = getDisplayAttribute();
        if (currentDisplayAttribute == null) return false;
        return currentDisplayAttribute.setGLColorList(gl, nullsFilled, shader);
	}

	public void itemStateChanged(ItemEvent e)
	{
		if(e.getItem() instanceof StsColorscale)
			currentModel.displayIfCursor3dObjectChanged(this);
	}

	public void drawTextureTileSurface(StsTextureTile tile, GL gl, int dir, boolean is3d)
	{
		byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
		byte volumeZDomain = getZDomain();

		if(projectZDomain == StsParameters.TD_TIME)
		{
			if(volumeZDomain == StsParameters.TD_TIME)
				drawTextureTileTimeSurface(tile, gl, is3d);
		}
		else if(projectZDomain == StsParameters.TD_DEPTH)
		{
			// volumeZDomain is TD_TIME
			StsModel model = StsObject.getCurrentModel();
			if(volumeZDomain == StsParameters.TD_DEPTH) // seismic already in depth, don't need to convert so draw as if in time
				drawTextureTileTimeSurface(tile, gl, is3d);
			else
			{
				StsSeismicVelocityModel velocityVolume = model.getProject().velocityModel;
				if(velocityVolume == null)return;
				drawTextureTileDepthSurface(velocityVolume, tile, gl, dir);
			}
		}
	}

	public void drawTextureTileTimeSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		if(is3d)
			tile.drawQuadSurface3d(gl, StsProject.TD_TIME);
		else
			tile.drawQuadSurface2d(gl);
	}

	public void drawTextureTileDepthSurface(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
	{
		float cursorXInc, cursorYInc;
		if(dir == StsCursor3d.ZDIR)
		{
			return;
		}

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
		double[] xyz = tile.xyzPlane[0];
		double x1 = xyz[0];
		double y1 = xyz[1];
		double t1 = xyz[2];
		int volumeRow = this.getNearestBoundedRowCoor((float)y1);
		int volumeCol = getNearestBoundedColCoor((float)x1);
		StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
		if(velocityVolume == null)return;
		float depthMin = velocityModel.depthDatum;
		int volumeRowInc = 0;
		int volumeColInc = 0;
		if(dir == StsCursor3d.XDIR)
		{
			cursorXInc = 0;
			cursorYInc = getYInc();
			volumeRowInc = 1;
		}
		else // dirNo == StsCursor3d.YDIR
		{
			cursorXInc = getXInc();
			cursorYInc = 0;
			volumeColInc = 1;
		}
		double tInc = getZInc();
		for(int row = tile.croppedRowMin + 1; row <= tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
		{
			double x0 = x1;
			double y0 = y1;
			x1 += cursorXInc;
			y1 += cursorYInc;

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;
			double t = t1 + tile.croppedColMin * tInc;

			for(int col = tile.croppedColMin; col <= tile.croppedColMax; col++, t += tInc, colTexCoor += dColTexCoor)
			{
				float v0 = velocityVolume.getValue(volumeRow, volumeCol, col);
				float z0 = (float)(v0 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor);
				gl.glVertex3d(x0, y0, z0);
				float v1 = velocityVolume.getValue(volumeRow + volumeRowInc, volumeCol + volumeColInc, col);
				float z1 = (float)(v1 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
				gl.glVertex3d(x1, y1, z1);
			}
			gl.glEnd();
			volumeRow += volumeRowInc;
			volumeCol += volumeColInc;
		}
	}

	public boolean canDisplayZDomain()
	{
		return currentModel.getProject().canDisplayZDomain(getZDomain());
	}

	public boolean isLineInVolume(String name)
	{
		for(int j = 0; j < lines.length; j++)
		{
			if(lines[j] == null)continue;
			if(lines[j].gatherFilename.equals(name))
				return true;
		}
		return false;
	}

	/*public boolean isTextureChanged(int dirNo, int nPlane)
	{
		if(dirNo == ZDIR) return false;
		if(velocityModel == null)return false;
		return !isPlaneOK(dirNo, nPlane);
	}*/

    public ByteBuffer computeVelocityPlane(int dir, int nPlane)
    {
        float[] floatPlane;
        if(dir == ZDIR) return null;
        if(velocityModel == null) return null;
        if(!velocityModel.checkRunInterpolation()) return null;
        if(!velocityModel.hasProfiles()) return null;
        StsPreStackVelocityModel3d velocityModel3d = (StsPreStackVelocityModel3d)velocityModel;
//        byte dataType = StsFilterProperties.VELOCITY;
        int rows = nRows;
        if(dir == YDIR)
            rows = nCols;

        switch(lineSetClass.stackOption)
        {
            case STACK_NONE:
                return null;
            case STACK_NEIGHBORS:
                floatPlane = velocityModel3d.computeFloatVelocityPlane(dir, nPlane);
//                floatPlane = StsSeismicFilter.filter(dataType, floatPlane, rows, nCroppedSlices, getZInc(), filterProperties,
//                         agcProperties, velocityModel3d.dataMin, velocityModel3d.dataMax);
//                floatPlane = getNeighborBuffer(dir, nPlane, floatPlane);
                break;
            case STACK_LINES:
                floatPlane = velocityModel3d.computeFloatVelocityPlane(dir, nPlane);
//                floatPlane = StsSeismicFilter.filter(dataType, floatPlane, rows, nCroppedSlices, getZInc(), filterProperties,
//                         agcProperties, velocityModel3d.dataMin, velocityModel3d.dataMax);
                break;
            default:
                return null;
        }
        
        // Convert to byte buffer
        if (floatPlane == null) return null;
        boolean autoSaturate = semblanceRangeProperties.autoSaturateColors;
        double dataMin = velocityModel3d.dataMin;
        double dataMax = velocityModel3d.dataMax;
        if (autoSaturate)
        {
            float[] minMax = StsMath.minMax(floatPlane);
            dataMin = minMax[0];
            dataMax = minMax[1];
            float[] otherPlane = null;
            if (dir == YDIR)
            {
                int nOtherPlane = currentModel.win3d.cursor3d.getCurrentGridCoordinate(XDIR);
                otherPlane = velocityModel3d.computeFloatVelocityPlane(XDIR, nOtherPlane);
            }
            else
            {
                int nOtherPlane = currentModel.win3d.cursor3d.getCurrentGridCoordinate(YDIR);
                otherPlane = velocityModel3d.computeFloatVelocityPlane(YDIR, nOtherPlane);
            }
            if (otherPlane != null)
            {
                minMax = StsMath.minMax(otherPlane);  //we want to use the same scalar for both planes!!
                dataMin = Math.min(minMax[0], dataMin);
                dataMax = Math.max(minMax[1], dataMax);
            }
            velocityModel3d.setColorScaleRange((float)dataMin, (float)dataMax, false); //this should update velocity colorbar to show new velocity ranges
        }
        double scale = 253.0/(dataMax - dataMin);  //using 254 was causing the max velocity to come out transparent - very obvious in interval-velocity mode SWC 11/17/09
        ByteBuffer buffer = ByteBuffer.allocateDirect(rows * nSlices);
        for(int n = 0; n < nSlices*rows; n++)
        {
            byte velocity = StsMath.unsignedIntToUnsignedByte((int)((floatPlane[n] - dataMin) * scale));
            buffer.put(velocity);
        }
        buffer.rewind();
        return buffer;
    }

    public StsPreStackVolume getPreStackVolume()
    {
        byte displayType = lineSetClass.getDisplayType();
        return outputVolumes[displayType];
    }

    public StsPreStackVolume getPreStackVolume(byte displayType)
    {
        if(outputVolumes == null) return null;
        return outputVolumes[displayType];
    }

    /** Prestack volumes (velocity, stack, semblance) can always be isVisible on X or Y (vertical) planes.
     *  If they haven't been computed for these planes, they will be computed in this display pass.
     *  In this case, we need only check that the dirCoordinate is legitimate.
     *  Stack and semblance can't be isVisible on Z (horizontal) plane unless volumes are already computed.
     *  This can be checked by the status of the planes: volume.isPlaneOk(dir, coor).
     *  If not ok, the displayType will be set to ATTRIBUTE if it is a z plane, otherwise set to NONE.
     */
    public boolean isVolumePlaneOk(byte displayType, int dirNo, float dirCoordinate)
    {
        if(displayType < 0 || displayType >= nVolumeTypes) return false;
        StsPreStackVolume volume = getPreStackVolume(displayType);
        if(volume == null) return false;
        int nPlane = getCursorPlaneIndex(dirNo, dirCoordinate);
        if(nPlane == -1) return false;
        if(dirNo != ZDIR) return true;
        return volume.isPlaneOK(dirNo, nPlane);
    }

     public void computePreStackVolume(byte displayType)
	 {
         outputVolumes[displayType].checkComputePreStackVolume();
	 }

	 public int getGatherIndex(int row, int col)
	 {
		 if (isInline)
			 return col + row * nCols;
		 else
			 return row + col * nRows;
	 }

	 public ByteBuffer computeSemblancePlane(int dir, int nPlane)
	 {
         StsPreStackVolume semblanceVolume = getSemblanceVolume();
         switch (lineSetClass.stackOption)
		 {
			 case STACK_NONE:
				 return null;
			 case STACK_NEIGHBORS:
				 return semblanceVolume.getByteData(dir, nPlane, true);
			 case STACK_LINES:
				 return semblanceVolume.getByteData(dir, nPlane, false);
			 default:
				 return null;
		 }
	 }

	 public ByteBuffer computeStackedPlane(int dir, int nPlane)
	 {
         StsPreStackVolume stackedVolume = getStackedVolume();
         switch (lineSetClass.stackOption)
		 {
			 case STACK_NONE:
				 return null;
			 case STACK_NEIGHBORS:
				 return stackedVolume.getByteData(dir, nPlane, true);
			 case STACK_LINES:
				 return stackedVolume.getByteData(dir, nPlane, false);
			 default:
				 return null;
		 }
	 }

	 public void objectPropertiesChanged(Object object)
	 {
		if(!agcOrFilterPropertiesChanged(object)) return;
        setVolumesPlaneOKFlagsFalse();
//        currentModel.viewObjectChanged(this);
     }

	public void objectChanged()
	{
		currentModel.changeCursor3dObject(this);
	}

    private float[] getNeighborBuffer(int dir, int nPlane, float[] lineBuffer)
    {
        int row = -1, col = -1;
        try
        {
            checkTransparentTrace(nSlices);
            float[] floats = new float[nSlices];
            boolean[][] isNeighbor = velocityModel.interpolation.isNeighbor;
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
			boolean[][] isNeighbor = velocityModel.interpolation.isNeighbor;
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

	public boolean isPlaneOK(byte displayType, int dir, int nPlane)
	{
        if(outputVolumes == null) return false;
        return outputVolumes[displayType].isPlaneOK(dir, nPlane);
	}

	public void setPlaneFalse(int row, int col)
	{
        for(int n = 0; n < nVolumeTypes; n++)
        {
            outputVolumes[n].setPlaneOK(XDIR, col, false);
            outputVolumes[n].setPlaneOK(YDIR, row, false);
        }
    }

    /** When this texture is no longer needed, unlock the corresponding mapBuffer plane so it can be cleared if
     *  additional space is needed.
     *  @param dir plane direction
     *  @param dirCoor plane coordinate
     */
    public void unlockPlanes(int dir, float dirCoor)
	{
        if(outputVolumes == null) return;
        int nPlane = getCursorPlaneIndex(dir, dirCoor);
        for(int n = 0; n < nVolumeTypes; n++)
            outputVolumes[n].unlockPlane(dir, nPlane);
    }

    public void setVolumesPlaneOKFlagsFalse()
	{
        if(outputVolumes == null) return;
        for(int n = 0; n < nVolumeTypes; n++)
            outputVolumes[n].setAllPlaneOKFlagsFalse();
    }

    public void updatePreStackVolumes(StsRadialInterpolation interpolation)
    {
       for(int row = 0; row < nRows; row++)
       {
           boolean rowChanged = interpolation.rowChanged[row];
           if(rowChanged)
           {
               for(int n = 0; n < nVolumeTypes; n++)
                    outputVolumes[n].setPlaneOK(YDIR, row, false);
           }
       }
       for(int col = 0; col < nCols; col++)
       {
           boolean colChanged = interpolation.colChanged[col];
           if(colChanged)
           {
               for(int n = 0; n < nVolumeTypes; n++)
                    outputVolumes[n].setPlaneOK(XDIR, col, false);
           }
       }
       boolean sliceChanged = interpolation.changed;
       if(sliceChanged)
       {
            for(int n = 0; n < nVolumeTypes; n++)
                outputVolumes[n].setDirPlanesOK(ZDIR, false);
       }
    }

    public void close()
    {
       boolean changed = false;
       if(outputVolumes == null) return;
       for(int n = 0; n < nVolumeTypes; n++)
            if(outputVolumes[n].closeChanged()) changed = true;
//       if(changed) dbFieldChanged("outputVolumeTimes", outputVolumeTimes);
    }
/*
	public void setPlaneOK(boolean ok)
	{
		for(int col = 0; col < nCols; col++)
		{
			stackedVolume.setPlaneOK(XDIR, col, ok);
			semblanceVolume.setPlaneOK(XDIR, col, ok);
		}
		for(int row = 0; row < nRows; row++)
		{
			stackedVolume.setPlaneOK(YDIR, row, ok);
			semblanceVolume.setPlaneOK(YDIR, row, ok);
		}
	}
*/
/*
    public FloatBuffer readRowPlaneFloatBufferData(int nPlane)
	{
		return fileMapRowFloatBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY).asFloatBuffer();
	}

	public ByteBuffer readRowPlaneByteBufferData(int nPlane)
	{
		return fileMapRowFloatBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY);
	}
*/
	public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, boolean axesFlipped)
	{
        if(lineSetClass.getShowVelStat() == true && velStat != null)
			velStat.drawOnCursor2d(glPanel3d, dirNo, axesFlipped);
	}
 /*
	public boolean hasBasemapChanged()
	{
		return cursorDisplayChanged[ZDIR];
	}

	public void setBasemapChanged(boolean value)
	{
		cursorDisplayChanged[ZDIR] = value;
	}
*/
    private boolean hasHandVels( Object[] preStackLineSet)
    {
        for( int pset=0; pset<preStackLineSet.length; pset++)
        {
            StsPreStackLineSet3d pLineSet = (StsPreStackLineSet3d)preStackLineSet[pset];
            if( pLineSet.getHandVelName() != null)
                return true;
        }
        return false;
    }

    public Object[] getAvailableModelsAndVolumesList()
	{
		Object[] objects;

        Object[] models = currentModel.getTrimmedList(StsPreStackVelocityModel3d.class);
		Object[] volumes = (Object[]) currentModel.getTrimmedList(StsSeismicVolume.class);
        objects = new Object[] {NO_MODEL};
        int nItems = models.length + volumes.length;
        if( nItems > 0)
        {
            objects = (Object[]) StsMath.arrayAddArray(objects, volumes);
			objects = (Object[]) StsMath.arrayAddArray(objects, models);
		}

        Object[] preStackLineSet = (Object[] )currentModel.getTrimmedList( StsPreStackLineSet3d.class);
        if( hasHandVels( preStackLineSet))
            objects = (Object[] )StsMath.arrayInsertElementBefore( objects, HAND_VEL, 1);

        return objects;
	}

	public Object[] getAvailableModelsList()
	{
		Object[] models = (Object[]) currentModel.getTrimmedList(StsPreStackVelocityModel3d.class);
		if (models.length == 0)
		{
			models = new Object[]
				{NO_MODEL};
		}
		return models;
	}

    /** PostStack3d is not NMOed.  If we have interpreted profiles, use them to NMO (flatten) gather;
	 *  Otherwise if we have an inputVelocityModel (typically a segy volume has been read), then use it to NMO.
	 */
	public float[] getNMOVelocities(int row, int col)
	{
		if(velocityModel == null)return null;
		float[] velocities = velocityModel.getVelocities(row, col);
		if(velocities != null)return velocities;
		if(inputVelocityModel == null)return null;
		return inputVelocityModel.getVolumeVelocities(row, col);
	}

	/** PostStack3d is NMOed.  If  we have an inputVelocityModel (typically a segy volume has been read), then use it to DNMO (unflatten) gather). */
	public float[] getDNMOVelocities(int row, int col)
	{
		if(inputVelocityModel == null)return null;
		float[] velocities = inputVelocityModel.getVelocities(row, col);
		if(velocities != null)return velocities;
		return inputVelocityModel.getVolumeVelocities(row, col);
	}

	public int compareProfileLocations(StsVelocityProfile profile0, StsVelocityProfile profile1)
	{
		if(isInline)
			return nCols * (profile0.row - profile1.row) + profile0.col - profile1.col;
		else
			return nRows * (profile0.col - profile1.col) + profile0.row - profile1.row;
	}

	public int getVolumeRowColIndex(int row , int col)
	{
		row = StsMath.minMax(row , 0 , nRows - 1);
		col = StsMath.minMax(col , 0 , nCols - 1);
		if (isInline)
		{
			return nCols * row + col;
		}
		else
		{
			return nRows * col + row;
		}
	}

	public String getGatherDescription(int row, int col, int nGatherTraces)
	{
		StsPreStackLine line = getDataLine(row, col);
		if (line == null)
		{
			return getEmptyGatherDescription(row, col);
		}
		else
		{
			return getFullGatherDescription(row, col, nGatherTraces);
		}
	}

    public String getFullGatherDescription(int row, int col, int nTraces)
    {
        if(isInline)
            return new String(getName() + " inline " + getRowNumFromRow(row) + " crossline " + getColNumFromCol(col) + " traces " + nTraces);
        else
			return new String(getName() + " inline " +  getColNumFromCol(col) + " crossline " + getRowNumFromRow(row) + " traces " + nTraces);
    }

    public String getEmptyGatherDescription(int row, int col)
    {
        if(isInline)
            return new String(getName() + " inline " + getRowNumFromRow(row) + " crossline " + getColNumFromCol(col) + " traces 0");
        else
			return new String(getName() + " inline " + getColNumFromCol(col) + " crossline " + getRowNumFromRow(row) + " traces 0");
    }

    public boolean getIsVisibleOnCursor()
	{
		return isVisible && ((StsPreStackLineSet3dClass)lineSetClass).getIsVisibleOnCursor();
	}

	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		 exportPanel = new StsSeismicExportPanel3d(currentModel, this, "Export Pre-Stack Seismic 3d", true);
		 return exportPanel;
	}

    public Class[] getViewClasses()
    {
        Class[] viewClasses = new Class[] { StsViewGather3d.class, StsViewSemblance3d.class, StsViewCVStacks.class, StsViewVVStacks.class };
        if(velocityModel != null && velocityModel.hasProfiles())
            viewClasses = (Class[])StsMath.arrayAddArray(viewClasses, new Class[] {  StsViewResidualSemblance3d.class } );
        return viewClasses;
    }

    static public void main(String[] args)
	{
		StsPreStackLineSet3d lineSet = new StsPreStackLineSet3d();
//		lineSet.analysisRowInc = 3;
//		lineSet.analysisColInc = 3;
		lineSet.nRows = 5;
		lineSet.nCols = 12;
		lineSet.lines = new StsPreStackLine[5];
		for (int n = 0; n < 5; n++)
		{
			StsPreStackLine3d line = new StsPreStackLine3d();
			lineSet.lines[n] = line;
		}
		int[] rowCol = new int[]
			{0, 0};
		for (int i = 0; i < 100; i++)
		{
			rowCol = lineSet.getNextRowCol(rowCol[0], rowCol[1]);
			System.out.println("next: " + rowCol[0] + " " + rowCol[1]);
		}
		rowCol = new int[]
			{0, 0};
		for (int i = 0; i < 100; i++)
		{
			rowCol = lineSet.getPrevRowCol(rowCol[0], rowCol[1]);
			System.out.println("prev: " + rowCol[0] + " " + rowCol[1]);
		}
	}
}

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
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 8, 2006
 * Time: 8:30:45 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackLineSet2d extends StsPreStackLineSet implements StsTreeObjectI, StsViewable
{
    /** temporary mapped fileBlocks where stacked byte display data is stored */
    transient public StsCubeFileBlocks fileMapBlocks2d;

	static private StsComboBoxFieldBean displayCultureBean;
	static public StsFieldBean[] displayFields2d = null;
	static public StsFieldBean[] propertyFields2d = null;

    public String getGroupname()
    {
        return StsSeismicBoundingBox.group2dPrestack;
    }

    public StsPreStackLineSet2d()
	{
//        System.out.println("new StsPreStackLineSet2d()");
    }

	public StsPreStackLineSet2d(boolean persistent)
	{
		super(persistent);
	}

	static public StsPreStackLineSet2d constructor(String name, StsFile[] files, StsModel model, StsProgressPanel panel, StsRotatedGridBoundingBox boundingBox)
	{
		try
		{
			return new StsPreStackLineSet2d(name, files, model, panel, boundingBox);
		}
		catch (FileNotFoundException fnfe)
		{
			return null;
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet2d.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}

	private StsPreStackLineSet2d(String name, StsFile[] files, StsModel model, StsProgressPanel panel, StsRotatedGridBoundingBox boundingBox) throws
		FileNotFoundException, StsException
	{
		super(false);

		setName(name);
//		toggleZDomain(StsProject.TD_TIME);
		int nFiles = files.length;
        for (int n = 0; n < nFiles; n++)
		{
			if (!files[n].exists())
			{
				panel.appendLine("\tFile not found " + files[n].getPathname());
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to find file " + files[n].getPathname() + ".");
				throw new FileNotFoundException();
			}
		}
//		preStackLines2d = new ArrayList(10);
//		preStackLines = StsObjectRefList.constructor(nFiles, 2, "preStackLines", this, true);
		boolean originSet = false;
		addUnrotatedBoundingBox(boundingBox);
		panel.initialize(files.length);
        stsDirectory = files[0].getDirectory();
        for (int n = 0; n < nFiles; n++)
		{
			StsPreStackLine line = StsPreStackLine2d.constructor(files[n], n, model, this);
			if (line == null)
			{
				continue;
			}
			if (!originSet)
			{
				checkSetOriginAndAngle(line.xOrigin, line.yOrigin, line.angle);
				isNMOed = line.isNMOed;

				byte zDomainByte = StsParameters.getZDomainFromString(line.zDomain);
				if (zDomainByte == StsParameters.TD_NONE)
				{
					zDomain = StsParameters.TD_TIME_STRING;
				}
				else
				{
					zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
				}
				originSet = true;
			}
			panel.appendLine("Loading line " + line.getName() + ": " + (n+1) + " of " + nFiles);

			// TODO:
			// dataMin & dataMax could be different for volume if lines were processed separately.
			// We really can't allow this as it makes the texture mapping enormously more expensive
			// to translate each unsigned byte on a byte-by-byte basis to the correct index into the
			// colorscale.  In this case, we are better off reading in floats, scaling them to bytes and
			// putting them in the texture array.
			dataMin = Math.min(dataMin, line.dataMin);
			dataMax = Math.max(dataMax, line.dataMax);
            lines = (StsPreStackLine[]) StsMath.arrayAddElement(lines, line);

            addUnrotatedBoundingBox(line);

			// haque until we figure this out....
			if (line.yInc < 0.0f)
			{
				line.yInc = -line.yInc;
			}

			setTraceOffsetRange(line);
            panel.setValue(n+1);
            panel.setDescription("Loaded line : " + (n+1) + " of " + nFiles);
            panel.appendLine("Completed loading line index:" + line.getIndex() + " from file " + line.getName());
//			checkMemoryStatus("after loading line" + n);
		}

		stsDirectory = files[0].getDirectory();
		isVisible = true;
		//       setVelStat(lineSet.getShowVelStat());

		setVolumeRange();

        analysisRowInc = 1;

//		sortLines2d();
		addToolbar();
		if (!addToProject())
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add volume " + name + " to project.");
			throw new StsException(StsException.WARNING, "Failed to add prestack 2d lines to project.");
		}
		setLinesXYs();
        currentModel.getProject().checkCursor3d();
        currentModel.getProject().setCursorDisplayXYAndGridCheckbox(true);
        initialize();
	}

    public StsPreStackVelocityModel constructVelocityModel()
    {
        return new StsPreStackVelocityModel2d(this);
    }

    public int getNColsForRow(int row)
    {
        if(row >= nRows) return 0;
        return lines[row].nCols;
    }

    private void setVolumeRange()
	{
		nRows = lines.length;
		for(int n = 0; n < nRows; n++)
        {
            nCols = Math.max(nCols, lines[n].nCols);
            nSlices = Math.max(nSlices, lines[n].nSlices);
        }
    }

	private void setLinesXYs()
	{
		for(int n = 0; n < lines.length; n++)
			((StsPreStackLine2d)lines[n]).setLineXYs();
    }

    public boolean textureChanged()
    {
		for(int n = 0; n < lines.length; n++)
			((StsPreStackLine2d)lines[n]).textureChanged();
        return true;
    }

    public boolean addToProject()
	{
//		adjustBoundingBox(100);
//		adjustRowColNumbering(this);

		if (!currentModel.getProject().addToProject(this, false))
		{
			return false;
		}
        return true;
	}

	public void setInputVelocityModel(Object object)
	{
		if (object instanceof StsPreStackVelocityModel2d)
		{
			dbFieldChanged("inputVelocityModel", inputVelocityModel, object);
		}
		else
		{
			dbFieldChanged("inputVelocityModel", inputVelocityModel, null);
		}
	}

    public boolean allocateVolumes(String mode, boolean loadFromFile)
	{
		if (fileMapBlocks2d != null) return true;

		try
		{
			StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
			setName(currentModel.getName() + "." + getName());
            // file used for mapped stacked data
            setRowCubeFilename(StsSeismicBoundingBox.group2dPrestack + ".stack." + getName());
			fileMapBlocks2d = new StsCubeFileBlocks(YDIR, nRows, nCols, nSlices, this, rowCubeFilename, "rw", 1, memoryManager);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolume.allocateVolumes) failed.", e, StsException.WARNING);
			return false;
		}
    }

    public int[] getRowColFromCoors(float x, float y)
    {
        float dSq = StsParameters.largeFloat;
        int nearestRow = -1;
        int nearestCol = -1;
        for(int row = 0; row < lines.length; row++)
        {
            StsPreStackLine line = lines[row];
            int nCols = line.nCols;
            float[] cdpX = line.cdpX;
            float[] cdpY = line.cdpY;
            for (int col = 0; col < nCols; col++)
			{
				float dx = (cdpX[col] - x);
                float dy = (cdpY[col] - y);
                float dSqNew = dx*dx + dy*dy;
                if(dSqNew < dSq)
                {
                    nearestRow = row;
                    nearestCol = col;
                    dSq = dSqNew;
                }
            }
        }
        return new int[] { nearestRow, nearestCol } ;
    }

    public int getVolumeRowColIndex(int row, int col)
	{
		row = StsMath.minMax(row, 0, nRows - 1);
		col = StsMath.minMax(col, 0, nCols - 1);
		return nCols * row + col;
	}

	/** for prestack 2d, each line is a row and col is gather position along the line  from 0 to nGathers-1 */
	public StsPreStackLine getDataLine(int row, int col)
	{
		if (lines == null)
		{
			return null;
		}

		int nLine = row;
		int nGather = col;

		int nLines = lines.length;
		if (nLine < 0 || nLine >= nLines)
		{
			return null;
		}

		StsPreStackLine line = lines[nLine];
		int compare = line.compareTo(nGather);
		if (compare == 0)
			return line;
		else
			return null;
	}

    public double[][] getLineXYCDPs(int row)
	{
        StsPreStackLine2d line = getLine(row);
        if (line == null) return null;
		return line.getLineXYCDPs();
	}

    private StsPreStackLine2d getLine(int row)
	{
 		if (lines == null) return null;
		if(row < 0 || row >= lines.length) return null;
		return (StsPreStackLine2d)lines[row];
	}

    public int getNSlices(int row)
	{
		StsPreStackLine2d line = getLine(row);
        if(line == null) return 0;
        return line.getNSlices();
    }

    public float getZMin(int row)
	{
		StsPreStackLine2d line = getLine(row);
        if(line == null) return 0;
        return line.getZMin();
    }

    public float getZInc(int row)
	{
		StsPreStackLine2d line = getLine(row);
        if(line == null) return 0;
        return line.getZInc();
    }

    public StsSeismicBoundingBox getLineBoundingBox(int row)
	{
		return getLine(row);
	}

   /*  JKF 5/17/07
	public void setPlaneOK(boolean ok)
	{
		if (lines == null)
		{
			return;
		}
		for (int n = 0; n < lines.length; n++)
		{
			( (StsPreStackLine2d) lines[n]).setPlaneOK(ok);
		}
	}*/

   /** Overrides method in StsRotatedGridBoundingBox.
     *  In 2d, row is the line number and col is the gather number
     */
    public float[] getXYCoor(int nLine, int nGather)
    {
        float x = ((StsPreStackLine2d)lines[nLine]).cdpX[nGather];
        float y = ((StsPreStackLine2d)lines[nLine]).cdpY[nGather];
        return new float[] { x, y };
    }

    /** Overrides method in StsRotatedGridBoundingBox.
     *  In 2d, row is the line number and col is the gather number
     */
    public void getXYCoors(int nLine, int nGather, float[] coors)
    {
        try
        {
            coors[0] = ((StsPreStackLine2d)lines[nLine]).cdpX[nGather];
            coors[1] = ((StsPreStackLine2d)lines[nLine]).cdpY[nGather];
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getXYCoors", " nLine: " + nLine + " nGather: " + nGather);
        }
    }


    public float[] getXYCoors(int nLine, int nGather)
    {
        float[] xy = new float[2];
        xy[0] = ((StsPreStackLine2d)lines[nLine]).cdpX[nGather];
        xy[1] = ((StsPreStackLine2d)lines[nLine]).cdpY[nGather];
        return xy;
    }

    /** Overrides method in StsRotatedGridBoundingBox.
     *  In 2d, row is the line number and col is the gather number
     */
    public float getXCoor(int nLine, int nGather)
    {
        return ((StsPreStackLine2d)lines[nLine]).cdpX[nGather];
    }
    /** Overrides method in StsRotatedGridBoundingBox.
     *  In 2d, row is the line number and col is the gather number
     */
    public float getYCoor(int nLine, int nGather)
    {
        return ((StsPreStackLine2d)lines[nLine]).cdpY[nGather];
    }

    /** For a 2d lineset, each line is a row and each row has a different number of columns. */
    public int[] adjustLimitRowCol(int nRow, int nCol)
    {
        nRow = StsMath.minMax(nRow, 0, nRows-1);
        nCol = StsMath.minMax(nCol, 0, lines[nRow].nCols-1);
        return new int[] { nRow, nCol };
    }

    public int adjustLimitCol(int col)
    {
        return StsMath.minMax(col, 0, nCols-1);
    }

    private boolean hasHandVels( Object[] preStackLineSet)
    {
        for( int pset=0; pset<preStackLineSet.length; pset++)
        {
            StsPreStackLineSet pLineSet = (StsPreStackLineSet2d)preStackLineSet[pset];
            for( int pline=0; pline<pLineSet.lines.length; pline++)
            {
                if( pLineSet.lines[pline].getHandVelName() != null)
                    return true;
            }
        }
        return false;
    }

    // TODO what can we do with volumes in 2d...just models?
    public Object[] getAvailableModelsAndVolumesList()
	{
		Object[] objects = new Object[] {NO_MODEL};

        Object[] lineSet = (Object[] )currentModel.getTrimmedList( StsSeismicLineSet.class);
        Object[] models  = (Object[] )currentModel.getTrimmedList( StsPreStackVelocityModel2d.class);
		Object[] volumes = (Object[] )currentModel.getTrimmedList( StsSeismicVolume.class);
        int nItems = lineSet.length + models.length + volumes.length;
        if (nItems > 0)
		{
            objects = (Object[] )StsMath.arrayAddArray( objects, lineSet);
            objects = (Object[] )StsMath.arrayAddArray( objects, models);
            objects = (Object[] )StsMath.arrayAddArray( objects, volumes);
        }

        Object[] preStackLineSet = (Object[] )currentModel.getTrimmedList( StsPreStackLineSet2d.class);
        if( hasHandVels( preStackLineSet))
            objects = (Object[] )StsMath.arrayInsertElementBefore( objects, HAND_VEL, 1);
        
        return objects;
	}

    public Object[] getAvailableModelsList()
	{
		Object[] models = currentModel.getTrimmedList(StsPreStackVelocityModel.class);
		if (models.length == 0)
		{
			models = new Object[] {NO_MODEL};
		}
		return models;
	}

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        if(!mouse.isButtonStateReleased(StsMouse.LEFT)) return false;
		// compute world coor line from front to back of view box
		double[][] mouseLine = ((StsGLPanel3d)glPanel).getViewLineAtMouse(mouse);
		// clip line to project boundingBox
		clipLine(mouseLine[0], mouseLine[1]);
		// for each line, compute intersection point and fraction along mouseLine (to determine closest)
		double[]  nearestXYZF = new double[4];
		double[] intersectionXYZF = new double[4];
		StsPreStackLine2d line = (StsPreStackLine2d)lines[0];
		int nNearestGather = line.getMouseSeismicLine2dIntersect(mouseLine, nearestXYZF);
		int nNearestLine = 0;
		for(int n = 1; n < lines.length; n++)
		{
			line = (StsPreStackLine2d)lines[n];
			int newNearestGather = line.getMouseSeismicLine2dIntersect(mouseLine, intersectionXYZF);
			if(newNearestGather >= 0 && intersectionXYZF[3] < nearestXYZF[3])
			{
				nNearestGather = newNearestGather;
				System.arraycopy(intersectionXYZF, 0, nearestXYZF, 0, 4);
				nNearestLine = n;
			}


		}
		if(nNearestGather < 0) return false;
		jumpToRowCol(new int[] { nNearestLine, nNearestGather }, ((StsGLPanel3d)glPanel).window);
		return true;
	}

    public void setCurrentLineTextureChanged()
    {
        if(currentLine != null) ((StsPreStackLine2d)currentLine).textureChanged();
    }

    public void display(StsGLPanel3d glPanel3d)
	{
        if(!isVisible) return;

        GL gl = glPanel3d.getGL();
		if (lines != null && lines.length > 0)
		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glDisable(GL.GL_BLEND);
			gl.glShadeModel(GL.GL_FLAT);

//			checkMemoryStatus("before display 2d lines");
			for (int i = 0; i < lines.length; i++)
			{
				( (StsPreStackLine2d) lines[i]).display2dLine(glPanel3d, gl);
			}
            drawGatherIn3d(glPanel3d);
//			checkMemoryStatus("after display 2d lines");
			gl.glEnable(GL.GL_LIGHTING);
		}
	}

	private void drawGatherIn3d(StsGLPanel3d glPanel3d)
	{
        StsSuperGather superGather = getSuperGather(glPanel3d.window);
        if(superGather.superGatherRow < 0) return;
        StsPreStackLine2d line2d = (StsPreStackLine2d)lines[superGather.superGatherRow];
        line2d.drawVerticalLineAtCDP(glPanel3d.getGL(), superGather.superGatherCol);
	}

    public String getEmptyGatherDescription(int row, int col)
    {
        return new String(stemname + " line " + row + " gather number " + col + " traces 0");
    }

    public String getFullGatherDescription(int row, int col, int nTraces)
    {
        return new String(stemname + " CDP " + getCDP(row, col) + " line " + row + " gather " + col + " traces " + nTraces);
    }

    
    public StsFieldBean[] getDisplayFields()
	{
		try
		{
			if (displayFields2d == null)
			{
				colorscalesBean = new StsComboBoxFieldBean(StsPreStackLineSet2d.class, "currentColorscale", "Colorscales");
				colorscaleBean = new StsEditableColorscaleFieldBean(StsPreStackLineSet2d.class, "currentColorscale");
				displayCultureBean = new StsComboBoxFieldBean(StsPreStackLineSet2d.class, "displayCulture", "S/R Attribute");

                StsIntFieldBean traceThresholdBean = new StsIntFieldBean(StsPreStackLineSet2d.class, "traceThreshold", 0, 100, "Minimum Traces Required:", true);
                traceThresholdBean.setToolTipText("Minimum number of traces required in a gather for velocity analysis");

                displayFields2d = new StsFieldBean[]
				{
					new StsBooleanFieldBean(StsPreStackLineSet2d.class, "isVisible", "Enable"),
					new StsBooleanFieldBean(StsPreStackLineSet2d.class, "readoutEnabled", "Mouse Readout"),
					new StsBooleanFieldBean(StsPreStackLineSet2d.class, "showSource", "Show Sources"),
					new StsBooleanFieldBean(StsPreStackLineSet2d.class, "showReceiver", "Show Receivers"),
                    new StsIntFieldBean(StsPreStackLineSet2d.class, "analysisColInc", 0, 100, "Col Analysis Increment:", true),
                    traceThresholdBean,
                    colorscalesBean = new StsComboBoxFieldBean(StsPreStackLineSet.class, "currentColorscale", "Colorscales", "colorscaleList"),
					colorscaleBean,
					displayCultureBean,
					new StsButtonFieldBean("Semblance Display Properties", null, StsPreStackLineSet.class, "semblanceDisplayPropertiesDialog"),
					new StsButtonFieldBean("Gather Display Properties", null, StsPreStackLineSet.class, "gatherDisplayPropertiesDialog")
				};
			}
			// colorscalesBean.setListItems(colorscales.getElements());
			// colorscalesBean.setSelectedItem(currentColorscale);
			// colorscaleBean.setValueObject(currentColorscale);
			// colorscaleBean.setHistogram(dataHist);
			displayCultureBean.setListItems(displayCultures);
			return displayFields2d;
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet2d.getDisplayFields() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public StsFieldBean[] getPropertyFields()
	{
		try
		{
			if (propertyFields2d == null)
			{
				propertyFields2d = new StsFieldBean[]
					{
					new StsStringFieldBean(StsPreStackLineSet2d.class, "name", true, "Name"),
					new StsStringFieldBean(StsPreStackLineSet2d.class, "zDomainString", false, "Z Domain"),
					new StsBooleanFieldBean(StsPreStackLineSet2d.class, "isNMOed", "Is NMOed"),
					new StsStringFieldBean(StsPreStackLineSet2d.class, "segyFilename", false, "SEGY Filename"),
					new StsStringFieldBean(StsPreStackLineSet2d.class, "segyFileDate", false, "SEGY creation date"),
					new StsDoubleFieldBean(StsPreStackLineSet2d.class, "xOrigin", false, "X Origin"),
					new StsDoubleFieldBean(StsPreStackLineSet2d.class, "yOrigin", false, "Y Origin"),
					new StsFloatFieldBean(StsPreStackLineSet2d.class, "zMin", false, "Min Z or T"),
					new StsFloatFieldBean(StsPreStackLineSet2d.class, "zMax", false, "Max Z or T"),
					new StsFloatFieldBean(StsPreStackLineSet2d.class, "dataMin", false, "Data Min"),
					new StsFloatFieldBean(StsPreStackLineSet2d.class, "dataMax", false, "Data Max"),
					new StsFloatFieldBean(StsPreStackLineSet2d.class, "dataAvg", false, "Data Avg"),
					new StsButtonFieldBean("Semblance AGC & filter properties", null, StsPreStackLineSet.class, "semblanceAgcAndFilterPropertiesDialog"),
					new StsButtonFieldBean("Gather AGC & filter properties", null, StsPreStackLineSet.class, "gatherAgcAndFilterPropertiesDialog")
				};
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet2d.getPropertyFields() failed.", e, StsException.WARNING);
			return null;
		}
		return propertyFields2d;
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		 exportPanel = new StsSeismicExportPanel2d(currentModel, this, "Export Pre-Stack Seismic 2d", true);
		 return exportPanel;
	}

    public void updatePreStackVolumes(StsRadialInterpolation interpolation)
    {
        
    }

    static public void main(String[] args)
	{
		StsPreStackLineSet2d lineSet = new StsPreStackLineSet2d();
//		lineSet.analysisRowInc = 3;
//		lineSet.analysisColInc = 3;
		lineSet.lines = new StsPreStackLine[5];
		for (int n = 0; n < 5; n++)
		{
			StsPreStackLine2d line = new StsPreStackLine2d();
			lineSet.lines[n] = line;
			line.nCols = n + 10;
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

	// TODO implement this  JKF 5/17/07
	public boolean isPlaneOK(byte displayType, int dir, int nPlane)
	{
		return false;
	}

	public void computePreStackVolume(byte displayType)
	{
	}

    public Class[] getViewClasses()
    {
        Class[] viewClasses = new Class[] { StsViewGather2d.class, StsViewSemblance2d.class, StsViewCVStacks.class, StsViewVVStacks.class };
        if(velocityModel != null && velocityModel.hasProfiles())
            viewClasses = (Class[])StsMath.arrayAddArray(viewClasses, new Class[] {  StsViewResidualSemblance2d.class } );
        return viewClasses;
    }

    /*	public StsGather2d setCurrentLineAndGather2d(StsGLPanel3d glPanel3d)
	 {
	  StsGather2d gather = getGather2d(glPanel3d);
	  if(lines2d == null)return null;
	  currentLine = lines2d[0];
	  currentModel.viewObjectChanged(this);
	  if(currentLine == null)return null;
//		setCurrentVelocityProfile(gather.row, gather.col);
	  return gather;
	 }

	 public StsGather2d getGather2d(StsGLPanel3d glPanel3d)
	 {
	  // if this is a parentWindow, it will be in the gatherList
	  for(int n = 0; n < gathers2d.length; n++)
	   if(gathers2d[n].glPanel3d == glPanel3d)return gathers2d[n];
	  // Might be an auxilary window: get gather from parentWindow
	  StsWin3dBase parentWindow = glPanel3d.window.parentWindow;
	  if(parentWindow != null && glPanel3d.window != parentWindow)
	   return getGather2d(parentWindow.glPanel3d);
	  // need to create a gather
	  StsGather2d gather = new StsGather2d(currentLine, glPanel3d);
	  gathers2d = (StsGather2d[])StsMath.arrayAddElement(gathers2d, gather);
	  return gather;
	 }

		public boolean setGLColorList(GL gl, boolean nullsFilled, int dirNo, boolean displayVelocity, int shader)
	 {
			StsColorList colorList = getSeismicColorList(displayVelocity);
			if(colorList == null) return false;
			colorList.setGLColorList( gl, nullsFilled, shader);
			return true;
		}

		private StsColorList getSeismicColorList(boolean displayVelocity)
	 {
	  if(displayVelocity)
	   return velocityColorList;
	  else
	   return seismicColorList;
	 }*/
    

    public RowColIterator getRowColIterator(int row, int col, int direction)
    {
        if(rowColIterator == null)
            rowColIterator = new RowColIterator2d(row, col, direction);
        else
            rowColIterator.initialize(row, col, direction);
        return rowColIterator;
    }
    
    public class RowColIterator2d extends com.Sts.DBTypes.StsPreStackLineSet.RowColIterator {
    	
    	int lineColInc;   //get these from initialization!!! SWC 8/4/09
    	int lineColStart;
    	int rowStart = 0;
    	int rowInc = 1;

    	RowColIterator2d(StsPreStackLineSet stsPreStackLineSet, int row, int col,int direction) {
    		stsPreStackLineSet.super(row, col, direction);
    		// TODO Auto-generated constructor stub
    	}

        RowColIterator2d(int row, int col, int direction)
        {
        	super(row, col, direction);
        }
        
        public void initialize(int r, int c, int direction)
        {
            int maxRow = lines.length - 1;
        	row = r;
            row = Math.min(row, maxRow);
            row = Math.max(0, row);
            if (lines != null && row < lines.length && row >= 0 && lines[row] != null)
            {
                lineColInc = ((StsPreStackLine2d) lines[row]).analysisColInc;
                lineColStart = (int) (((StsPreStackLine2d) lines[row]).analysisColStart - lines[row].colNumMin); //units of analysisColStart are column index, not sequential column that we need
                lineColStart = (int) StsMath.minMax(lineColStart, 0 , lines[row].nCols);
            }
            if (c < lineColStart)
            {
                c = lineColStart;
            }
            if (direction > 0)
            {
                col = StsMath.intervalRoundDown(c, lineColStart, lineColInc);
            }
            else
            {
                col = StsMath.intervalRoundUp(c, lineColStart, lineColInc);
            }
            endRow = maxRow;
            endRow = StsMath.intervalRoundDown(endRow, rowStart, rowInc);
            endCol = getNColsForRow(row)-1;
            endCol = StsMath.intervalRoundDown(endCol, lineColStart, lineColInc);
            row = StsMath.minMax(this.row, rowStart, endRow);
            col = StsMath.minMax(this.col, lineColStart, endCol);
            startRow = row;
            startCol = col;
            this.direction = direction;
            hasNext = true;
        }

        public Object next()
        {
            if(!hasNext) return null;
            int[] rowCol = new int[] { row, col };
            if(direction > 0)
            {
            	if (col < lineColStart) {
                	col = lineColStart;
                }
                col += lineColInc;
                if(col > endCol)
                {
                	row += rowInc;
                	if(row > endRow)
                        row = rowStart;
                	if (lines != null && lines.length > row && lines[row] != null) {
                		lineColInc = ((StsPreStackLine2d)lines[row]).analysisColInc;
                		lineColStart = (int) (((StsPreStackLine2d)lines[row]).analysisColStart - lines[row].colNumMin);
                	}
                    col = lineColStart;
                    endCol = getNColsForRow(row)-1;
                    endCol = StsMath.intervalRoundDown(endCol, lineColStart, lineColInc);
                }
            }
            else
            {
                col -= lineColInc;
                if(col < lineColStart)
                {
                    row -= rowInc;
                    if(row < rowStart)
                        row = endRow;
                    if (lines != null && lines.length > row && lines[row] != null) {
                		lineColInc = ((StsPreStackLine2d)lines[row]).analysisColInc;
                		lineColStart = (int) (((StsPreStackLine2d)lines[row]).analysisColStart - lines[row].colNumMin);
                	}
                    endCol = getNColsForRow(row)-1;
                    endCol = StsMath.intervalRoundDown(endCol, lineColStart, lineColInc);
                    col = endCol;
                }
            }
            if (col < 0 || col > lines[row].nCols) {
                // System.out.println("StsPreStackLineSet2d.RowColIterator2d: Bad column = " + col );
                col = StsMath.minMax(col, 0, lines[row].nCols);
            }
            hasNext = row != startRow || col != startCol;
            return rowCol;
        }

    }
}
